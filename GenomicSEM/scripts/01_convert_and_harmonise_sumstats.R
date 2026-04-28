#!/usr/bin/env Rscript
# =============================================================================
# Script 01: Convert + harmonise GWAS summary stats to HG19 (GRCh37) with triple outputs
# -----------------------------------------------------------------------------
# VERSION: 2026-01-29-v3
# CHANGES:
#   - Added FUMA output: <TRAIT>.hg19.fuma.sumstats.tsv.gz (minimal schema, all SNPs)
#   - FULL output now uses CHR:BP as SNP identifier (guarantees rsID↔BP consistency)
#   - HM3 output schema reduced to minimal: SNP,CHR,BP,A1,A2,BETA,SE,P,N
#   - Removed N_CASES/N_CONTROLS from HM3 files (FUMA compatibility)
# -----------------------------------------------------------------------------
# PURPOSE
#   - Produce THREE hg19 outputs per trait:
#       (A) FULL  : <TRAIT>.hg19.full.sumstats.tsv.gz   (archival, all columns)
#       (A2) FUMA : <TRAIT>.hg19.fuma.sumstats.tsv.gz   (FUMA upload, minimal schema)
#       (B) HM3   : <TRAIT>.hg19.hm3.sumstats.tsv.gz    (LDSC/GenomicSEM only)
#
# KEY GUARANTEES
#   - Liftover hg38 -> hg19 using UCSC chain:
#       * drops unmapped or multi-mapped variants
#       * FAIL if (unmapped + multi-mapped) > 5%
#   - FULL output uses CHR:BP as SNP identifier:
#       * Guarantees rsID↔BP consistency for FUMA/MAGMA
#       * Avoids silent SNP drops from dbSNP mismatches
#   - HM3 output ensures rsIDs for LDSC by filling missing SNP ids using a
#     HapMap3 (HM3) position->rsID map built from eur_w_ld_chr/{chr}.l2.ldscore.gz
#     NOTE: This rsID fill is acceptable for LDSC alignment, but should NOT be used
#     for fine-grained locus biology. Use FULL for biology.
#
# INPUTS
#   --manifest  config/manifest.tsv
#   --chain     data/reference/liftover/hg38ToHg19.over.chain
#   --ld_ref    data/reference/eur_w_ld_chr   (expects 1..22.l2.ldscore.gz)
#   --outdir    data/processed/hg19
#   --trait     optional, process a single trait
#
# MANIFEST required columns:
#   trait,family,original_file,build,sep,has_header,header_override,
#   col_chr,col_bp,col_snp,col_a1,col_a2,effect_type,col_effect,col_se,
#   p_type,col_p,col_log10p
#
# CANONICAL OUTPUT COLUMNS:
#   FULL: SNP, CHR, BP, A1, A2, BETA, SE, P, N, N_CASES, N_CONTROLS (archival)
#   FUMA: SNP, CHR, BP, A1, A2, BETA, SE, P, N  (minimal schema, all SNPs)
#   HM3:  SNP, CHR, BP, A1, A2, BETA, SE, P, N  (minimal schema, HM3 SNPs only)
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cat("\n[FAIL] Missing package: data.table\n", file=stderr()); quit(status=1)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    cat("\n[FAIL] Missing package: jsonlite\n", file=stderr()); quit(status=1)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE)) {
    cat("\n[FAIL] Missing Bioconductor packages: rtracklayer, GenomicRanges, IRanges\n", file=stderr())
    quit(status=1)
  }
  library(data.table)
  library(jsonlite)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
})

setDTthreads(0L)

die    <- function(msg) { cat("\n[FAIL] ", msg, "\n\n", sep="", file=stderr()); quit(status=1, save="no") }
warn   <- function(msg) cat("[WARN] ", msg, "\n", sep="")
info   <- function(msg) cat("[INFO] ", msg, "\n", sep="")
assert <- function(cond, msg) if (!isTRUE(cond)) die(msg)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i[length(i)] == length(args)) die(paste0("Missing value after ", flag))
  args[i[length(i)] + 1]
}
normalize_safe <- function(p) normalizePath(p, winslash="/", mustWork=FALSE)

# ------------------------------------------------------------------------------
# CLI args
# ------------------------------------------------------------------------------
manifest_path <- normalize_safe(get_arg("--manifest", "config/manifest.tsv"))
chain_path    <- normalize_safe(get_arg("--chain",    "data/reference/liftover/hg38ToHg19.over.chain"))
out_root      <- normalize_safe(get_arg("--outdir",   "data/processed/hg19"))
ld_ref_dir    <- normalize_safe(get_arg("--ld_ref",   "data/reference/eur_w_ld_chr"))
only_trait    <- get_arg("--trait", NULL)

assert(file.exists(manifest_path), paste0("Manifest not found: ", manifest_path))
assert(file.exists(chain_path),    paste0("Chain not found: ", chain_path))
assert(dir.exists(ld_ref_dir),     paste0("LD ref dir not found: ", ld_ref_dir))
dir.create(out_root, recursive=TRUE, showWarnings=FALSE)

# ------------------------------------------------------------------------------
# Load manifest
# ------------------------------------------------------------------------------
manifest <- fread(manifest_path, sep="\t", header=TRUE, quote="", data.table=FALSE,
                  na.strings=c("", "NA"), encoding="UTF-8")

required_cols <- c(
  "trait","family","original_file","build","sep",
  "has_header","header_override",
  "col_chr","col_bp","col_snp","col_a1","col_a2",
  "effect_type","col_effect","col_se",
  "p_type","col_p","col_log10p"
)
missing_req <- setdiff(required_cols, names(manifest))
assert(length(missing_req) == 0, paste0("Manifest missing required columns: ", paste(missing_req, collapse=", ")))

if (!is.null(only_trait)) {
  manifest <- manifest[manifest$trait == only_trait, , drop=FALSE]
  assert(nrow(manifest) == 1, paste0("--trait requested but not found uniquely in manifest: ", only_trait))
}

sep_to_fread <- function(sep) {
  sep <- toupper(trimws(sep))
  if (sep == "TAB")   return("\t")
  if (sep == "COMMA") return(",")
  if (sep == "SPACE") return(" ")
  die(paste0("Unknown sep: ", sep, " (expected TAB/COMMA/SPACE)"))
}

# ------------------------------------------------------------------------------
# Liftover chain (cached)
# ------------------------------------------------------------------------------
info("Loading liftover chain")
CHAIN <- rtracklayer::import.chain(chain_path)

liftover_hg38_to_hg19 <- function(chr_vec, bp_vec) {
  chr_clean <- as.character(chr_vec)
  chr_clean <- gsub("^chr", "", chr_clean, ignore.case=TRUE)
  ok <- chr_clean %in% as.character(1:22) & !is.na(bp_vec) & bp_vec > 0
  keep0 <- which(ok)
  if (length(keep0) == 0) return(list(keep_idx=integer(0), chr=character(0), bp=integer(0), fail_rate=1))
  
  gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", chr_clean[keep0]),
    ranges   = IRanges::IRanges(start=bp_vec[keep0], end=bp_vec[keep0])
  )
  lifted <- rtracklayer::liftOver(gr, CHAIN)
  map_len <- lengths(lifted)
  
  unique_idx <- which(map_len == 1)
  fail_rate <- sum(map_len == 0 | map_len > 1) / length(map_len)
  
  if (length(unique_idx) == 0) return(list(keep_idx=integer(0), chr=character(0), bp=integer(0), fail_rate=fail_rate))
  
  gr2 <- unlist(lifted[unique_idx], use.names=FALSE)
  list(
    keep_idx = keep0[unique_idx],
    chr = gsub("^chr","", as.character(GenomicRanges::seqnames(gr2)), ignore.case=TRUE),
    bp  = as.integer(GenomicRanges::start(gr2)),
    fail_rate = fail_rate
  )
}

# ------------------------------------------------------------------------------
# HapMap3 position -> rsID mapping from {chr}.l2.ldscore.gz
# ------------------------------------------------------------------------------
ldscore_path_for_chr <- function(chr) {
  f <- file.path(ld_ref_dir, paste0(chr, ".l2.ldscore.gz"))
  if (!file.exists(f)) die(paste0("Missing LD-score file for chr ", chr, ": ", f))
  f
}

info("Building HM3 position→rsID map from eur_w_ld_chr/*.l2.ldscore.gz")
HM3_POS <- rbindlist(lapply(1:22, function(chr) {
  f <- ldscore_path_for_chr(chr)
  x <- fread(f, showProgress=FALSE)
  
  nms <- names(x)
  pick <- function(cands) {
    hit <- cands[cands %in% nms]
    if (length(hit)==0) return(NA_character_)
    hit[1]
  }
  c_chr <- pick(c("CHR","Chrom","chrom"))
  c_snp <- pick(c("SNP","snp","RSID","rsid"))
  c_bp  <- pick(c("BP","bp","POS","pos","Position","position"))
  
  assert(!is.na(c_chr) && !is.na(c_snp) && !is.na(c_bp),
         paste0("LD-score file missing needed columns. Found: ", paste(nms, collapse=", "), " in ", f))
  
  out <- x[, .(
    CHR  = suppressWarnings(as.integer(get(c_chr))),
    BP   = suppressWarnings(as.integer(get(c_bp))),
    RSID = as.character(get(c_snp))
  )]
  
  out <- out[CHR == chr & !is.na(BP) & BP > 0 & !is.na(RSID)]
  out <- out[grepl("^rs[0-9]+[a-zA-Z0-9_]*$", RSID, ignore.case=TRUE)]
  out
}))

setorder(HM3_POS, CHR, BP)
HM3_POS <- HM3_POS[, .SD[1], by=.(CHR, BP)]
setkey(HM3_POS, CHR, BP)
info(paste0("HM3 position map: ", format(nrow(HM3_POS), big.mark=","), " SNPs"))

restore_rsids_hm3 <- function(dt, trait) {
  # Ensure proper types without nuking values
  dt[, CHR := suppressWarnings(as.integer(as.character(CHR)))]
  dt[, BP  := suppressWarnings(as.integer(as.character(BP)))]
  
  if (!("SNP" %in% names(dt))) dt[, SNP := NA_character_]
  dt[, SNP := as.character(SNP)]
  
  # IMPORTANT: Fill rsIDs for:
  #   (a) missing SNP (NA, "", "NA")
  #   (b) non-rsID SNPs (e.g., "10:100000625:G:A" - CHR:BP:A1:A2 format from FULL)
  # This ensures we try to replace fallback SNP IDs with real HM3 rsIDs
  is_rs <- !is.na(dt$SNP) & grepl("^rs[0-9]+[a-zA-Z0-9_]*$", dt$SNP, ignore.case=TRUE)
  need <- which(is.na(dt$SNP) | dt$SNP == "" | dt$SNP == "NA" | !is_rs)
  
  if (length(need) > 0) {
    matched <- HM3_POS[dt[need], on=.(CHR, BP), nomatch=NA, .(RSID)]
    dt[need, SNP := matched$RSID]
  }
  
  # Count final rsID-matched variants
  is_rs2 <- !is.na(dt$SNP) & grepl("^rs[0-9]+[a-zA-Z0-9_]*$", dt$SNP, ignore.case=TRUE)
  n_matched <- sum(is_rs2)
  
  info(paste0(trait, ": HM3 rsID available ", format(n_matched, big.mark=","), " / ",
              format(nrow(dt), big.mark=","), " variants (", round(100*n_matched/nrow(dt), 2), "%)"))
  dt
}

# ------------------------------------------------------------------------------
# Per trait processing
# ------------------------------------------------------------------------------
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop=FALSE]
  trait <- as.character(row$trait)
  infile <- as.character(row$original_file)
  build  <- tolower(trimws(as.character(row$build)))
  sep    <- sep_to_fread(as.character(row$sep))
  
  has_header <- as.logical(row$has_header)
  if (is.na(has_header)) has_header <- TRUE
  
  header_override <- as.character(row$header_override)
  if (is.na(header_override) || header_override == "" || header_override == "NA") header_override <- NULL
  
  # Extract N-related columns/constants from manifest
  col_n <- as.character(row$col_n)
  col_n_case <- as.character(row$col_n_case)
  col_n_control <- as.character(row$col_n_control)
  N_const <- as.numeric(row$N_const)
  N_case_const <- as.numeric(row$N_case_const)
  N_control_const <- as.numeric(row$N_control_const)
  
  # Clean up empty/NA strings and NAs
  if (is.na(col_n) || col_n == "" || col_n == "NA") col_n <- NA_character_
  if (is.na(col_n_case) || col_n_case == "" || col_n_case == "NA") col_n_case <- NA_character_
  if (is.na(col_n_control) || col_n_control == "" || col_n_control == "NA") col_n_control <- NA_character_
  
  info(paste0("=== TRAIT: ", trait, " | build=", build, " ==="))
  assert(file.exists(infile), paste0(trait, ": input file not found: ", infile))
  assert(build %in% c("hg19","hg38"), paste0(trait, ": build must be hg19 or hg38"))
  
  outdir <- file.path(out_root, trait)
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  
  out_sumstats_full <- file.path(outdir, paste0(trait, ".hg19.full.sumstats.tsv.gz"))
  out_sumstats_hm3  <- file.path(outdir, paste0(trait, ".hg19.hm3.sumstats.tsv.gz"))
  out_report        <- file.path(outdir, paste0(trait, ".conversion_report.json"))
  done_flag         <- file.path(outdir, ".done")
  
  if (file.exists(done_flag) && file.exists(out_sumstats_full) && file.exists(out_sumstats_hm3) && file.exists(out_report)) {
    info(paste0(trait, ": already done, skipping"))
    next
  }
  
  # Read file
  if (!is.null(header_override)) {
    # Remove surrounding quotes that may have been added during TSV parsing
    header_override <- gsub('^"|"$', '', header_override)
    override_cols <- strsplit(header_override, "[[:space:],]+")[[1]]
    override_cols <- override_cols[override_cols != ""]
    info(paste0(trait, ": Header override (", length(override_cols), " cols)"))
    dt <- fread(infile, sep=sep, header=FALSE, showProgress=FALSE, skip=ifelse(has_header, 1, 0))
    if (ncol(dt) >= length(override_cols)) {
      setnames(dt, seq_along(override_cols), override_cols)
    } else {
      die(paste0(trait, ": header_override has ", length(override_cols), " columns but file has only ", ncol(dt)))
    }
  } else if (has_header) {
    dt <- fread(infile, sep=sep, header=TRUE, showProgress=FALSE)
  } else {
    die(paste0(trait, ": has_header=FALSE but no header_override provided"))
  }
  
  info(paste0(trait, ": Read ", format(nrow(dt), big.mark=","), " rows"))
  
  # ==========================================================================
  # FIX: Strip # prefix from first column name (common REGENIE/PLINK artifact)
  # Some tools output files where the first column has a leading # (e.g., "#chrom")
  # This causes column matching to fail when manifest specifies "chrom"
  # ==========================================================================
  first_col <- names(dt)[1]
  if (startsWith(first_col, "#")) {
    new_name <- sub("^#", "", first_col)
    setnames(dt, first_col, new_name)
    info(paste0(trait, ": Renamed column '", first_col, "' → '", new_name, "'"))
  }
  
  col_chr <- sub("^#", "", as.character(row$col_chr))
  col_bp  <- sub("^#", "", as.character(row$col_bp))
  col_a1  <- sub("^#", "", as.character(row$col_a1))
  col_a2  <- sub("^#", "", as.character(row$col_a2))
  col_snp <- sub("^#", "", as.character(row$col_snp))
  col_eff <- sub("^#", "", as.character(row$col_effect))
  col_se  <- sub("^#", "", as.character(row$col_se))
  
  p_type  <- toupper(trimws(as.character(row$p_type)))
  col_p   <- sub("^#", "", as.character(row$col_p))
  col_l10 <- sub("^#", "", as.character(row$col_log10p))
  
  effect_type <- toupper(trimws(as.character(row$effect_type)))
  
  # ==========================================================================
  # CRITICAL FIX: Force BP column to character IMMEDIATELY to preserve
  # precision. When fread reads large genomic positions (9+ digits) as
  # numeric/double, it loses leading digits due to floating-point precision.
  # By converting to character immediately, we lock in the exact values
  # before any numeric conversion.
  # ==========================================================================
  if (col_bp %in% names(dt)) {
    dt[[col_bp]] <- as.character(dt[[col_bp]])
    info(paste0(trait, ": BP column locked as character (precision preservation)"))
  }
  
  # Verify required columns
  required_input_cols <- c(col_chr, col_bp, col_a1, col_a2, col_eff, col_se)
  if (p_type == "LOG10P") {
    required_input_cols <- c(required_input_cols, col_l10)
  } else {
    required_input_cols <- c(required_input_cols, col_p)
  }
  missing_cols <- setdiff(required_input_cols, names(dt))
  if (length(missing_cols) > 0) {
    die(paste0(trait, ": Missing columns: ", paste(missing_cols, collapse=", "),
               "\nAvailable: ", paste(names(dt), collapse=", ")))
  }
  
  CHR <- as.character(dt[[col_chr]])
  CHR <- gsub("^chr","", CHR, ignore.case=TRUE)
  CHR <- suppressWarnings(as.integer(CHR))
  
  BP  <- suppressWarnings(as.integer(dt[[col_bp]]))
  A1  <- toupper(as.character(dt[[col_a1]]))
  A2  <- toupper(as.character(dt[[col_a2]]))
  
  SNP_RAW <- NA_character_
  if (!is.na(col_snp) && nzchar(col_snp) && (col_snp %in% names(dt))) {
    SNP_RAW <- as.character(dt[[col_snp]])
  }
  
  EFF <- suppressWarnings(as.numeric(dt[[col_eff]]))
  SE  <- suppressWarnings(as.numeric(dt[[col_se]]))
  
  if (p_type == "LOG10P") {
    LOG10P <- suppressWarnings(as.numeric(dt[[col_l10]]))
    P <- 10^(-LOG10P)
  } else {
    P <- suppressWarnings(as.numeric(dt[[col_p]]))
  }
  
  # Extract N / N_CASES / N_CONTROLS
  # Priority: column from file > constant from manifest > NA
  N <- NA_integer_
  N_CASES <- NA_integer_
  N_CONTROLS <- NA_integer_
  
  # Try to get total N
  if (!is.na(col_n) && nzchar(col_n) && (col_n %in% names(dt))) {
    N <- suppressWarnings(as.integer(dt[[col_n]]))
    N <- ifelse(is.na(N), NA_integer_, N)
  } else if (!is.na(N_const)) {
    N <- as.integer(N_const)
  }
  
  # Try to get case/control N
  if (!is.na(col_n_case) && nzchar(col_n_case) && (col_n_case %in% names(dt))) {
    N_CASES <- suppressWarnings(as.integer(dt[[col_n_case]]))
    N_CASES <- ifelse(is.na(N_CASES), NA_integer_, N_CASES)
  } else if (!is.na(N_case_const)) {
    N_CASES <- as.integer(N_case_const)
  }
  
  if (!is.na(col_n_control) && nzchar(col_n_control) && (col_n_control %in% names(dt))) {
    N_CONTROLS <- suppressWarnings(as.integer(dt[[col_n_control]]))
    N_CONTROLS <- ifelse(is.na(N_CONTROLS), NA_integer_, N_CONTROLS)
  } else if (!is.na(N_control_const)) {
    N_CONTROLS <- as.integer(N_control_const)
  }
  
  if (effect_type == "OR") {
    bad_or <- which(is.na(EFF) | !is.finite(EFF) | EFF <= 0)
    if (length(bad_or) > 0) warn(paste0(trait, ": ", length(bad_or), " invalid OR values will be filtered"))
    BETA <- log(EFF)
  } else {
    BETA <- EFF
  }
  
  canon <- data.table(SNP_RAW=SNP_RAW, CHR=CHR, BP=BP, A1=A1, A2=A2, BETA=BETA, SE=SE, P=P, N=N, N_CASES=N_CASES, N_CONTROLS=N_CONTROLS)
  n_in <- nrow(canon)
  
  # Basic QC filters (SNPs only; drop indels for LDSC stability)
  canon <- canon[
    CHR %in% 1:22 &
      !is.na(BP) & BP > 0 &
      !is.na(A1) & !is.na(A2) &
      nchar(A1)==1L & nchar(A2)==1L &
      A1 %chin% c("A","C","G","T") &
      A2 %chin% c("A","C","G","T") &
      A1 != A2 &
      is.finite(BETA) &
      is.finite(SE) & SE > 0 &
      is.finite(P) & P > 0 & P <= 1
  ]
  
  info(paste0(trait, ": After QC: ", format(nrow(canon), big.mark=","), " variants"))
  assert(nrow(canon) > 0, paste0(trait, ": no rows after QC filtering"))
  
  report <- list(
    trait = trait,
    input_file = infile,
    build_input = build,
    build_output = "hg19",
    n_input_rows = n_in,
    n_after_qc = nrow(canon)
  )
  
  # Liftover if hg38
  if (build == "hg38") {
    info(paste0(trait, ": Liftover hg38 → hg19"))
    lo <- liftover_hg38_to_hg19(canon$CHR, canon$BP)
    report$liftover_fail_rate <- lo$fail_rate
    
    info(paste0(trait, ": Liftover kept ", format(length(lo$keep_idx), big.mark=","),
                " variants (fail rate: ", round(lo$fail_rate*100, 2), "%)"))
    
    if (lo$fail_rate > 0.05) die(paste0(trait, ": liftover failure rate >5%"))
    
    canon <- canon[lo$keep_idx]
    canon[, CHR := as.integer(lo$chr)]
    canon[, BP  := as.integer(lo$bp)]
    
    # Verify lifted CHR values are still valid (1-22)
    bad_chr <- which(!(canon$CHR %in% 1:22))
    if (length(bad_chr) > 0) {
      warn(paste0(trait, ": ", length(bad_chr), " variants lifted to invalid CHR; removing"))
      canon <- canon[canon$CHR %in% 1:22]
    }
  } else {
    report$liftover_fail_rate <- NA_real_
  }
  
  # Build-consistent hg19 variant key (for FULL fallback when rsID absent)
  canon[, VARID19 := paste0(CHR, ":", BP, ":", A1, ":", A2)]
  
  # ---------------------------------------------------------------------------
  # Output A (FULL): MAGMA/FUMA-safe file
  # - Use CHR:BP as SNP identifier (guarantees rsID↔BP consistency for FUMA)
  # - Do NOT HM3-filter
  # - Retains N_CASES/N_CONTROLS for archival purposes
  # ---------------------------------------------------------------------------
  full <- copy(canon)
  
  # Force CHR:BP for all SNPs - guarantees position consistency for FUMA/MAGMA
  full[, SNP := paste0(CHR, ":", BP)]
  
  setorder(full, SNP, P)
  full <- full[, .SD[1], by=SNP]
  
  fwrite(full[,.(SNP,CHR,BP,A1,A2,BETA,SE,P,N,N_CASES,N_CONTROLS)],
         out_sumstats_full, sep="\t", quote=FALSE, compress="gzip")
  
  report$n_full <- nrow(full)
  info(paste0(trait, ": Wrote FULL hg19 file (archival): ", format(nrow(full), big.mark=","), " variants"))
  
  # ---------------------------------------------------------------------------
  # Output A2 (FUMA): Minimal schema for FUMA upload
  # - Same SNPs as FULL (not HM3-filtered)
  # - Minimal columns: SNP, CHR, BP, A1, A2, BETA, SE, P, N
  # ---------------------------------------------------------------------------
  out_sumstats_fuma <- file.path(outdir, paste0(trait, ".hg19.fuma.sumstats.tsv.gz"))
  fwrite(full[,.(SNP,CHR,BP,A1,A2,BETA,SE,P,N)],
         out_sumstats_fuma, sep="\t", quote=FALSE, compress="gzip")
  info(paste0(trait, ": Wrote FUMA hg19 file (minimal schema): ", format(nrow(full), big.mark=","), " variants"))
  
  # ---------------------------------------------------------------------------
  # Output B (HM3): LDSC/GenomicSEM-safe file
  # - HARD HM3 SUBSET: Keep ONLY variants that exist in HapMap3 reference
  # - This ensures LDSC/GenomicSEM assumptions are met (~1.1-1.3M variants per trait)
  # ---------------------------------------------------------------------------
  hm3 <- copy(canon)
  
  # CRITICAL: Join with HM3_POS to subset to HM3 membership
  # This is a right join: keep only variants in BOTH hm3 AND HM3_POS
  # nomatch=0 means drop variants not found in HM3_POS
  info(paste0(trait, ": Subsetting to HM3 reference (", format(nrow(HM3_POS), big.mark=","), " SNPs)"))
  hm3 <- HM3_POS[hm3, on=.(CHR,BP), nomatch=0]
  
  # Force SNP to HM3 rsID (guaranteed valid HM3 rsID from the join)
  hm3[, SNP := RSID]
  
  report$n_hm3_matched <- nrow(hm3)
  
  if (nrow(hm3) == 0) die(paste0(trait, ": 0 HM3 variants matched - file coordinates don't overlap with HM3 reference (check build)"))
  
  setorder(hm3, SNP, P)
  hm3 <- hm3[, .SD[1], by=SNP]
  
  fwrite(hm3[,.(SNP,CHR,BP,A1,A2,BETA,SE,P,N)],
         out_sumstats_hm3, sep="\t", quote=FALSE, compress="gzip")
  
  # By-chromosome splits (HM3 file)
  by_chr_dir <- file.path(outdir, paste0(trait, ".hg19.hm3.by_chr"))
  dir.create(by_chr_dir, recursive=TRUE, showWarnings=FALSE)
  for (chr in 1:22) {
    chr_dt <- hm3[CHR == chr, .(SNP,CHR,BP,A1,A2,BETA,SE,P,N)]
    fwrite(chr_dt, file.path(by_chr_dir, paste0("chr", chr, ".tsv.gz")),
           sep="\t", quote=FALSE, compress="gzip")
  }
  info(paste0(trait, ": Wrote HM3 by_chr files"))
  
  report$n_final_hm3 <- nrow(hm3)
  
  writeLines(toJSON(report, pretty=TRUE, auto_unbox=TRUE, null="null"), out_report)
  writeLines("OK", done_flag)
  
  info(paste0(trait, ": DONE ✅  FULL=", format(report$n_full, big.mark=","), " ; HM3=", format(report$n_final_hm3, big.mark=",")))
}

info("All traits processed ✅")
quit(status=0, save="no")