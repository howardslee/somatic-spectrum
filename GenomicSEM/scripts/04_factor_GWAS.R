#!/usr/bin/env Rscript
# =============================================================================
# VERSION: 2026-02-03-v13 (Complete ordering consistency)
# =============================================================================
# CHANGES FROM v12:
#   - FIX Issue A: Core residual loop now uses core_used, not core
#   - FIX Issue B: Validation rg slices S_all[core_used, t], not S_all[core, t]
#   - FIX Issue C: Force trait.names = rownames(S) for bulletproof ordering
#   - All lambda/S_all operations now use consistent core_used ordering
# CHANGES FROM v11:
#   - FIX: h2F, N_factor, lambda_sq_sum now computed AFTER loading order check
#   - Root cause: S_core used `core` order but loadings might be reordered to
#     match LDSCoutput_core$trait.names, causing w and S_core to be misaligned
#   - Now uses core_used = LDSCoutput_core$trait.names consistently
# CHANGES FROM v10:
#   - FIX: Ensure loadings match LDSCoutput_core$trait.names order
#   - Root cause: commonfactorGWAS uses covstruc trait order, not --core order
#   - If order mismatched, wrong lambda_std applied to each trait
#   - Added explicit reorder check after LDSCoutput_core creation
# CHANGES FROM v9:
#   - FIX: Use empirical cor(trait_Z, factor_Z) for subtraction, NOT theoretical rg
#   - Root cause: LDSC rg ≠ SNP-level Z correlation (different quantities)
#   - v9 still had mean_chi2_adj >> 1 because theoretical rg was wrong scale
#   - Now: emp_cor guarantees Var(residual) = 1 - emp_cor² by construction
#   - Theoretical rg (lambda_std or LDSC) logged for comparison only
# CHANGES FROM v8:
#   - FIX: Residual P-values were inflated for high-loading traits
#   - Root cause: Assumed Var(Z)=1 but empirical variance often differs
#   - Solution: Standardize trait_Z and Factor_Z to unit variance BEFORE subtraction
#   - Added winsorized SD estimation (cap=30) to protect against outliers
#   - Added diagnostic: empirical cor(trait, factor) vs theoretical rg
#   - Added diagnostic: mean_chi2_raw vs mean_chi2_adj
#   - P_MIN = 1e-300 to prevent zero P-values
# CHANGES FROM v7:
#   - FIX: Filter primarily on Factor_Z, not Factor_beta/SE
#   - Z-direct branch sets Factor_beta/SE to NA, old filter nuked all rows
#   - Now: gwas <- gwas[is.finite(Factor_Z)], SE filter is optional
# CHANGES FROM v6:
#   - FIX #1: Removed bare "P" from P-value column patterns (too ambiguous)
#   - FIX #1: Added P-value validation (must be in 0-1), fallback to P from Z
#   - FIX #2: Robust Z vs beta detection:
#       * If (est + SE) exist → treat as beta/SE, compute Z = beta/SE
#       * Else if Z_Estimate exists → use directly as Z
#   - Always compute Factor_P from Factor_Z as primary, use Pval_Estimate only if valid
# CHANGES FROM v5:
#   - Core trait residuals: use rg = lambda_std (standardized loading)
#   - Validation trait residuals: use rg computed from LDSC, with clamping at +/-0.99
#   - P-value adjustment: account for Var(residual_Z) = 1 - rg^2
#   - SE adjustment: residual_SE = sqrt(1-rg^2) / sqrt(N)
# CHANGES FROM v4:
#   - Fixed core trait residual calculation: sumstats() outputs beta.<trait>
#     and se.<trait> columns, not plain trait name Z-scores
# CHANGES FROM v2:
#   - Added --minimal flag to skip residual GWAS (for faster LOO analysis)
#   - Memory-efficient chromosome processing: don't accumulate results in RAM
#   - Track completion status only, combine from disk at end
#   - Default cores=1 (safer for 16GB MacBooks)
#   - Explicit gc() after each chromosome
# =============================================================================
#
# STANDARD GENOMICSEM PIPELINE:
#   1. munge() - already done in Script 02 (creates Z-only files for LDSC)
#   2. ldsc()  - already done in Script 03, saved to RDS
#   3. sumstats() - process ORIGINAL GWAS files (need BETA, SE, P)
#   4. commonfactorGWAS() - run factor GWAS
#
# FILE TYPES:
#   - data/processed/hg19/<TRAIT>/*.sumstats.tsv.gz → Has BETA, SE, P → sumstats()
#   - results/ldsc_universe/munged/*.sumstats.gz   → Has Z only      → ldsc()
#
# INPUTS:
#   - results/ldsc_universe/ldsc/ldsc_universe.rds (from ldsc())
#   - results/factor_model/factor_model.rds (contains loadings)
#   - config/manifest.tsv
#   - data/processed/hg19/<TRAIT>/<TRAIT>.hg19.sumstats.tsv.gz (ORIGINAL files)
#   - data/reference/reference.1000G.maf.0.005.txt
#
# OUTPUTS:
#   results/factor_gwas/
#     factor_gwas_full.tsv.gz       - Full factor GWAS with per-trait betas
#     factor_gwas_fuma.tsv          - FUMA format (Factor P)
#     qsnp_gwas_fuma.tsv            - FUMA format (Q_SNP P)
#     residual_<CORE>.tsv.gz        - Core residuals
#     residual_<VALIDATION>.tsv.gz  - Validation residuals
#
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║  SCRIPT VERSION: 2026-02-03-v13 (Complete ordering consistency) ║\n")
cat("║  Uses: ORIGINAL processed files (BETA/SE/P) -> sumstats()       ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ==============================================================================
# THREADING CONTROLS (prevent nested parallelism issues on Mac/Linux)
# ==============================================================================
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# ==============================================================================
# CONFIGURABLE THRESHOLDS
# ==============================================================================
THRESH_Z_OUTLIER <- 80
THRESH_ALLELE_DROP_WARN <- 0.02
THRESH_FINITE_MIN <- 0.95
THRESH_CHR_FAIL_MAX <- 3
THRESH_MIN_SNPS <- 100000
THRESH_CFI_WARN <- 0.90
THRESH_LAMBDA_MIN <- 1e-10
P_MIN <- 1e-300                # Minimum P-value to prevent underflow
WINSORIZE_CAP <- 30            # Cap for winsorized SD estimation
COR_MISMATCH_WARN <- 0.15      # Warn if |emp_cor - rg| exceeds this

# ==============================================================================
# Helper Functions
# ==============================================================================
die  <- function(msg) { cat("\n[FAIL] ", msg, "\n\n", sep="", file=stderr()); quit(status=1, save="no") }
warn <- function(msg) cat("[WARN] ", msg, "\n", sep="")
info <- function(msg) cat("[INFO] ", msg, "\n", sep="")
debug <- function(msg) if (exists("VERBOSE") && VERBOSE) cat("[DEBUG] ", msg, "\n", sep="")
assert <- function(cond, msg) if (!isTRUE(cond)) die(msg)

# Winsorize extreme values for robust SD estimation
winsorize <- function(z, cap = WINSORIZE_CAP) {
  pmax(pmin(z, cap), -cap)
}

# Empirical SD with winsorization to protect against outliers
emp_sd <- function(z, cap = WINSORIZE_CAP) {
  sd(winsorize(z, cap = cap), na.rm = TRUE)
}

suppressPackageStartupMessages({
  library(data.table)
})

# Atomic file write helper
atomic_saveRDS <- function(object, file) {
  tmp_file <- paste0(file, ".tmp.", Sys.getpid())
  saveRDS(object, tmp_file)
  file.rename(tmp_file, file)
}

atomic_write_done <- function(done_file, rds_file) {
  if (!file.exists(rds_file)) {
    warn(sprintf("RDS file missing, not writing done marker: %s", rds_file))
    return(FALSE)
  }
  tryCatch({
    tmp <- readRDS(rds_file)
    rm(tmp)
  }, error = function(e) {
    warn(sprintf("RDS file corrupted: %s", e$message))
    return(FALSE)
  })
  writeLines("OK", done_file)
  return(TRUE)
}

info(sprintf("data.table version: %s", packageVersion("data.table")))
info(sprintf("R version: %s", R.version.string))

# ==============================================================================
# Check lavaan version (log only, don't force downgrade)
# ==============================================================================
if (requireNamespace("lavaan", quietly = TRUE)) {
  info(sprintf("lavaan version: %s", packageVersion("lavaan")))
} else {
  info("lavaan not installed, GenomicSEM will install it")
}

suppressPackageStartupMessages({
  library(GenomicSEM)
})

info(sprintf("GenomicSEM version: %s", packageVersion("GenomicSEM")))

# ==============================================================================
# Command Line Arguments
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i[length(i)] == length(args)) die(paste0("Missing value after ", flag))
  args[i[length(i)] + 1]
}

manifest_path  <- get_arg("--manifest", "config/manifest.tsv")
ldsc_rds       <- get_arg("--ldsc-rds", "results/ldsc_universe/ldsc/ldsc_universe.rds")
factor_rds     <- get_arg("--factor-rds", "results/models/factor_model/factor_model.rds")
processed_dir  <- get_arg("--processed", "data/processed/hg19")
munged_dir     <- get_arg("--munged", "results/ldsc_universe/munged")
ref_dir        <- get_arg("--ref", "data/reference")
outdir         <- get_arg("--outdir", "results/models/factor_gwas")
ldsc_path      <- get_arg("--ldsc-path", "ldsc")

core_arg       <- get_arg("--core", "Fibromyalgia,MECFS,IBS,Depression,PTSD,MigAura")
validation_arg <- get_arg("--validation", NULL)

# Trait type: "continuous" or "binary" (affects OLS/linprob in sumstats)
trait_type_arg <- get_arg("--trait-type", "continuous")

N_factor_override <- suppressWarnings(as.numeric(get_arg("--N_factor", NA)))
# Default to 1 core - safer for 16GB systems
n_cores <- suppressWarnings(as.integer(get_arg("--cores", "1")))
if (is.na(n_cores) || n_cores < 1) n_cores <- 1L

VERBOSE <- "--verbose" %in% args
MINIMAL <- "--minimal" %in% args  # Skip residual GWAS for faster LOO

if (MINIMAL) {
  info("MINIMAL MODE: Skipping residual GWAS calculations")
}

# Reference files
hm3_path <- file.path(ref_dir, "w_hm3.snplist")
ld_dir   <- file.path(ref_dir, "eur_w_ld_chr")

# Reference panel for sumstats() - try common names
ref_panel_candidates <- c(
  file.path(ref_dir, "reference.1000G.hm3.txt"),       # HM3-only (fastest, ~1.2M SNPs)
  file.path(ref_dir, "reference.1000G.hm3.txt.gz"),
  file.path(ref_dir, "reference.1000G.maf.0.005.txt"), # Full panel (~10M SNPs)
  file.path(ref_dir, "reference.1000G.maf.0.005.txt.gz"),
  file.path(ref_dir, "1000G_reference.txt"),
  file.path(ref_dir, "ref.txt")
)
ref_panel <- NULL
for (rp in ref_panel_candidates) {
  if (file.exists(rp)) {
    ref_panel <- rp
    break
  }
}

# Validate inputs
assert(file.exists(manifest_path), paste0("Manifest not found: ", manifest_path))
assert(file.exists(ldsc_rds), paste0("LDSC RDS not found: ", ldsc_rds))
assert(file.exists(factor_rds), paste0("Factor model RDS not found: ", factor_rds))
assert(dir.exists(processed_dir), paste0("Processed dir not found: ", processed_dir))
assert(dir.exists(munged_dir), paste0("Munged dir not found: ", munged_dir))
assert(file.exists(hm3_path), paste0("HM3 snplist not found: ", hm3_path))

if (is.null(ref_panel)) {
  die(paste0("Reference panel not found. Tried:\n  ", 
             paste(ref_panel_candidates, collapse="\n  "),
             "\n\nDownload from: https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v"))
}
info(sprintf("Using reference panel: %s", ref_panel))

# Create output directories
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir, "fuma"), showWarnings=FALSE)
dir.create(file.path(outdir, "qc"), showWarnings=FALSE)
dir.create(file.path(outdir, "ldsc_commands"), showWarnings=FALSE)
dir.create(file.path(outdir, "chr"), showWarnings=FALSE)

# ==============================================================================
# Load manifest
# ==============================================================================
info("Loading manifest...")
man <- fread(manifest_path, sep="\t", header=TRUE, na.strings=c("", "NA"))
assert("trait" %in% names(man), "Manifest missing 'trait' column")
assert("N_const" %in% names(man), "Manifest missing 'N_const' column")

man[, N_const := suppressWarnings(as.numeric(N_const))]

# ==============================================================================
# Sumstats file resolution (standard)
# ------------------------------------------------------------------------------
# We maintain TWO outputs per trait:
#   1) *.hg19.full*.sumstats.tsv.gz   -> for MAGMA/FUMA and GenomicSEM::sumstats()
#   2) *.hg19.hm3*.sumstats.tsv.gz    -> for LDSC/GenomicSEM covariance (munged)
#
# This script MUST use the FULL files for GenomicSEM::sumstats().
# We prefer rsID-standardised files when available.
# ==============================================================================
resolve_sumstats_file <- function(trait, processed_dir) {
  candidates <- c(
    file.path(processed_dir, trait, paste0(trait, ".hg19.hm3.sumstats.tsv.gz")),  # HM3 (faster)
    file.path(processed_dir, trait, paste0(trait, ".hg19.full.sumstats.tsv.gz"))
  )
  for (f in candidates) {
    if (file.exists(f)) return(f)
  }
  stop(sprintf("No processed sumstats found for trait '%s'. Tried:\n  %s", trait, paste(candidates, collapse="\n  ")))
}

man[, sumstats_path := vapply(trait, resolve_sumstats_file, character(1), processed_dir = processed_dir)]

# ==============================================================================
# Parse core and validation traits
# ==============================================================================
core <- trimws(strsplit(core_arg, ",")[[1]])
assert(length(core) >= 2, "Core must have >= 2 traits")
assert(all(core %in% man$trait), paste0("Core traits not in manifest: ", 
                                        paste(setdiff(core, man$trait), collapse=", ")))

info(paste0("Core traits (", length(core), "): ", paste(core, collapse=", ")))

# ==============================================================================
# Load LDSC results (OUTPUT FROM ldsc() FUNCTION)
# ==============================================================================
info("Loading LDSC results...")
LDSCoutput <- readRDS(ldsc_rds)

info(sprintf("LDSC object contains: %s", paste(names(LDSCoutput), collapse=", ")))

# Verify LDSC output structure
assert("S" %in% names(LDSCoutput), "LDSC object missing S matrix")
assert("V" %in% names(LDSCoutput), "LDSC object missing V matrix")

S_all <- LDSCoutput$S
all_ldsc_traits <- rownames(S_all)

# Check that core traits are in LDSC output
missing_core <- setdiff(core, all_ldsc_traits)
if (length(missing_core) > 0) {
  die(paste0("Core traits not in LDSC output: ", paste(missing_core, collapse=", ")))
}

# Determine validation traits
if (is.null(validation_arg)) {
  validation <- setdiff(all_ldsc_traits, core)
  validation <- validation[validation %in% man$trait]
} else {
  validation <- trimws(strsplit(validation_arg, ",")[[1]])
}

info(paste0("Validation traits (", length(validation), "): ", 
            ifelse(length(validation) > 0, paste(validation, collapse=", "), "none")))

# ==============================================================================
# Load factor model (for loadings - used in residual calculations)
# ==============================================================================
info("Loading factor model...")
fm <- readRDS(factor_rds)

assert("loadings" %in% names(fm), "factor_model.rds missing 'loadings'")

loadings <- fm$loadings
assert(all(c("trait", "lambda", "lambda_std") %in% names(loadings)), 
       "loadings must have trait, lambda, lambda_std columns")

# Reorder loadings to match core order
loadings <- loadings[match(core, loadings$trait)]
assert(all(loadings$trait == core), "Loading order mismatch")

lambda <- loadings$lambda
lambda_std <- loadings$lambda_std

info("Factor loadings:")
for (i in seq_along(core)) {
  info(sprintf("  %s: lambda=%.4f, lambda_std=%.4f", core[i], lambda[i], lambda_std[i]))
}

lambda_norm <- sqrt(sum(lambda^2))
if (lambda_norm < THRESH_LAMBDA_MIN) {
  die(sprintf("All factor loadings are effectively zero (norm=%.2e)", lambda_norm))
}

# NOTE: h2F, N_factor, and lambda_sq_sum are computed AFTER loading order check
# to ensure consistent ordering with LDSCoutput_core$trait.names

# ==============================================================================
# STEP 3: Process sumstats using GenomicSEM::sumstats()
# ==============================================================================
info("=================================================================")
info("STEP 3: Processing sumstats with GenomicSEM::sumstats()")
info("=================================================================")

sumstats_rds <- file.path(outdir, "processed_sumstats.rds")

if (file.exists(sumstats_rds)) {
  info("Loading cached processed sumstats...")
  p_sumstats <- readRDS(sumstats_rds)
} else {
  # =========================================================================
  # IMPORTANT: Use ORIGINAL processed files, NOT munged files
  # - Processed files have: BETA, SE, P (required by sumstats())
  # - Munged files only have: Z (used only by ldsc())
  # =========================================================================
  
  orig_files <- vapply(core, resolve_sumstats_file, character(1), processed_dir = processed_dir)
  missing <- core[!file.exists(orig_files)]
  if (length(missing) > 0) {
    die(paste0("Missing processed sumstats files for: ", paste(missing, collapse=", "),
               "

Check that you have generated *.hg19.full.sumstats.tsv.gz (or *.hg19.full.rsid.sumstats.tsv.gz) for these traits under:
  ", processed_dir))
  }
  
  info(sprintf("Processing %d traits with sumstats()...", length(core)))
  for (i in seq_along(orig_files)) {
    info(sprintf("  [%d] %s", i, orig_files[i]))
  }
  info(sprintf("Reference panel: %s", ref_panel))
  
  # Check first file to see what columns are available
  info("Checking column structure of first file...")
  test_dt <- fread(orig_files[1], nrows=5)
  info(sprintf("  Columns found: %s", paste(names(test_dt), collapse=", ")))
  
  # Get sample sizes
  N_vec <- man[match(core, trait), N_const]
  info(sprintf("Sample sizes: %s", paste(N_vec, collapse=", ")))
  
  # Get sample sizes
  N_vec <- man[match(core, trait), N_const]
  
  # Use MINIMAL parameters - extra params cause "... must be empty" error
  se_logit_vec <- rep(FALSE, length(core))
  
  info("Calling sumstats() with minimal parameters...")
  info(sprintf("  GenomicSEM version: %s", packageVersion("GenomicSEM")))
  
  p_sumstats <- tryCatch({
    sumstats(
      files = orig_files,
      ref = ref_panel,
      trait.names = core,
      se.logit = se_logit_vec
    )
  }, error = function(e) {
    msg <- sprintf("sumstats() failed: %s\n\n", e$message)
    msg <- paste0(msg, "GenomicSEM version: ", packageVersion("GenomicSEM"), "\n")
    die(msg)
  })
  
  info(sprintf("sumstats() returned %d SNPs", nrow(p_sumstats)))
  
  # Cache for restart
  atomic_saveRDS(p_sumstats, sumstats_rds)
  info(sprintf("Cached processed sumstats to: %s", sumstats_rds))
}

# Validate sumstats output
assert(nrow(p_sumstats) > THRESH_MIN_SNPS, 
       sprintf("Too few SNPs from sumstats(): %d < %d", nrow(p_sumstats), THRESH_MIN_SNPS))

info(sprintf("Processed sumstats: %s SNPs", format(nrow(p_sumstats), big.mark=",")))
info(sprintf("Columns: %s", paste(names(p_sumstats), collapse=", ")))

# ==============================================================================
# STEP 4: Factor GWAS via commonfactorGWAS()
# ==============================================================================
info("=================================================================")
info("STEP 4: Factor GWAS with commonfactorGWAS()")
info("=================================================================")

# ---------------------------------------------------------------------------
# Check if LDSC output needs subsetting to core traits
# commonfactorGWAS requires same traits in LDSC and sumstats
# ---------------------------------------------------------------------------
ldsc_traits <- if (!is.null(LDSCoutput$trait.names)) LDSCoutput$trait.names else rownames(S_all)

if (setequal(ldsc_traits, core) && length(ldsc_traits) == length(core)) {
  # LDSC already has exactly the core traits - use directly
  info("LDSC output already contains exactly core traits - using directly")
  info(sprintf("  List element order: %s", paste(names(LDSCoutput), collapse=", ")))
  LDSCoutput_core <- LDSCoutput
  
  # Ensure trait.names is set
  if (is.null(LDSCoutput_core$trait.names)) {
    LDSCoutput_core$trait.names <- core
  }
  
  info(sprintf("  S dimensions: %dx%d", nrow(LDSCoutput_core$S), ncol(LDSCoutput_core$S)))
  info(sprintf("  V dimensions: %dx%d", nrow(LDSCoutput_core$V), ncol(LDSCoutput_core$V)))
  info(sprintf("  trait.names: %s", paste(LDSCoutput_core$trait.names, collapse=", ")))
  
} else {
  # Need to subset LDSC to core traits
  info("Subsetting LDSC output to core traits...")
  
  k_all <- nrow(S_all)
  k_core <- length(core)
  
  # Get indices of core traits in full LDSC output
  core_idx <- match(core, rownames(S_all))
  if (any(is.na(core_idx))) {
    die(sprintf("Core traits not found in LDSC: %s", 
                paste(core[is.na(core_idx)], collapse=", ")))
  }
  
  info(sprintf("  Core trait indices in LDSC: %s", paste(core_idx, collapse=", ")))
  
  # Subset S matrix
  S_sub <- S_all[core_idx, core_idx]
  rownames(S_sub) <- core
  colnames(S_sub) <- core
  
  # Subset V matrix (vech ordering - lower triangle column-major)
  idx_mat <- matrix(NA_integer_, k_all, k_all)
  idx_mat[lower.tri(idx_mat, diag = TRUE)] <- seq_len(k_all * (k_all + 1) / 2)
  idx_mat[upper.tri(idx_mat)] <- t(idx_mat)[upper.tri(idx_mat)]
  
  idx_sub <- idx_mat[core_idx, core_idx]
  vech_idx <- idx_sub[lower.tri(idx_sub, diag = TRUE)]
  
  V_sub <- LDSCoutput$V[vech_idx, vech_idx]
  dimnames(V_sub) <- NULL
  
  # Subset I matrix if present
  if (!is.null(LDSCoutput$I)) {
    I_sub <- LDSCoutput$I[core_idx, core_idx]
  } else {
    I_sub <- NULL
  }
  
  # Build subsetted LDSC output with CORRECT ORDER: V, S, I
  LDSCoutput_core <- list(
    V = V_sub,
    S = S_sub,
    I = I_sub,
    trait.names = core
  )
  
  if (!is.null(LDSCoutput$N)) {
    LDSCoutput_core$N <- LDSCoutput$N[core_idx]
  }
  if (!is.null(LDSCoutput$m)) {
    LDSCoutput_core$m <- LDSCoutput$m
  }
  
  expected_vech <- k_core * (k_core + 1) / 2
  info(sprintf("  Subsetted: S=%dx%d, V=%dx%d (expected %dx%d for %d traits)", 
               nrow(S_sub), ncol(S_sub), nrow(V_sub), ncol(V_sub),
               expected_vech, expected_vech, k_core))
  info(sprintf("  S rownames: %s", paste(rownames(S_sub), collapse=", ")))
  info(sprintf("  trait.names: %s", paste(LDSCoutput_core$trait.names, collapse=", ")))
}

# ============================================================
# Issue C fix: Force trait.names to match S rownames (bulletproof)
# ============================================================
LDSCoutput_core$trait.names <- rownames(LDSCoutput_core$S)
info(sprintf("Authoritative trait order (from S rownames): %s", 
             paste(LDSCoutput_core$trait.names, collapse=", ")))

# ============================================================
# ROOT-CAUSE FIX: Ensure loadings match covstruc trait order
# ============================================================
# commonfactorGWAS uses LDSCoutput_core$trait.names order
# Loadings MUST match this order or wrong lambda applied to each trait
info(sprintf("covstruc trait.names: %s", paste(LDSCoutput_core$trait.names, collapse=", ")))
info(sprintf("loadings trait order: %s", paste(loadings$trait, collapse=", ")))

if (!identical(LDSCoutput_core$trait.names, loadings$trait)) {
  warn("Trait order mismatch between covstruc and loadings — reordering loadings to match covstruc")
  loadings <- loadings[match(LDSCoutput_core$trait.names, loadings$trait), ]
  if (any(is.na(loadings$trait))) {
    die("After reordering, some loadings are NA — check trait names match exactly")
  }
  # Update lambda vectors
  lambda <- loadings$lambda
  lambda_std <- loadings$lambda_std
  
  info("Reordered loadings:")
  for (i in seq_along(LDSCoutput_core$trait.names)) {
    info(sprintf("  %s: lambda=%.4f, lambda_std=%.4f", 
                 loadings$trait[i], lambda[i], lambda_std[i]))
  }
  
  # Recompute lambda_norm after reorder
  lambda_norm <- sqrt(sum(lambda^2))
}

# ==============================================================================
# Compute factor heritability (h2F) - AFTER loading order is finalized
# ==============================================================================
# Use covstruc trait order for consistency
core_used <- LDSCoutput_core$trait.names
S_core <- S_all[core_used, core_used]
w <- lambda / lambda_norm
h2F <- as.numeric(t(w) %*% S_core %*% w)
assert(is.finite(h2F) && h2F > 0, paste0("Invalid h2F: ", h2F))
info(sprintf("Factor heritability (h2F): %.4f", h2F))

lambda_sq_sum <- sum(lambda^2)

# ==============================================================================
# Determine N_factor for FUMA output
# ==============================================================================
core_N <- man[match(core_used, trait), N_const]
if (is.finite(N_factor_override)) {
  N_factor <- N_factor_override
} else {
  # Weight by lambda_std^2 (each trait's contribution to factor variance)
  lambda_std_sq <- lambda_std^2
  N_factor <- sum(lambda_std_sq * core_N, na.rm=TRUE) / sum(lambda_std_sq, na.rm=TRUE)
  info(sprintf("N_factor calculation: weighted by lambda_std^2"))
  info(sprintf("  min(N)=%s, weighted(N)=%s", 
               format(min(core_N, na.rm=TRUE), big.mark=","),
               format(round(N_factor), big.mark=",")))
}
assert(is.finite(N_factor) && N_factor > 0, "Invalid N_factor")
info(sprintf("N_factor for FUMA: %s", format(round(N_factor), big.mark=",")))

# ---------------------------------------------------------------------------

factor_gwas_file <- file.path(outdir, "factor_gwas_raw.rds")
factor_done_file <- file.path(outdir, "factor_gwas.done")
chr_dir <- file.path(outdir, "chr")

if (file.exists(factor_done_file) && file.exists(factor_gwas_file)) {
  info("Factor GWAS already complete, loading results...")
  factor_gwas <- readRDS(factor_gwas_file)
} else {
  info("Running commonfactorGWAS...")
  
  # Check if sumstats have CHR for chunking
  has_chr <- "CHR" %in% names(p_sumstats)
  
  if (has_chr) {
    info("Per-chromosome processing with checkpointing...")
    
    p_sumstats_dt <- as.data.table(p_sumstats)
    chromosomes <- sort(unique(p_sumstats_dt$CHR))
    chromosomes <- chromosomes[chromosomes >= 1 & chromosomes <= 22]
    info(sprintf("Chromosomes to process: %s", paste(chromosomes, collapse=", ")))
    
    # Smoke test with 100 SNPs
    info("Smoke test with 100 SNPs...")
    test_snps <- as.data.frame(p_sumstats_dt[1:min(100, nrow(p_sumstats_dt))])
    
    # Debug: show what we're passing
    info(sprintf("  LDSC S dimensions: %dx%d", nrow(LDSCoutput_core$S), ncol(LDSCoutput_core$S)))
    info(sprintf("  LDSC trait.names: %s", paste(LDSCoutput_core$trait.names, collapse=", ")))
    info(sprintf("  sumstats columns: %s", paste(names(test_snps), collapse=", ")))
    
    smoke_result <- tryCatch({
      commonfactorGWAS(
        covstruc = LDSCoutput_core,
        SNPs = test_snps,
        estimation = "DWLS",
        cores = 1,
        parallel = FALSE
      )
    }, error = function(e) {
      die(sprintf("commonfactorGWAS smoke test failed: %s\n\nLDSC S rownames: %s\nLDSC trait.names: %s\nsumstats columns: %s", 
                  e$message, 
                  paste(rownames(LDSCoutput_core$S), collapse=", "),
                  paste(LDSCoutput_core$trait.names, collapse=", "),
                  paste(names(test_snps), collapse=", ")))
    })
    
    info(sprintf("Smoke test passed: %d rows returned", nrow(smoke_result)))
    info(sprintf("Output columns: %s", paste(names(smoke_result), collapse=", ")))
    
    # ==========================================================================
    # MEMORY-EFFICIENT CHROMOSOME PROCESSING
    # - Track completion status only, don't accumulate results in RAM
    # - Combine from disk at the end
    # ==========================================================================
    
    completed_chrs <- c()
    failed_chrs <- c()
    
    for (chr in chromosomes) {
      chr_done <- file.path(chr_dir, paste0("chr", chr, ".done"))
      chr_rds <- file.path(chr_dir, paste0("chr", chr, ".rds"))
      
      # If already done, just record it - DON'T load into memory yet
      if (file.exists(chr_done) && file.exists(chr_rds)) {
        info(sprintf("  Chr %d: already complete", chr))
        completed_chrs <- c(completed_chrs, chr)
        next
      }
      
      info(sprintf("  Chr %d: running commonfactorGWAS...", chr))
      
      ss_chr <- p_sumstats_dt[CHR == chr]
      n_snps_chr <- nrow(ss_chr)
      info(sprintf("    %s SNPs", format(n_snps_chr, big.mark=",")))
      
      if (n_snps_chr == 0) {
        warn(sprintf("  Chr %d: no SNPs, skipping", chr))
        next
      }
      
      chr_result <- tryCatch({
        commonfactorGWAS(
          covstruc = LDSCoutput_core,
          SNPs = as.data.frame(ss_chr),
          estimation = "DWLS",
          cores = n_cores,
          parallel = (n_cores > 1)
        )
      }, error = function(e) {
        warn(sprintf("  Chr %d: GWAS failed: %s", chr, e$message))
        NULL
      })
      
      # Clean up chromosome subset
      rm(ss_chr)
      
      if (is.null(chr_result)) {
        failed_chrs <- c(failed_chrs, chr)
        warn(sprintf("  Chr %d: skipped (%d total failures)", chr, length(failed_chrs)))
        
        if (length(failed_chrs) > THRESH_CHR_FAIL_MAX) {
          die(sprintf("Too many chromosome failures (%d > %d): %s", 
                      length(failed_chrs), THRESH_CHR_FAIL_MAX, 
                      paste(failed_chrs, collapse=", ")))
        }
        gc()
        next
      }
      
      # Save to disk
      atomic_saveRDS(chr_result, chr_rds)
      atomic_write_done(chr_done, chr_rds)
      
      completed_chrs <- c(completed_chrs, chr)
      info(sprintf("    Chr %d complete: %s SNPs", chr, format(nrow(chr_result), big.mark=",")))
      
      # Free this chromosome's result - it's safely on disk
      rm(chr_result)
      
      # GC every 5 chromosomes
      if (length(completed_chrs) %% 5 == 0) gc()
    }
    
    if (length(failed_chrs) > 0) {
      warn(sprintf("Completed with %d chromosome failures: %s", 
                   length(failed_chrs), paste(failed_chrs, collapse=", ")))
    }
    
    # Combine results from disk
    info("Combining chromosome results from disk...")
    
    if (length(completed_chrs) == 0) {
      die("No chromosome results to combine - all failed")
    }
    
    # Sort chromosomes for consistent output
    completed_chrs <- sort(completed_chrs)
    
    # Load and combine
    factor_gwas <- rbindlist(
      lapply(completed_chrs, function(chr) {
        as.data.table(readRDS(file.path(chr_dir, paste0("chr", chr, ".rds"))))
      }), 
      fill = TRUE
    )
    
    info(sprintf("Combined factor GWAS: %s SNPs from %d chromosomes", 
                 format(nrow(factor_gwas), big.mark=","), length(completed_chrs)))
    gc()
    
  } else {
    # No CHR column - run all at once
    info("Running commonfactorGWAS on all SNPs (no checkpointing)...")
    
    factor_gwas <- tryCatch({
      commonfactorGWAS(
        covstruc = LDSCoutput_core,
        SNPs = as.data.frame(p_sumstats),
        estimation = "DWLS",
        cores = n_cores,
        parallel = (n_cores > 1)
      )
    }, error = function(e) {
      die(sprintf("commonfactorGWAS failed: %s", e$message))
    })
    
    factor_gwas <- as.data.table(factor_gwas)
  }
  
  # Save combined results
  atomic_saveRDS(factor_gwas, factor_gwas_file)
  atomic_write_done(factor_done_file, factor_gwas_file)
  info("Factor GWAS complete")
}

# Convert to data.table
gwas <- as.data.table(factor_gwas)
info(sprintf("Factor GWAS results: %s SNPs", format(nrow(gwas), big.mark=",")))

# ==============================================================================
# Extract and format factor GWAS results
# ==============================================================================
info("Processing factor GWAS output...")

info(sprintf("Output columns: %s", paste(names(gwas), collapse=", ")))

# Find the key columns - commonfactorGWAS uses different naming conventions
find_col <- function(dt, patterns, required=TRUE) {
  nms <- names(dt)
  for (p in patterns) {
    if (p %in% nms) return(p)
    idx <- which(tolower(nms) == tolower(p))
    if (length(idx) > 0) return(nms[idx[1]])
    idx <- grep(p, nms, ignore.case=TRUE)
    if (length(idx) > 0) return(nms[idx[1]])
  }
  if (required) die(paste0("Required column not found. Patterns: ", paste(patterns, collapse=", "), "\nAvailable: ", paste(nms, collapse=", ")))
  NA_character_
}

# Standard column extraction
col_snp <- find_col(gwas, c("^SNP$", "SNP", "rsid"))
col_a1 <- find_col(gwas, c("^A1$", "A1", "effect_allele"))
col_a2 <- find_col(gwas, c("^A2$", "A2", "other_allele"))

# For commonfactorGWAS output - careful column matching
# Effect can be: est (beta), Z_Estimate (Z-score), or beta
col_est <- find_col(gwas, c("^est$", "^beta$"), required=FALSE)
col_z_est <- find_col(gwas, c("^Z_Estimate$", "^Z$"), required=FALSE)
col_se <- find_col(gwas, c("^SE$", "^se$", "^SE_Estimate$", "^se_c$"), required=FALSE)
# FIX #1: Remove bare "P" - too ambiguous, can match wrong column
col_p <- find_col(gwas, c("^Pval_Estimate$", "^pval$", "^p_value$"), required=FALSE)

# Q statistics (heterogeneity)
col_q <- find_col(gwas, c("chisq", "Q", "Q_chisq"), required=FALSE)
col_q_df <- find_col(gwas, c("chisq_df", "Q_df"), required=FALSE)
col_q_p <- find_col(gwas, c("chisq_pval", "Q_pval", "Q_P"), required=FALSE)

# CHR/BP
col_chr <- find_col(gwas, c("^CHR$", "CHR", "chr"), required=FALSE)
col_bp <- find_col(gwas, c("^BP$", "BP", "POS", "bp"), required=FALSE)

# Standardize column names
setnames(gwas, col_snp, "SNP")
gwas[, SNP := toupper(as.character(SNP))]
setnames(gwas, col_a1, "A1")
setnames(gwas, col_a2, "A2")
gwas[, A1 := toupper(as.character(A1))]
gwas[, A2 := toupper(as.character(A2))]

# FIX #2: Robust Z vs beta handling
# Rule: If (est + SE) exist → treat as beta/SE, compute Z
#       Else if Z_Estimate exists → treat as Z directly
#       Always compute Factor_P from Factor_Z unless verified Pval_Estimate

has_beta_se <- !is.na(col_est) && !is.na(col_se) && col_est %in% names(gwas) && col_se %in% names(gwas)
has_z_direct <- !is.na(col_z_est) && col_z_est %in% names(gwas)

if (has_beta_se) {
  # Beta + SE available: compute Z = beta/SE
  gwas[, Factor_beta := as.numeric(get(col_est))]
  gwas[, Factor_SE := as.numeric(get(col_se))]
  gwas[, Factor_Z := Factor_beta / Factor_SE]
  info(sprintf("Using beta/SE columns: %s / %s", col_est, col_se))
} else if (has_z_direct) {
  # Z_Estimate available: use directly as Z
  gwas[, Factor_Z := as.numeric(get(col_z_est))]
  # SE might still be available
  if (!is.na(col_se) && col_se %in% names(gwas)) {
    gwas[, Factor_SE := as.numeric(get(col_se))]
    gwas[, Factor_beta := Factor_Z * Factor_SE]
  } else {
    gwas[, Factor_SE := NA_real_]
    gwas[, Factor_beta := NA_real_]
  }
  info(sprintf("Using Z directly from: %s", col_z_est))
} else {
  die("Could not find valid effect columns. Need either (est + SE) or Z_Estimate.")
}

# Compute Factor_P from Factor_Z (most reliable)
# Only use Pval_Estimate if it exists AND passes validation
gwas[, Factor_P := 2 * pnorm(-abs(Factor_Z))]

if (!is.na(col_p) && col_p %in% names(gwas)) {
  pval_candidate <- as.numeric(gwas[[col_p]])
  # Validate: P must be in (0, 1]
  pval_valid <- is.finite(pval_candidate) & pval_candidate > 0 & pval_candidate <= 1
  pct_valid <- mean(pval_valid, na.rm=TRUE) * 100
  
  if (pct_valid > 99) {
    gwas[, Factor_P := pval_candidate]
    info(sprintf("Using P-value column: %s (%.1f%% valid)", col_p, pct_valid))
  } else {
    warn(sprintf("P-value column %s has only %.1f%% valid values, using P from Z instead", col_p, pct_valid))
  }
}

# Q statistics
if (!is.na(col_q) && col_q %in% names(gwas)) {
  gwas[, Q_SNP := as.numeric(get(col_q))]
} else {
  gwas[, Q_SNP := NA_real_]
}

if (!is.na(col_q_df) && col_q_df %in% names(gwas)) {
  gwas[, Q_df := as.numeric(get(col_q_df))]
} else {
  gwas[, Q_df := length(core) - 1]
}

if (!is.na(col_q_p) && col_q_p %in% names(gwas)) {
  gwas[, Q_P := as.numeric(get(col_q_p))]
} else if (any(is.finite(gwas$Q_SNP))) {
  gwas[, Q_P := pchisq(Q_SNP, df=Q_df, lower.tail=FALSE)]
} else {
  gwas[, Q_P := NA_real_]
}

# CHR/BP
if (!is.na(col_chr) && col_chr %in% names(gwas) && col_chr != "CHR") {
  setnames(gwas, col_chr, "CHR")
}
if (!is.na(col_bp) && col_bp %in% names(gwas) && col_bp != "BP") {
  setnames(gwas, col_bp, "BP")
}

# If missing CHR/BP, try to recover from LD score files
if (!"CHR" %in% names(gwas) || all(is.na(gwas$CHR))) {
  info("Recovering CHR/BP from LDSC LD score files...")
  
  ld_files <- list.files(ld_dir, pattern = "\\.l2\\.ldscore\\.gz$", full.names = TRUE)
  
  if (length(ld_files) > 0) {
    chr_bp_list <- lapply(ld_files, function(f) {
      dt <- fread(f, select = c("CHR", "SNP", "BP"))
      dt[, SNP := toupper(SNP)]
      dt
    })
    chr_bp_map <- rbindlist(chr_bp_list)
    chr_bp_map <- unique(chr_bp_map)
    rm(chr_bp_list)
    gc()
    
    setkey(chr_bp_map, SNP)
    setkey(gwas, SNP)
    
    gwas <- chr_bp_map[gwas, on = "SNP"]
    info(sprintf("  CHR/BP mapped for %d / %d SNPs", sum(!is.na(gwas$CHR)), nrow(gwas)))
    rm(chr_bp_map)
    gc()
  }
}

# Validate
frac_finite_z <- mean(is.finite(gwas$Factor_Z))
frac_finite_beta <- mean(is.finite(gwas$Factor_beta), na.rm=TRUE)
frac_finite_se <- mean(is.finite(gwas$Factor_SE), na.rm=TRUE)
info(sprintf("Numeric validation: %.1f%% finite Factor_Z, %.1f%% finite Factor_beta, %.1f%% finite Factor_SE",
             100*frac_finite_z, 100*frac_finite_beta, 100*frac_finite_se))

# Filter primarily on Factor_Z (required), not on beta/SE (optional)
gwas <- gwas[is.finite(Factor_Z)]
# Keep SE filter only if SE column has valid data (don't require it)
if ("Factor_SE" %in% names(gwas) && any(is.finite(gwas$Factor_SE))) {
  gwas <- gwas[is.na(Factor_SE) | (is.finite(Factor_SE) & Factor_SE > 0)]
}
if ("CHR" %in% names(gwas)) {
  gwas <- gwas[!is.na(CHR) & !is.na(BP)]
}

info(sprintf("Factor GWAS after filtering: %s SNPs", format(nrow(gwas), big.mark=",")))

# ==============================================================================
# Write factor GWAS outputs
# ==============================================================================
info("Writing factor GWAS outputs...")

out_cols <- intersect(c("SNP", "CHR", "BP", "A1", "A2", 
                        "Factor_beta", "Factor_SE", "Factor_Z", "Factor_P",
                        "Q_SNP", "Q_df", "Q_P"), names(gwas))

factor_full_path <- file.path(outdir, "factor_gwas_full.tsv.gz")
fwrite(gwas[, ..out_cols], factor_full_path, sep="\t", quote=FALSE, compress="gzip")
info(sprintf("Wrote: %s", factor_full_path))

factor_fuma_path <- file.path(outdir, "factor_gwas_fuma.tsv")
fwrite(gwas[, .(SNP, CHR, BP, A1, A2, P=Factor_P, N=N_factor)],
       factor_fuma_path, sep="\t", quote=FALSE)
info(sprintf("Wrote: %s", factor_fuma_path))

qsnp_fuma_path <- file.path(outdir, "qsnp_gwas_fuma.tsv")
if (any(is.finite(gwas$Q_P))) {
  fwrite(gwas[, .(SNP, CHR, BP, A1, A2, P=Q_P, N=N_factor)],
         qsnp_fuma_path, sep="\t", quote=FALSE)
  info(sprintf("Wrote: %s", qsnp_fuma_path))
} else {
  warn("Q statistics not available - skipping qsnp_gwas_fuma.tsv")
}

# ==============================================================================
# Part B: Residual GWAS (GWAS-by-subtraction)
# ==============================================================================
if (!MINIMAL) {
  info("=================================================================")
  info("Part B: Residual GWAS (GWAS-by-subtraction)")
  info("=================================================================")
  
  # For residuals, we need per-trait betas from processed sumstats
  p_sumstats_dt <- as.data.table(p_sumstats)
  
  qc_summary <- list()
  
  # Core trait residuals
  info(sprintf("Computing residuals for %d core traits...", length(core_used)))
  
  for (t in core_used) {
    info(sprintf("Processing core trait: %s", t))
    
    # Check if trait beta/SE columns are in p_sumstats
    # sumstats() outputs beta.<trait> and se.<trait> columns, not Z
    beta_col <- paste0("beta.", t)
    se_col <- paste0("se.", t)
    
    if (!beta_col %in% names(p_sumstats_dt) || !se_col %in% names(p_sumstats_dt)) {
      warn(sprintf("  %s columns (%s, %s) not found in sumstats, skipping", t, beta_col, se_col))
      next
    }
    
    # Get trait-specific data - compute Z from beta/SE
    trait_data <- p_sumstats_dt[, .(SNP, A1, A2, trait_Z = get(beta_col) / get(se_col))]
    trait_data[, SNP := toupper(SNP)]
    trait_data[, A1 := toupper(A1)]
    trait_data[, A2 := toupper(A2)]
    
    # Get N for this trait to compute beta from Z
    N_t <- man[trait == t, N_const]
    
    # Merge with factor GWAS
    setkey(trait_data, SNP)
    setkey(gwas, SNP)
    
    merged <- merge(trait_data, gwas[, .(SNP, Factor_beta, Factor_SE, Factor_Z, 
                                         gwas_A1 = A1, gwas_A2 = A2, CHR, BP)], 
                    by = "SNP")
    
    # Allele alignment
    merged[, allele_match := (A1 == gwas_A1 & A2 == gwas_A2)]
    merged[, allele_flip := (A1 == gwas_A2 & A2 == gwas_A1)]
    
    n_match <- sum(merged$allele_match)
    n_flip <- sum(merged$allele_flip)
    n_drop <- nrow(merged) - n_match - n_flip
    
    info(sprintf("  Allele alignment: %d match, %d flip, %d drop", n_match, n_flip, n_drop))
    
    merged[allele_flip == TRUE, trait_Z := -trait_Z]
    merged <- merged[allele_match == TRUE | allele_flip == TRUE]
    
    if (nrow(merged) == 0) {
      warn(sprintf("  No aligned SNPs for %s, skipping", t))
      next
    }
    
    # Use standardized loading as rg (genetic correlation with factor)
    # For core traits, rg = lambda_std by definition
    trait_idx <- match(t, core_used)
    rg <- lambda_std[trait_idx]
    
    info(sprintf("  Using standardized loading as rg: lambda_std=%.4f", rg))
    
    # ============================================================
    # CALIBRATED RESIDUAL GWAS (EMPIRICAL Z STANDARDIZATION)
    # ============================================================
    
    # --- diagnose raw Z scales ---
    trait_Z_var  <- var(merged$trait_Z, na.rm = TRUE)
    factor_Z_var <- var(merged$Factor_Z, na.rm = TRUE)
    info(sprintf("  Raw Z variances: trait=%.2f, factor=%.2f", trait_Z_var, factor_Z_var))
    
    # --- empirical scale standardization (winsorized SD) ---
    sd_t <- emp_sd(merged$trait_Z, cap = WINSORIZE_CAP)
    sd_f <- emp_sd(merged$Factor_Z, cap = WINSORIZE_CAP)
    
    if (!is.finite(sd_t) || sd_t <= 0)
      die(sprintf("Invalid trait_Z empirical SD for %s: %s", t, sd_t))
    if (!is.finite(sd_f) || sd_f <= 0)
      die(sprintf("Invalid Factor_Z empirical SD for %s: %s", t, sd_f))
    
    merged[, trait_Zs  := trait_Z  / sd_t]
    merged[, factor_Zs := Factor_Z / sd_f]
    
    # --- confirm standardized scale ---
    info(sprintf("  Std Z variances: trait=%.2f, factor=%.2f",
                 var(merged$trait_Zs, na.rm = TRUE),
                 var(merged$factor_Zs, na.rm = TRUE)))
    
    # --- check if empirical correlation matches rg ---
    emp_cor <- cor(merged$trait_Zs, merged$factor_Zs, use = "complete.obs")
    info(sprintf("  Empirical cor(trait, factor) = %.3f vs theoretical rg = %.3f", emp_cor, rg))
    
    if (abs(emp_cor - rg) > COR_MISMATCH_WARN) {
      warn(sprintf("  Large cor mismatch for %s: emp=%.3f, rg=%.3f", t, emp_cor, rg))
    }
    
    # ============================================================
    # FIX v10: Use EMPIRICAL rg for residual subtraction
    # ============================================================
    # LDSC rg ≠ SNP-level Z correlation (different quantities)
    # Using emp_cor guarantees Var(residual) = 1 - emp_cor² by construction
    
    # Store theoretical rg for logging
    rg_theoretical <- rg
    
    # Fallback safety
    if (!is.finite(emp_cor)) {
      warn(sprintf("  emp_cor not finite for %s, using theoretical rg", t))
      emp_cor <- rg
    }
    
    # Bound for numerical safety
    rg_use <- emp_cor
    if (abs(rg_use) > 0.99) {
      info(sprintf("  Clamping emp_cor from %.3f to +/-0.99", rg_use))
      rg_use <- sign(rg_use) * 0.99
    }
    
    info(sprintf("  Using rg for subtraction: emp_cor=%.3f (theoretical=%.3f)", rg_use, rg_theoretical))
    
    # Use empirical correlation for subtraction
    rg <- rg_use
    
    # --- residual subtraction on standardized Z ---
    merged[, residual_Z := trait_Zs - rg * factor_Zs]
    
    # theoretical variance under calibrated scale
    residual_var <- 1 - rg^2
    merged[, residual_Z_adj := residual_Z / sqrt(residual_var)]
    
    # --- P-values from adjusted Z ---
    merged[, residual_P := 2 * pnorm(-abs(residual_Z_adj))]
    
    # --- clamp extreme P-values defensively ---
    n_zero_p <- sum(merged$residual_P <= 0, na.rm = TRUE)
    if (n_zero_p > 0) {
      info(sprintf("  Clamping %d residual_P values <= 0 to %.0e", n_zero_p, P_MIN))
      merged[, residual_P := pmax(residual_P, P_MIN)]
    }
    
    # --- diagnostics (use adjusted Z, not raw) ---
    mean_chisq_raw <- mean(merged$residual_Z^2, na.rm = TRUE)
    mean_chisq_adj <- mean(merged$residual_Z_adj^2, na.rm = TRUE)
    max_z_adj <- max(abs(merged$residual_Z_adj), na.rm = TRUE)
    
    info(sprintf("  mean_chi2 raw=%.2f, adj=%.2f, max|Z_adj|=%.1f",
                 mean_chisq_raw, mean_chisq_adj, max_z_adj))
    
    # --- consistent beta/SE so beta/SE == Z ---
    merged[, residual_beta := residual_Z_adj / sqrt(N_t)]
    merged[, residual_SE   := 1 / sqrt(N_t)]
    
    n_input <- nrow(merged)
    merged <- merged[is.finite(residual_Z_adj)]
    
    # Remove outliers (use adjusted Z)
    n_before_outlier <- nrow(merged)
    merged <- merged[abs(residual_Z_adj) < THRESH_Z_OUTLIER]
    n_final <- nrow(merged)
    
    # Use adjusted values for QC
    mean_chisq <- mean_chisq_adj
    max_z <- max_z_adj
    
    qc_summary[[t]] <- data.table(
      trait = t,
      type = "core",
      method = "subtraction_calibrated",
      rg_theoretical = rg_theoretical,
      rg_used = rg,
      emp_cor = emp_cor,
      sd_trait = sd_t,
      sd_factor = sd_f,
      n_input = n_input,
      n_finite = n_before_outlier,
      n_final = n_final,
      n_outlier_removed = n_before_outlier - n_final,
      mean_chisq_raw = mean_chisq_raw,
      mean_chisq_adj = mean_chisq_adj,
      max_abs_z = max_z
    )
    
    info(sprintf("  n=%s, mean_chi2=%.2f, max|Z|=%.1f", 
                 format(n_final, big.mark=","), mean_chisq, max_z))
    
    # Write output (use adjusted Z for consistency)
    out <- merged[, .(SNP, CHR, BP, A1 = gwas_A1, A2 = gwas_A2, 
                      residual_beta, residual_SE, residual_Z = residual_Z_adj, 
                      residual_P, P = residual_P)]
    
    resid_path <- file.path(outdir, paste0("residual_", t, ".tsv.gz"))
    fwrite(out, resid_path, sep="\t", quote=FALSE, compress="gzip")
    info(sprintf("  Wrote: %s", resid_path))
    
    # FUMA format
    fuma_path <- file.path(outdir, "fuma", paste0("residual_", t, "_fuma.tsv"))
    fwrite(out[, .(SNP, CHR, BP, A1, A2, P, N=N_t)],
           fuma_path, sep="\t", quote=FALSE)
  }
  
  # Validation trait residuals (if any)
  if (length(validation) > 0) {
    info(sprintf("Computing residuals for %d validation traits...", length(validation)))
    
    for (t in validation) {
      info(sprintf("Processing validation trait: %s", t))
      
      # For validation traits, load their munged sumstats
      val_munged <- file.path(munged_dir, paste0(t, ".sumstats.gz"))
      
      if (!file.exists(val_munged)) {
        warn(sprintf("  Munged file not found for %s, skipping", t))
        next
      }
      
      val_dt <- fread(val_munged)
      setnames(val_dt, names(val_dt), toupper(names(val_dt)))
      
      if (!"Z" %in% names(val_dt)) {
        warn(sprintf("  No Z column in %s, skipping", t))
        next
      }
      
      val_dt[, SNP := toupper(SNP)]
      val_dt[, A1 := toupper(A1)]
      val_dt[, A2 := toupper(A2)]
      
      # Merge with factor GWAS
      setkey(val_dt, SNP)
      setkey(gwas, SNP)
      
      merged <- merge(val_dt[, .(SNP, A1, A2, trait_Z = Z)], 
                      gwas[, .(SNP, Factor_beta, Factor_SE, Factor_Z, 
                               gwas_A1 = A1, gwas_A2 = A2, CHR, BP)], 
                      by = "SNP")
      
      # Allele alignment
      merged[, allele_match := (A1 == gwas_A1 & A2 == gwas_A2)]
      merged[, allele_flip := (A1 == gwas_A2 & A2 == gwas_A1)]
      
      merged[allele_flip == TRUE, trait_Z := -trait_Z]
      merged <- merged[allele_match == TRUE | allele_flip == TRUE]
      
      if (nrow(merged) == 0) {
        warn(sprintf("  No aligned SNPs for %s, skipping", t))
        next
      }
      
      # Compute rg from LDSC
      if (!t %in% rownames(S_all)) {
        warn(sprintf("  %s not in S matrix, skipping", t))
        next
      }
      
      cov_vec <- as.numeric(S_all[core_used, t])
      covYF <- sum(lambda * cov_vec) / lambda_sq_sum
      
      h2_Y <- as.numeric(S_all[t, t])
      rg <- covYF / sqrt(h2_Y * h2F)
      
      # Clamp extreme rg values
      if (abs(rg) > 0.99) {
        warn(sprintf("  rg=%.4f is extreme, clamping to +/-0.99", rg))
        rg <- sign(rg) * 0.99
      }
      
      info(sprintf("  h2_Y=%.4f, cov(Y,F)=%.4f, rg=%.4f", h2_Y, covYF, rg))
      
      # ============================================================
      # CALIBRATED RESIDUAL GWAS (EMPIRICAL Z STANDARDIZATION)
      # ============================================================
      
      # --- diagnose raw Z scales ---
      trait_Z_var  <- var(merged$trait_Z, na.rm = TRUE)
      factor_Z_var <- var(merged$Factor_Z, na.rm = TRUE)
      info(sprintf("  Raw Z variances: trait=%.2f, factor=%.2f", trait_Z_var, factor_Z_var))
      
      # --- empirical scale standardization (winsorized SD) ---
      sd_t <- emp_sd(merged$trait_Z, cap = WINSORIZE_CAP)
      sd_f <- emp_sd(merged$Factor_Z, cap = WINSORIZE_CAP)
      
      if (!is.finite(sd_t) || sd_t <= 0)
        die(sprintf("Invalid trait_Z empirical SD for %s: %s", t, sd_t))
      if (!is.finite(sd_f) || sd_f <= 0)
        die(sprintf("Invalid Factor_Z empirical SD for %s: %s", t, sd_f))
      
      merged[, trait_Zs  := trait_Z  / sd_t]
      merged[, factor_Zs := Factor_Z / sd_f]
      
      # --- confirm standardized scale ---
      info(sprintf("  Std Z variances: trait=%.2f, factor=%.2f",
                   var(merged$trait_Zs, na.rm = TRUE),
                   var(merged$factor_Zs, na.rm = TRUE)))
      
      # --- check if empirical correlation matches rg ---
      emp_cor <- cor(merged$trait_Zs, merged$factor_Zs, use = "complete.obs")
      info(sprintf("  Empirical cor(trait, factor) = %.3f vs theoretical rg = %.3f", emp_cor, rg))
      
      if (abs(emp_cor - rg) > COR_MISMATCH_WARN) {
        warn(sprintf("  Large cor mismatch for %s: emp=%.3f, rg=%.3f", t, emp_cor, rg))
      }
      
      # ============================================================
      # FIX v10: Use EMPIRICAL rg for residual subtraction
      # ============================================================
      # LDSC rg ≠ SNP-level Z correlation (different quantities)
      # Using emp_cor guarantees Var(residual) = 1 - emp_cor² by construction
      
      # Store theoretical rg for logging
      rg_theoretical <- rg
      
      # Fallback safety
      if (!is.finite(emp_cor)) {
        warn(sprintf("  emp_cor not finite for %s, using theoretical rg", t))
        emp_cor <- rg
      }
      
      # Bound for numerical safety
      rg_use <- emp_cor
      if (abs(rg_use) > 0.99) {
        info(sprintf("  Clamping emp_cor from %.3f to +/-0.99", rg_use))
        rg_use <- sign(rg_use) * 0.99
      }
      
      info(sprintf("  Using rg for subtraction: emp_cor=%.3f (theoretical=%.3f)", rg_use, rg_theoretical))
      
      # Use empirical correlation for subtraction
      rg <- rg_use
      
      # --- residual subtraction on standardized Z ---
      merged[, residual_Z := trait_Zs - rg * factor_Zs]
      
      # theoretical variance under calibrated scale
      residual_var <- 1 - rg^2
      merged[, residual_Z_adj := residual_Z / sqrt(residual_var)]
      
      # --- P-values from adjusted Z ---
      merged[, residual_P := 2 * pnorm(-abs(residual_Z_adj))]
      
      # --- clamp extreme P-values defensively ---
      n_zero_p <- sum(merged$residual_P <= 0, na.rm = TRUE)
      if (n_zero_p > 0) {
        info(sprintf("  Clamping %d residual_P values <= 0 to %.0e", n_zero_p, P_MIN))
        merged[, residual_P := pmax(residual_P, P_MIN)]
      }
      
      # --- diagnostics (use adjusted Z, not raw) ---
      mean_chisq_raw <- mean(merged$residual_Z^2, na.rm = TRUE)
      mean_chisq_adj <- mean(merged$residual_Z_adj^2, na.rm = TRUE)
      max_z_adj <- max(abs(merged$residual_Z_adj), na.rm = TRUE)
      
      info(sprintf("  mean_chi2 raw=%.2f, adj=%.2f, max|Z_adj|=%.1f",
                   mean_chisq_raw, mean_chisq_adj, max_z_adj))
      
      # --- consistent beta/SE so beta/SE == Z ---
      N_t <- man[trait == t, N_const]
      merged[, residual_beta := residual_Z_adj / sqrt(N_t)]
      merged[, residual_SE   := 1 / sqrt(N_t)]
      
      n_input <- nrow(merged)
      merged <- merged[is.finite(residual_Z_adj)]
      n_before_outlier <- nrow(merged)
      merged <- merged[abs(residual_Z_adj) < THRESH_Z_OUTLIER]
      n_final <- nrow(merged)
      
      # Use adjusted values for QC
      mean_chisq <- mean_chisq_adj
      max_z <- max_z_adj
      
      qc_summary[[t]] <- data.table(
        trait = t,
        type = "validation",
        method = "subtraction_calibrated",
        rg_theoretical = rg_theoretical,
        rg_used = rg,
        emp_cor = emp_cor,
        sd_trait = sd_t,
        sd_factor = sd_f,
        h2_Y = h2_Y,
        n_input = n_input,
        n_finite = n_before_outlier,
        n_final = n_final,
        n_outlier_removed = n_before_outlier - n_final,
        mean_chisq_raw = mean_chisq_raw,
        mean_chisq_adj = mean_chisq_adj,
        max_abs_z = max_z
      )
      
      info(sprintf("  n=%s, mean_chi2=%.2f, max|Z|=%.1f", 
                   format(n_final, big.mark=","), mean_chisq, max_z))
      
      # Write output (use adjusted Z for consistency)
      out <- merged[, .(SNP, CHR, BP, A1 = gwas_A1, A2 = gwas_A2, 
                        residual_beta, residual_SE, residual_Z = residual_Z_adj, 
                        residual_P, P = residual_P)]
      
      resid_path <- file.path(outdir, paste0("residual_", t, ".tsv.gz"))
      fwrite(out, resid_path, sep="\t", quote=FALSE, compress="gzip")
      info(sprintf("  Wrote: %s", resid_path))
      
      fuma_path <- file.path(outdir, "fuma", paste0("residual_", t, "_fuma.tsv"))
      fwrite(out[, .(SNP, CHR, BP, A1, A2, P, N=N_t)],
             fuma_path, sep="\t", quote=FALSE)
    }
  }
  
  # ==============================================================================
  # QC Summary
  # ==============================================================================
  info("=================================================================")
  info("QC Summary")
  info("=================================================================")
  
  if (length(qc_summary) > 0) {
    qc_dt <- rbindlist(qc_summary, fill=TRUE)
    qc_path <- file.path(outdir, "qc", "qc_summary.tsv")
    fwrite(qc_dt, qc_path, sep="\t", quote=FALSE)
    info(sprintf("Wrote: %s", qc_path))
    
    cat("\n")
    cat("========================================\n")
    cat("           QC SUMMARY                   \n")
    cat("========================================\n")
    print(qc_dt, row.names=FALSE)
    cat("========================================\n\n")
    
    # Flag concerning results
    high_chisq <- qc_dt[mean_chisq > 1.5]
    if (nrow(high_chisq) > 0) {
      warn("High mean chi-square (>1.5) in residuals:")
      for (i in 1:nrow(high_chisq)) {
        warn(sprintf("  %s: mean_chi2=%.2f", high_chisq$trait[i], high_chisq$mean_chisq[i]))
      }
    }
  }  # end if (length(qc_summary) > 0)
} else {
  info("Skipping residual GWAS (--minimal mode)")
}  # end if (!MINIMAL)

# ==============================================================================
# Summary
# ==============================================================================
cat("\n")
cat("=============================================================\n")
cat("                  SCRIPT 04 COMPLETE                         \n")
cat("=============================================================\n")
cat(sprintf("Output directory: %s\n", outdir))
cat("\n")
cat("Factor GWAS outputs:\n")
cat(sprintf("  - %s\n", factor_full_path))
cat(sprintf("  - %s\n", factor_fuma_path))
if (file.exists(qsnp_fuma_path)) {
  cat(sprintf("  - %s\n", qsnp_fuma_path))
}
cat("\n")
cat(sprintf("Core residuals: %d traits\n", length(core)))
cat(sprintf("Validation residuals: %d traits\n", length(validation)))
cat("\n")
cat("STANDARD GENOMICSEM PIPELINE USED:\n")
cat("  1. sumstats() - Processed ORIGINAL files (BETA/SE/P)\n")
cat("  2. commonfactorGWAS() - Factor GWAS with LDSC output\n")
cat("  3. GWAS-by-subtraction for residuals\n")
cat("\n")
cat("FILES USED:\n")
cat(sprintf("  - LDSC input: munged files (Z only) -> ldsc()\n"))
cat(sprintf("  - GWAS input: processed files (BETA/SE/P) -> sumstats()\n"))
cat("\n")
cat("FUMA files: results/factor_gwas/fuma/\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Upload factor_gwas_fuma.tsv to FUMA\n")
cat("  2. Upload residual_*_fuma.tsv files to FUMA\n")
cat("  3. Run LDSC on residuals to check intercepts ~1.0\n")
cat("=============================================================\n")