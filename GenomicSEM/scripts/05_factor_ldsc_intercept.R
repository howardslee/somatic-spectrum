#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicSEM)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i)) return(default)
  if (i[length(i)] == length(args)) stop("Missing value after ", flag)
  args[i[length(i)] + 1]
}

model_dir <- get_arg("--model-dir", "results/models/somatic5")
ld_dir <- get_arg("--ld", "data/reference/eur_w_ld_chr")
outdir <- get_arg("--outdir", file.path(model_dir, "ldsc_factor"))
trait_name <- get_arg("--trait-name", "Somatic5_factor")

factor_full <- get_arg("--factor-full", file.path(model_dir, "factor_gwas_full.tsv.gz"))
factor_fuma <- get_arg("--factor-fuma", file.path(model_dir, "factor_gwas_fuma.tsv"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(factor_full))
stopifnot(file.exists(factor_fuma))
stopifnot(dir.exists(ld_dir))

message("Reading factor GWAS: ", factor_full)
gwas <- fread(factor_full, select = c("SNP", "A1", "A2", "Factor_Z"))
setnames(gwas, "Factor_Z", "Z")
gwas[, SNP := sub("^RS", "rs", SNP)]

message("Reading factor effective N: ", factor_fuma)
fuma_n <- fread(factor_fuma, select = "N")
n_factor <- unique(na.omit(fuma_n$N))
if (length(n_factor) != 1L) {
  stop("Expected exactly one non-missing N in factor FUMA file; found ", length(n_factor))
}

gwas[, N := as.numeric(n_factor)]
gwas <- gwas[is.finite(Z) & is.finite(N) & !is.na(SNP) & !is.na(A1) & !is.na(A2)]
setcolorder(gwas, c("SNP", "N", "Z", "A1", "A2"))

factor_munged <- file.path(outdir, paste0(trait_name, ".sumstats.gz"))
message("Writing munged factor LDSC input: ", factor_munged)
fwrite(gwas, factor_munged, sep = "\t", quote = FALSE, compress = "gzip")

ldsc_prefix <- file.path(outdir, trait_name)
message("Running GenomicSEM::ldsc()")
ldsc_out <- ldsc(
  traits = factor_munged,
  sample.prev = NA,
  population.prev = NA,
  ld = ld_dir,
  wld = ld_dir,
  trait.names = trait_name,
  ldsc.log = ldsc_prefix,
  chr = 22
)

saveRDS(ldsc_out, file.path(outdir, paste0(trait_name, "_ldsc.rds")))

log_path <- paste0(ldsc_prefix, "_ldsc.log")
log_lines <- readLines(log_path, warn = FALSE)
pick <- function(pattern) {
  hit <- grep(pattern, log_lines, value = TRUE)
  if (length(hit)) hit[1] else NA_character_
}

summary <- data.table(
  trait = trait_name,
  n_factor = n_factor,
  n_snps_input = nrow(gwas),
  mean_chi2_line = pick("^Mean Chi\\^2"),
  lambda_gc_line = pick("^Lambda GC"),
  intercept_line = pick("^Intercept:"),
  ratio_line = pick("^Ratio:"),
  h2_line = pick("^Total Observed Scale h2:"),
  log_path = log_path,
  rds_path = file.path(outdir, paste0(trait_name, "_ldsc.rds")),
  sumstats_path = factor_munged
)

summary_path <- file.path(outdir, paste0(trait_name, "_ldsc_summary.tsv"))
fwrite(summary, summary_path, sep = "\t", quote = FALSE)

message("Wrote summary: ", summary_path)
print(summary)
