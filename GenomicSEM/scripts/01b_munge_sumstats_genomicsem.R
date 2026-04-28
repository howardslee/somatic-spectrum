#!/usr/bin/env Rscript
# =============================================================================
# Script 01b: Munge GWAS summary statistics using GenomicSEM (R-native)
#
# Replaces the Python LDSC munge_sumstats.py with GenomicSEM's built-in munge()
#
# Input:  data/processed/hg19/<TRAIT>/<TRAIT>.hg19.hm3.sumstats.tsv.gz
# Output: results/ldsc_universe/munged/<TRAIT>.sumstats.gz
#
# Usage:
#   Rscript 01b_munge_sumstats_genomicsem.R
#   Rscript 01b_munge_sumstats_genomicsem.R --trait ADHD
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
})

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
IN_ROOT  <- "data/processed/hg19"
OUT_DIR  <- "results/ldsc_universe/munged"
HM3_FILE <- "data/reference/w_hm3.snplist"

# Parse command line args
args <- commandArgs(trailingOnly = TRUE)
only_trait <- NULL
if ("--trait" %in% args) {
  idx <- which(args == "--trait")
  if (idx < length(args)) {
    only_trait <- args[idx + 1]
  }
}

# Create output directory
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Check HM3 file exists
if (!file.exists(HM3_FILE)) {
  stop("[ERROR] HM3 SNP list not found: ", HM3_FILE)
}

cat("[INFO] Using HM3 SNP list:", HM3_FILE, "\n")
cat("[INFO] Input root:", IN_ROOT, "\n")
cat("[INFO] Output dir:", OUT_DIR, "\n\n")

# ------------------------------------------------------------------------------
# Find all trait input files
# ------------------------------------------------------------------------------
input_files <- list.files(
  path = IN_ROOT,
  pattern = "\\.hg19\\.hm3\\.sumstats\\.tsv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(input_files) == 0) {
  stop("[ERROR] No input files found matching *.hg19.hm3.sumstats.tsv.gz in ", IN_ROOT)
}

# Extract trait names from paths
traits <- sapply(input_files, function(f) {
  basename(dirname(f))
})
names(input_files) <- traits

# Filter to single trait if requested
if (!is.null(only_trait)) {
  if (!only_trait %in% traits) {
    stop("[ERROR] Trait not found: ", only_trait)
  }
  input_files <- input_files[only_trait]
  traits <- only_trait
}

cat("[INFO] Found", length(input_files), "traits to process\n\n")

# ------------------------------------------------------------------------------
# Process each trait
# ------------------------------------------------------------------------------
fails <- 0

for (i in seq_along(input_files)) {
  trait <- names(input_files)[i]
  infile <- input_files[i]
  
  cat("-------------------------------------------------------------------------------\n")
  cat("[INFO] Trait:", trait, "\n")
  cat("[INFO] Input:", infile, "\n")
  
  # Check if already done
  outfile <- file.path(OUT_DIR, paste0(trait, ".sumstats.gz"))
  if (file.exists(outfile)) {
    cat("[INFO] Output already exists, skipping:", outfile, "\n")
    next
  }
  
  # Read header to check columns
  header <- names(fread(infile, nrows = 0))
  cat("[INFO] Columns:", paste(header, collapse = ", "), "\n")
  
  # Determine N handling
  # GenomicSEM munge() needs either:
  #   - N parameter (single value)
  #   - N column in file (auto-detected)
  #   - N.cas and N.con for case/control
  
  has_n <- "N" %in% header
  has_ncase <- "N_CASES" %in% header
  has_ncon <- "N_CONTROLS" %in% header
  
  if (has_ncase && has_ncon) {
    cat("[INFO] Using N_CASES / N_CONTROLS columns\n")
  } else if (has_n) {
    cat("[INFO] Using N column\n")
  } else {
    cat("[ERROR] No N column found (expected N or N_CASES/N_CONTROLS)\n")
    fails <- fails + 1
    next
  }
  
  # Run GenomicSEM munge
  # Note: munge() writes directly to working directory with trait name
  # We'll run it and then move the output
  
  tryCatch({
    # GenomicSEM munge expects specific column names
    # Our files have: SNP, CHR, BP, A1, A2, BETA, SE, P, N, N_CASES, N_CONTROLS
    # munge() auto-detects these standard names
    
    munge(
      files = infile,
      hm3 = HM3_FILE,
      trait.names = trait,
      N = NULL,  # Use N from file
      info.filter = 0.9,
      maf.filter = 0.01
    )
    
    # munge() creates <trait>.sumstats.gz in current directory
    # Move to output directory
    temp_out <- paste0(trait, ".sumstats.gz")
    if (file.exists(temp_out)) {
      file.rename(temp_out, outfile)
      cat("[INFO] Created:", outfile, "\n")
      
      # Report SNP count
      n_snps <- nrow(fread(outfile, select = 1))
      cat("[INFO] Retained SNPs:", format(n_snps, big.mark = ","), "\n")
    } else {
      cat("[ERROR] munge() did not create expected output\n")
      fails <- fails + 1
    }
    
    # Clean up log file
    log_file <- paste0(trait, ".log")
    if (file.exists(log_file)) {
      file.rename(log_file, file.path(OUT_DIR, log_file))
    }
    
  }, error = function(e) {
    cat("[ERROR] munge() failed:", conditionMessage(e), "\n")
    fails <<- fails + 1
  })
}

cat("-------------------------------------------------------------------------------\n")
if (fails > 0) {
  cat("[DONE] Completed with", fails, "failures\n")
  quit(status = 1)
} else {
  cat("[DONE] All traits munged successfully into:", OUT_DIR, "\n")
}