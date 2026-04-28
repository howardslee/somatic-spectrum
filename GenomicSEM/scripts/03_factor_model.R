#!/usr/bin/env Rscript
# =============================================================================
# VERSION: 2026-01-24-v3 (Standardised pipeline notes)
# Script 03: Fit factor model to core traits (uses ldsc_universe.rds from LDSC munged HM3 inputs)
# FIXED: Added fill=TRUE to first rbind in pair selection logic
# =============================================================================
# FIXES IMPLEMENTED:
#  0) NEW: If S has no row/col names, attempt to recover trait names from ldsc object
#  1) Single tolerance-aware PD/conditioning check (no duplicates)
#  2) Guarded condition number (no divide-by-zero)
#  3) Loadings extraction explicit (lhs=="F") + reordered to match `core`
#  4) Residual correlations labeled "approx" (observed - implied)
#  5) NEW: Formal residual coupling test for high-residual pairs
#  6) FIXED: Added fill=TRUE to first rbind in pair selection to avoid column mismatch
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicSEM)
  library(lavaan)  # For lavTestLRT and parameterEstimates
})

die  <- function(msg) { cat("\n[FAIL] ", msg, "\n\n", sep="", file=stderr()); quit(status=1, save="no") }
warn <- function(msg) cat("[WARN] ", msg, "\n", sep="")
info <- function(msg) cat("[INFO] ", msg, "\n", sep="")
assert <- function(cond, msg) if (!isTRUE(cond)) die(msg)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) die(paste0("Missing value after ", flag))
  args[i + 1]
}

# ------------------------------------------------------------------------------
# Arguments / Paths
# ------------------------------------------------------------------------------
ldsc_rds    <- get_arg("--ldsc-rds", "results/ldsc_universe/ldsc/ldsc_universe.rds")
outdir      <- get_arg("--outdir",   "results/models/factor_model")
core_traits <- get_arg("--core",     "Fibromyalgia,MECFS,IBS,Depression,PTSD,MigAura")
resid_threshold <- as.numeric(get_arg("--resid-threshold", "0.10"))
top_pairs   <- as.integer(get_arg("--top-pairs", "0"))  # Test top N misfit pairs (0 = use threshold only)
apriori_pairs <- get_arg("--apriori-pairs", "")  # Comma-separated pairs to always test, e.g. "IBS-MigAura,Fibromyalgia-MECFS"
test_all_pairs <- "--test-all-pairs" %in% args  # Test all pairs with FDR correction
only_coupling <- "--only-coupling" %in% args

core <- trimws(strsplit(core_traits, ",")[[1]])
info(paste0("Core traits: ", paste(core, collapse=", ")))
info(paste0("Residual coupling threshold: ", resid_threshold))
if (top_pairs > 0) info(paste0("Testing top ", top_pairs, " misfit pairs"))
if (nzchar(apriori_pairs)) info(paste0("A priori pairs to test: ", apriori_pairs))
if (test_all_pairs) info("Testing ALL pairs with FDR correction")

# ------------------------------------------------------------------------------
# --only-coupling mode: load existing factor_model.rds and run coupling tests
# ------------------------------------------------------------------------------
if (only_coupling) {
  info("Running in --only-coupling mode: loading existing factor model...")
  
  rds_path <- file.path(outdir, "factor_model.rds")
  assert(file.exists(rds_path), paste0("factor_model.rds not found at ", rds_path, ". Run full script first."))
  
  factor_model <- readRDS(rds_path)
  
  # Extract what we need
  fit        <- factor_model$fit
  S_core     <- factor_model$S_core
  V_core     <- factor_model$V_core
  core       <- factor_model$core_traits
  n_core     <- length(core)
  loadings_table <- factor_model$loadings
  fit_table  <- factor_model$fit_stats
  resid_table <- factor_model$residual_corr
  model      <- factor_model$model
  
  info(paste0("Loaded existing model: ", model))
  info(paste0("Core traits: ", paste(core, collapse=", ")))
  
  # Calculate n_fail from existing fit_table for summary
  n_pass <- sum(fit_table$pass == "PASS", na.rm=TRUE)
  n_marginal <- sum(fit_table$pass == "MARGINAL", na.rm=TRUE)
  n_fail <- sum(fit_table$pass == "FAIL", na.rm=TRUE)
  
  # Rebuild covstruc for residual coupling tests
  # Load full LDSC to get the structure usermodel expects
  assert(file.exists(ldsc_rds), paste0("LDSC RDS not found: ", ldsc_rds))
  ldsc <- readRDS(ldsc_rds)
  
  covstruc <- ldsc
  covstruc$S <- S_core
  covstruc$V <- V_core
  rownames(covstruc$S) <- colnames(covstruc$S) <- core
  covstruc$S_Stand <- cov2cor(S_core)
  rownames(covstruc$S_Stand) <- colnames(covstruc$S_Stand) <- core
  
  # Subset I if needed (same logic as main script)
  if (!is.null(covstruc$I)) {
    all_traits <- rownames(ldsc$S)
    n_all <- length(all_traits)
    core_idx <- match(core, all_traits)
    I <- covstruc$I
    k_all <- n_all
    expected_len <- k_all + (k_all * (k_all - 1)) / 2
    
    if (length(I) == expected_len) {
      I_mat <- matrix(NA_real_, k_all, k_all)
      diag(I_mat) <- I[1:k_all]
      idx <- k_all + 1
      for (jj in 1:(k_all - 1)) {
        for (ii in (jj + 1):k_all) {
          I_mat[ii, jj] <- I[idx]
          I_mat[jj, ii] <- I[idx]
          idx <- idx + 1
        }
      }
      I_core_mat <- I_mat[core_idx, core_idx, drop = FALSE]
      k_core <- n_core
      I_core <- numeric(k_core + (k_core * (k_core - 1)) / 2)
      I_core[1:k_core] <- diag(I_core_mat)
      idx <- k_core + 1
      for (jj in 1:(k_core - 1)) {
        for (ii in (jj + 1):k_core) {
          I_core[idx] <- I_core_mat[ii, jj]
          idx <- idx + 1
        }
      }
      covstruc$I <- I_core
    }
  }
  
  # Identify high-residual pairs
  high_resid <- resid_table[abs(approx_residual) > resid_threshold]
  if (nrow(high_resid) > 0) {
    info(paste0("High residual pairs to test: ", 
                paste(paste0(high_resid$trait1, "-", high_resid$trait2), collapse=", ")))
  } else {
    info("No high-residual pairs detected.")
  }
  
  # Helper function needed for coupling tests
  pick_col <- function(df, candidates) {
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  
  # Ensure output directory exists
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  
  # Jump to residual coupling tests (the rest of the script from that section)
  # We'll use a flag to skip the main model fitting
  skip_to_coupling <- TRUE
} else {
  skip_to_coupling <- FALSE
}

if (!skip_to_coupling) {
  # BEGIN MAIN SCRIPT (only runs if not --only-coupling)
  
  # ------------------------------------------------------------------------------
  # Load LDSC results
  # ------------------------------------------------------------------------------
  assert(file.exists(ldsc_rds), paste0("LDSC RDS not found: ", ldsc_rds))
  info(paste0("Loading: ", ldsc_rds))
  
  ldsc <- readRDS(ldsc_rds)
  S <- ldsc$S
  V <- ldsc$V
  
  assert(!is.null(S) && is.matrix(S), "ldsc$S is missing or not a matrix")
  assert(!is.null(V) && is.matrix(V), "ldsc$V is missing or not a matrix")
  
  # ------------------------------------------------------------------------------
  # NEW: Recover trait names if dimnames are missing
  # ------------------------------------------------------------------------------
  recover_trait_names <- function(ldsc_obj, S_mat) {
    k <- nrow(S_mat)
    cand <- NULL
    
    # Common locations GenomicSEM / user pipelines store these
    if (!is.null(ldsc_obj$trait.names)) cand <- ldsc_obj$trait.names
    if (is.null(cand) && !is.null(ldsc_obj$traits)) {
      # sometimes a vector, sometimes a named list
      if (is.character(ldsc_obj$traits)) cand <- ldsc_obj$traits
      if (is.list(ldsc_obj$traits) && !is.null(names(ldsc_obj$traits))) cand <- names(ldsc_obj$traits)
    }
    if (is.null(cand) && !is.null(ldsc_obj$input) && !is.null(ldsc_obj$input$trait.names)) cand <- ldsc_obj$input$trait.names
    
    if (!is.null(cand)) {
      cand <- as.character(cand)
      cand <- cand[nzchar(cand)]
      if (length(cand) == k) return(cand)
    }
    NULL
  }
  
  if (is.null(rownames(S)) || is.null(colnames(S)) || any(!nzchar(rownames(S))) || any(!nzchar(colnames(S)))) {
    warn("S matrix has missing/empty row/column names. Attempting to recover trait names from ldsc object...")
    tn <- recover_trait_names(ldsc, S)
    if (is.null(tn)) {
      die(paste0(
        "S matrix has no usable row/column names AND trait names could not be recovered from ldsc object.\n",
        "Fix Script 02 to set row/colnames(S) before saveRDS(), or re-run Script 02.\n",
        "As a last resort, you can manually add names by editing the RDS or adding a trait list to Script 03."
      ))
    }
    rownames(S) <- colnames(S) <- tn
    info(paste0("Recovered and assigned ", length(tn), " trait names to S."))
  }
  
  all_traits <- rownames(S)
  info(paste0("Traits in S matrix: ", paste(all_traits, collapse=", ")))
  
  missing <- setdiff(core, all_traits)
  if (length(missing) > 0) {
    die(paste0("Core traits not found in S matrix: ", paste(missing, collapse=", ")))
  }
  
  # ------------------------------------------------------------------------------
  # Subset S and V to core traits
  # ------------------------------------------------------------------------------
  info("Subsetting S and V matrices to core traits...")
  
  S_core <- S[core, core]
  n_core <- length(core)
  
  get_v_idx <- function(i, j, n) {
    if (i < j) { tmp <- i; i <- j; j <- tmp }
    as.integer((j - 1) * n - j * (j - 1) / 2 + i)
  }
  
  core_idx <- match(core, all_traits)
  n_all <- length(all_traits)
  
  n_core_v <- n_core * (n_core + 1) / 2
  core_v_indices <- integer(n_core_v)
  
  k <- 1
  for (j in 1:n_core) {
    for (i in j:n_core) {
      orig_i <- core_idx[i]
      orig_j <- core_idx[j]
      core_v_indices[k] <- get_v_idx(orig_i, orig_j, n_all)
      k <- k + 1
    }
  }
  
  V_core <- V[core_v_indices, core_v_indices]
  
  info(paste0("S_core dimensions: ", nrow(S_core), " x ", ncol(S_core)))
  info(paste0("V_core dimensions: ", nrow(V_core), " x ", ncol(V_core)))
  
  # ------------------------------------------------------------------------------
  # PD / conditioning check (tolerance-aware; single block)
  # ------------------------------------------------------------------------------
  info("Checking S_core positive-(semi)definiteness and conditioning...")
  
  eigs <- eigen(S_core, symmetric=TRUE, only.values=TRUE)$values
  min_eigen <- min(eigs)
  max_eigen <- max(eigs)
  
  tol <- 1e-6
  cond_number <- max(abs(eigs)) / max(tol, min(abs(eigs)))
  
  info(paste0("Eigenvalue range: ", format(min_eigen, digits=6), " to ", format(max_eigen, digits=6)))
  info(paste0("Condition number (guarded): ", format(cond_number, digits=4)))
  
  if (min_eigen < -1e-4) {
    die(paste0("S_core is not positive semi-definite (min eigenvalue = ",
               format(min_eigen, digits=6), "). Check inputs / trait set."))
  } else if (min_eigen < -tol) {
    warn(paste0("S_core has small negative eigenvalue (", format(min_eigen, digits=6),
                ")—likely sampling noise. Proceeding with caution."))
  } else if (min_eigen <= tol) {
    warn(paste0("S_core is near-singular (min eigenvalue = ", format(min_eigen, digits=6),
                "). Results may be unstable."))
  } else {
    info(paste0("S_core is positive definite (min eigenvalue = ", format(min_eigen, digits=6), ")."))
  }
  
  if (cond_number > 1000) {
    warn(paste0("S_core is poorly conditioned (condition number ≈ ", format(cond_number, digits=4),
                "). Estimates may be unstable; consider alternative models."))
  }
  
  # ------------------------------------------------------------------------------
  # Output directory
  # ------------------------------------------------------------------------------
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  
  # ------------------------------------------------------------------------------
  # Define and fit model
  # ------------------------------------------------------------------------------
  info("Fitting single-factor model...")
  
  model <- paste0("F =~ ", paste(core, collapse=" + "))
  info(paste0("Model: ", model))
  # ------------------------------------------------------------------------------
  # Build GenomicSEM-compatible covstruc (do NOT hand-roll minimal list)
  # ------------------------------------------------------------------------------
  info("Building covstruc for GenomicSEM (subset of LDSC object)...")
  
  # Start from full ldsc object so usermodel() gets expected fields
  covstruc <- ldsc
  covstruc$S <- S_core
  covstruc$V <- V_core
  
  # Ensure S has dimnames
  rownames(covstruc$S) <- colnames(covstruc$S) <- core
  
  # Ensure S_Stand exists and is named
  covstruc$S_Stand <- cov2cor(S_core)
  rownames(covstruc$S_Stand) <- colnames(covstruc$S_Stand) <- core
  
  # ---- Subset intercepts I if present ----
  # GenomicSEM ldsc() returns I as:
  #   - first k entries: univariate intercepts
  #   - then k*(k-1)/2 entries: bivariate intercepts (pairwise)
  # usermodel() may assume I matches S/V ordering.
  if (!is.null(covstruc$I)) {
    I <- covstruc$I
    k_all <- n_all
    expected_len <- k_all + (k_all * (k_all - 1)) / 2
    
    if (length(I) == expected_len) {
      # Build full intercept matrix from I
      I_mat <- matrix(NA_real_, k_all, k_all)
      diag(I_mat) <- I[1:k_all]
      
      idx <- k_all + 1
      for (jj in 1:(k_all - 1)) {
        for (ii in (jj + 1):k_all) {
          I_mat[ii, jj] <- I[idx]
          I_mat[jj, ii] <- I[idx]
          idx <- idx + 1
        }
      }
      
      # Subset to core traits
      I_core_mat <- I_mat[core_idx, core_idx, drop = FALSE]
      
      # Re-pack into GenomicSEM I vector format for k_core traits
      k_core <- n_core
      I_core <- numeric(k_core + (k_core * (k_core - 1)) / 2)
      I_core[1:k_core] <- diag(I_core_mat)
      
      idx <- k_core + 1
      for (jj in 1:(k_core - 1)) {
        for (ii in (jj + 1):k_core) {
          I_core[idx] <- I_core_mat[ii, jj]
          idx <- idx + 1
        }
      }
      
      covstruc$I <- I_core
    } else {
      warn(paste0("I present but unexpected length (", length(I),
                  "). Leaving I unmodified; if usermodel fails, rerun Script 02 with intercepts intact."))
    }
  }
  
  # ------------------------------------------------------------------------------
  # Fit model using usermodel
  # ------------------------------------------------------------------------------
  info("Fitting single-factor model...")
  
  fit <- tryCatch({
    usermodel(
      covstruc = covstruc,
      estimation = "DWLS",
      model = model,
      std.lv = TRUE
    )
  }, error = function(e) {
    die(paste0("Model fitting failed: ", conditionMessage(e)))
  })
  
  info("Model fitted successfully")
  
  
  # ------------------------------------------------------------------------------
  # Fit statistics
  # ------------------------------------------------------------------------------
  info("Extracting fit statistics...")
  
  modelfit <- fit$modelfit
  
  # GenomicSEM returns: chisq, df, p_chisq, AIC, CFI, SRMR
  # TLI and RMSEA are NOT provided by GenomicSEM
  
  # Safe extraction helper
  safe_get <- function(x, name) {
    val <- x[[name]]
    if (is.null(val) || length(val) == 0) return(NA_real_)
    as.numeric(val[1])
  }
  
  chisq_val <- safe_get(modelfit, "chisq")
  df_val <- safe_get(modelfit, "df")
  p_chisq_val <- safe_get(modelfit, "p_chisq")
  aic_val <- safe_get(modelfit, "AIC")
  cfi_val <- safe_get(modelfit, "CFI")
  srmr_val <- safe_get(modelfit, "SRMR")
  
  fit_table <- data.table(
    metric = c("chisq", "df", "p_chisq", "CFI", "SRMR", "AIC"),
    value = c(chisq_val, df_val, p_chisq_val, cfi_val, srmr_val, aic_val),
    threshold = c(NA, NA, "> 0.05", "> 0.95", "< 0.08", NA),
    pass = c(
      NA, NA,
      ifelse(is.na(p_chisq_val), NA, ifelse(p_chisq_val > 0.05, "PASS", "FAIL")),
      ifelse(is.na(cfi_val), NA, ifelse(cfi_val > 0.95, "PASS", ifelse(cfi_val > 0.90, "MARGINAL", "FAIL"))),
      ifelse(is.na(srmr_val), NA, ifelse(srmr_val < 0.08, "PASS", ifelse(srmr_val < 0.10, "MARGINAL", "FAIL"))),
      NA
    )
  )
  
  fwrite(fit_table, file.path(outdir, "model_fit.tsv"), sep="\t")
  info("Saved model_fit.tsv")
  
  n_pass <- sum(fit_table$pass == "PASS", na.rm=TRUE)
  n_marginal <- sum(fit_table$pass == "MARGINAL", na.rm=TRUE)
  n_fail <- sum(fit_table$pass == "FAIL", na.rm=TRUE)
  
  if (n_fail == 0 && n_marginal == 0) {
    info("Overall: EXCELLENT FIT - all criteria met")
  } else if (n_fail == 0) {
    info(paste0("Overall: ACCEPTABLE FIT - ", n_marginal, " marginal criteria"))
  } else {
    warn(paste0("Overall: POOR FIT - ", n_fail, " failed criteria. Consider alternative models."))
  }
  
  # ------------------------------------------------------------------------------
  # Extract factor loadings (robust to character/factor columns)
  # ------------------------------------------------------------------------------
  info("Extracting factor loadings...")
  
  results <- fit$results
  
  # Make sure results is a data.frame-like object
  if (is.null(results) || nrow(results) == 0) {
    die("fit$results is empty. usermodel() returned no results table.")
  }
  
  # Helper: coerce a column to numeric if it exists
  as_num_if_present <- function(df, col) {
    if (col %in% names(df)) {
      df[[col]] <- suppressWarnings(as.numeric(as.character(df[[col]])))
    }
    df
  }
  
  # Coerce common numeric columns used downstream
  for (cc in c("Unstand_Est","Unstand_SE","STD_All","STD_LV","Stand_Est","Stand_SE")) {
    results <- as_num_if_present(results, cc)
  }
  
  # Coerce op/lhs/rhs to character (for reliable filtering)
  for (cc in c("op","lhs","rhs")) {
    if (cc %in% names(results)) results[[cc]] <- as.character(results[[cc]])
  }
  
  # Sanity check required columns exist
  need_cols <- c("op","lhs","rhs")
  missing_cols <- setdiff(need_cols, names(results))
  if (length(missing_cols) > 0) {
    die(paste0("fit$results missing required columns: ", paste(missing_cols, collapse=", ")))
  }
  
  # Filter to factor loadings for latent factor "F"
  loadings_raw <- results[results$op == "=~" & results$lhs == "F", , drop = FALSE]
  
  if (nrow(loadings_raw) == 0) {
    die("No factor loadings found: expected rows with op='=~' and lhs='F' in fit$results.")
  }
  
  # Prefer Unstand_Est / Unstand_SE; fall back if naming differs
  pick_col <- function(df, candidates) {
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  
  col_est <- pick_col(loadings_raw, c("Unstand_Est","Stand_Est","STD_LV","STD_All"))
  col_se  <- pick_col(loadings_raw, c("Unstand_SE","Stand_SE"))
  
  if (is.na(col_est)) die("Could not find an estimate column in loadings_raw (checked Unstand_Est/Stand_Est/STD_LV/STD_All).")
  if (is.na(col_se))  warn("Could not find an SE column in loadings_raw (checked Unstand_SE/Stand_SE). lambda_z/p will be NA.")
  
  # Reorder to match core trait order
  loading_order <- match(core, loadings_raw$rhs)
  if (any(is.na(loading_order))) {
    die(paste0("Missing loadings for traits: ", paste(core[is.na(loading_order)], collapse=", ")))
  }
  loadings_raw <- loadings_raw[loading_order, , drop = FALSE]
  
  # Build table
  lambda     <- loadings_raw[[col_est]]
  lambda_se  <- if (!is.na(col_se)) loadings_raw[[col_se]] else rep(NA_real_, nrow(loadings_raw))
  
  # Standardised loading: prefer STD_All if present, otherwise NA
  lambda_std <- if ("STD_All" %in% names(loadings_raw)) loadings_raw$STD_All else rep(NA_real_, nrow(loadings_raw))
  
  # If STD_All is missing/NA, we can approximate from cov2cor(S_core) only if you want;
  # for now we keep NA to avoid false precision.
  lambda_z <- ifelse(is.finite(lambda) & is.finite(lambda_se) & lambda_se > 0, lambda / lambda_se, NA_real_)
  lambda_p <- ifelse(is.finite(lambda_z), 2 * pnorm(-abs(lambda_z)), NA_real_)
  
  loadings_table <- data.table(
    trait      = loadings_raw$rhs,
    lambda     = lambda,
    lambda_se  = lambda_se,
    lambda_std = lambda_std,
    lambda_z   = lambda_z,
    lambda_p   = lambda_p
  )
  
  fwrite(loadings_table, file.path(outdir, "loadings.tsv"), sep = "\t")
  info("Saved loadings.tsv")
  
  # Print loadings summary
  cat("\n")
  cat("========================================\n")
  cat("          FACTOR LOADINGS               \n")
  cat("========================================\n")
  cat(sprintf("%-15s %10s %10s %10s %12s\n", "Trait", "Lambda", "SE", "Std", "P"))
  cat("----------------------------------------\n")
  for (i in 1:nrow(loadings_table)) {
    cat(sprintf("%-15s %10.4f %10.4f %10.4f %12s\n",
                loadings_table$trait[i],
                loadings_table$lambda[i],
                loadings_table$lambda_se[i],
                loadings_table$lambda_std[i],
                ifelse(is.na(loadings_table$lambda_p[i]), "NA", format(loadings_table$lambda_p[i], scientific=TRUE, digits=3))
    ))
  }
  cat("========================================\n\n")
  
  # ------------------------------------------------------------------------------
  # Residual correlations (approx)
  # ------------------------------------------------------------------------------
  info("Computing residual correlations (approx)...")
  
  lambda_std <- loadings_table$lambda_std
  if (any(lambda_std^2 > 1.01)) warn("Some standardised loadings have |lambda| > 1 (check model / Heywood).")
  
  residual_var <- 1 - lambda_std^2
  if (any(residual_var < -0.01)) warn("Negative residual variances detected (Heywood case risk).")
  
  model_implied <- outer(lambda_std, lambda_std)
  diag(model_implied) <- 1
  
  observed <- cov2cor(S_core)
  
  resid_list <- list()
  kk <- 1
  for (i in 1:(n_core-1)) {
    for (j in (i+1):n_core) {
      resid_list[[kk]] <- data.table(
        trait1 = core[i],
        trait2 = core[j],
        observed_rg = observed[i, j],
        model_implied_rg = model_implied[i, j],
        approx_residual = observed[i, j] - model_implied[i, j]
      )
      kk <- kk + 1
    }
  }
  resid_table <- rbindlist(resid_list)
  resid_table <- resid_table[order(-abs(approx_residual))]
  
  fwrite(resid_table, file.path(outdir, "residual_correlations.tsv"), sep="\t")
  info("Saved residual_correlations.tsv (approx residuals: observed_rg - lambda_i*lambda_j)")
  
  high_resid <- resid_table[abs(approx_residual) > resid_threshold]
  if (nrow(high_resid) > 0) {
    warn(paste0("High residual correlations (|r| > ", resid_threshold, "): ",
                paste(paste0(high_resid$trait1, "-", high_resid$trait2), collapse=", ")))
    warn("May indicate need for correlated residuals or a 2-factor model")
  }
  
} # END if (!skip_to_coupling)

# ------------------------------------------------------------------------------
# Two-pass pair selection for residual coupling tests
# ------------------------------------------------------------------------------
info("Selecting pairs for residual coupling tests...")

# Parse a priori pairs (format: "IBS-MigAura,Fibromyalgia-MECFS")
apriori_list <- list()
if (nzchar(apriori_pairs)) {
  pairs_split <- trimws(strsplit(apriori_pairs, ",")[[1]])
  for (p in pairs_split) {
    traits <- trimws(strsplit(p, "-")[[1]])
    if (length(traits) == 2 && all(traits %in% core)) {
      # Ensure consistent ordering (alphabetical or as in resid_table)
      t1 <- traits[1]; t2 <- traits[2]
      apriori_list[[length(apriori_list) + 1]] <- c(t1, t2)
    } else {
      warn(sprintf("Skipping invalid a priori pair: %s", p))
    }
  }
  if (length(apriori_list) > 0) {
    info(sprintf("A priori pairs to test: %d", length(apriori_list)))
  }
}

# Build list of pairs to test
pairs_to_test <- data.table(trait1 = character(), trait2 = character(), 
                            approx_residual = numeric(), selection_reason = character())

if (test_all_pairs) {
  # Test all pairs
  pairs_to_test <- copy(resid_table)
  pairs_to_test[, selection_reason := "all_pairs"]
  info(sprintf("Testing ALL %d pairs (will apply FDR correction)", nrow(pairs_to_test)))
  
} else {
  # Pass 1: threshold-based selection
  if (nrow(high_resid) > 0) {
    thresh_pairs <- copy(high_resid)
    thresh_pairs[, selection_reason := "threshold"]
    pairs_to_test <- rbind(pairs_to_test, thresh_pairs, fill = TRUE)  # FIXED: Added fill=TRUE
  }
  
  # Pass 1b: top N misfit pairs
  if (top_pairs > 0) {
    top_n <- head(resid_table[order(-abs(approx_residual))], top_pairs)
    top_n[, selection_reason := "top_misfit"]
    pairs_to_test <- rbind(pairs_to_test, top_n, fill = TRUE)
  }
  
  # Pass 1c: a priori pairs
  for (ap in apriori_list) {
    t1 <- ap[1]; t2 <- ap[2]
    # Find in resid_table (check both orderings)
    match_row <- resid_table[(trait1 == t1 & trait2 == t2) | (trait1 == t2 & trait2 == t1)]
    if (nrow(match_row) > 0) {
      match_row[, selection_reason := "apriori"]
      pairs_to_test <- rbind(pairs_to_test, match_row, fill = TRUE)
    }
  }
  
  # Deduplicate (keep first occurrence with its reason)
  pairs_to_test[, pair_key := paste(pmin(trait1, trait2), pmax(trait1, trait2), sep = "_")]
  pairs_to_test <- pairs_to_test[!duplicated(pair_key)]
  pairs_to_test[, pair_key := NULL]
}

n_pairs_to_test <- nrow(pairs_to_test)
if (n_pairs_to_test > 0) {
  info(sprintf("Pairs selected for formal testing: %d", n_pairs_to_test))
  for (i in 1:n_pairs_to_test) {
    info(sprintf("  %s - %s (residual=%.3f, reason=%s)", 
                 pairs_to_test$trait1[i], pairs_to_test$trait2[i],
                 pairs_to_test$approx_residual[i], pairs_to_test$selection_reason[i]))
  }
} else {
  info("No pairs selected for formal testing.")
}

# ------------------------------------------------------------------------------
# Formal residual coupling tests
# ------------------------------------------------------------------------------
info("Running residual coupling tests...")

resid_coupling_results <- list()

if (n_pairs_to_test > 0) {
  cat("\n")
  cat("========================================\n")
  cat("     RESIDUAL COUPLING TESTS            \n")
  cat("========================================\n")
  cat("Testing whether pair-specific shared liability exists\n")
  cat("beyond the core NES factor.\n")
  cat("----------------------------------------\n\n")
  
  for (row_idx in 1:n_pairs_to_test) {
    t1 <- pairs_to_test$trait1[row_idx]
    t2 <- pairs_to_test$trait2[row_idx]
    approx_r <- pairs_to_test$approx_residual[row_idx]
    selection_reason <- pairs_to_test$selection_reason[row_idx]
    
    info(sprintf("Testing %s ~~ %s (approx residual = %.3f, reason = %s)...", 
                 t1, t2, approx_r, selection_reason))
    
    # Build model with freed residual covariance
    model_resid <- paste0(model, "\n", t1, " ~~ ", t2)
    
    fit_resid <- tryCatch({
      usermodel(
        covstruc = covstruc,
        estimation = "DWLS",
        model = model_resid,
        std.lv = TRUE
      )
    }, error = function(e) {
      warn(sprintf("  Failed to fit residual model for %s ~~ %s: %s", t1, t2, conditionMessage(e)))
      NULL
    })
    
    if (!is.null(fit_resid)) {
      # Initialize result variables
      resid_est <- NA_real_
      resid_se <- NA_real_
      resid_z <- NA_real_
      resid_p <- NA_real_
      chisq_diff <- NA_real_
      df_diff <- NA_real_
      p_diff <- NA_real_
      aic_base <- NA_real_
      aic_freed <- NA_real_
      bic_base <- NA_real_
      bic_freed <- NA_real_
      delta_aic <- NA_real_
      delta_bic <- NA_real_
      lrt_method <- "none"
      
      # --- Try lavaan functions first, fall back to manual extraction ---
      
      # 1. Extract residual covariance estimate
      pe <- tryCatch(parameterEstimates(fit_resid), error = function(e) NULL)
      if (!is.null(pe)) {
        rc <- pe[pe$op == "~~" & 
                   ((pe$lhs == t1 & pe$rhs == t2) | (pe$lhs == t2 & pe$rhs == t1)), ]
        if (nrow(rc) >= 1) {
          resid_est <- rc$est[1]
          resid_se <- rc$se[1]
          if (!is.na(resid_est) && !is.na(resid_se) && resid_se > 0) {
            resid_z <- resid_est / resid_se
            resid_p <- 2 * pnorm(-abs(resid_z))
          }
        }
        lrt_method <- "lavaan"
      } else {
        # Fallback: extract from fit$results
        resid_row <- fit_resid$results[
          fit_resid$results$op == "~~" & 
            ((fit_resid$results$lhs == t1 & fit_resid$results$rhs == t2) |
               (fit_resid$results$lhs == t2 & fit_resid$results$rhs == t1)), 
        ]
        if (nrow(resid_row) > 0) {
          est_col <- pick_col(resid_row, c("Unstand_Est", "Stand_Est", "STD_All"))
          se_col <- pick_col(resid_row, c("Unstand_SE", "Stand_SE"))
          if (!is.na(est_col)) resid_est <- as.numeric(resid_row[[est_col]][1])
          if (!is.na(se_col)) resid_se <- as.numeric(resid_row[[se_col]][1])
          if (!is.na(resid_est) && !is.na(resid_se) && resid_se > 0) {
            resid_z <- resid_est / resid_se
            resid_p <- 2 * pnorm(-abs(resid_z))
          }
        }
        lrt_method <- "manual"
      }
      
      # 2. Likelihood ratio test (chi-square difference)
      lrt <- tryCatch(lavTestLRT(fit, fit_resid), error = function(e) NULL)
      if (!is.null(lrt) && nrow(lrt) >= 2) {
        # Second row is the comparison (base vs freed)
        chisq_diff <- as.numeric(lrt[2, "Chisq diff"])
        df_diff <- as.numeric(lrt[2, "Df diff"])
        p_diff <- as.numeric(lrt[2, "Pr(>Chisq)"])
        lrt_method <- "lavTestLRT"
      } else {
        # Fallback: manual chi-square difference
        chisq_base_tmp <- fit$modelfit$chisq
        chisq_freed_tmp <- fit_resid$modelfit$chisq
        df_base_tmp <- fit$modelfit$df
        df_freed_tmp <- fit_resid$modelfit$df
        
        chisq_diff <- chisq_base_tmp - chisq_freed_tmp
        df_diff <- df_base_tmp - df_freed_tmp
        
        if (df_diff <= 0) {
          warn(sprintf("  Unexpected df difference (%d). Setting to 1.", df_diff))
          df_diff <- 1
        }
        p_diff <- pchisq(chisq_diff, df = df_diff, lower.tail = FALSE)
        lrt_method <- paste0(lrt_method, "+manual_chisq")
      }
      
      # 3. AIC/BIC - try lavaan functions first
      aic_base <- tryCatch(AIC(fit), error = function(e) fit$modelfit$AIC)
      aic_freed <- tryCatch(AIC(fit_resid), error = function(e) fit_resid$modelfit$AIC)
      bic_base <- tryCatch(BIC(fit), error = function(e) NA_real_)
      bic_freed <- tryCatch(BIC(fit_resid), error = function(e) NA_real_)
      
      delta_aic <- aic_freed - aic_base  # Negative = freed model better
      delta_bic <- if (!is.na(bic_base) && !is.na(bic_freed)) bic_freed - bic_base else NA_real_
      
      # Store raw chi-square values for reference
      chisq_base <- if (!is.null(fit$modelfit$chisq)) fit$modelfit$chisq else NA_real_
      chisq_freed <- if (!is.null(fit_resid$modelfit$chisq)) fit_resid$modelfit$chisq else NA_real_
      
      # Determine significance and interpretation
      is_significant <- !is.na(p_diff) && p_diff < 0.05
      is_aic_better <- !is.na(delta_aic) && delta_aic < -2
      is_bic_better <- !is.na(delta_bic) && delta_bic < -2
      
      interpretation <- if (is_significant && is_aic_better) {
        "SIGNIFICANT: Pair shares liability beyond factor"
      } else if (is_significant) {
        "MARGINAL: χ² significant but AIC not improved"
      } else {
        "NOT SIGNIFICANT: Factor fully explains covariance"
      }
      
      # Store results (significance will be recalculated with FDR if needed)
      resid_coupling_results[[length(resid_coupling_results) + 1]] <- data.table(
        trait1 = t1,
        trait2 = t2,
        selection_reason = selection_reason,
        approx_residual = approx_r,
        resid_cov_est = resid_est,
        resid_cov_se = resid_se,
        resid_cov_z = resid_z,
        resid_cov_p = resid_p,
        chisq_base = chisq_base,
        chisq_freed = chisq_freed,
        chisq_diff = chisq_diff,
        df_diff = df_diff,
        p_chisq_diff = p_diff,
        aic_base = aic_base,
        aic_freed = aic_freed,
        delta_aic = delta_aic,
        bic_base = bic_base,
        bic_freed = bic_freed,
        delta_bic = delta_bic,
        lrt_method = lrt_method,
        significant_nominal = is_significant,
        interpretation_nominal = interpretation
      )
      
      # Print results (nominal, FDR correction applied after all tests)
      cat(sprintf("  %s ~~ %s:\n", t1, t2))
      cat(sprintf("    Residual covariance: %.4f (SE=%.4f, Z=%.2f, p=%.4f)\n", 
                  resid_est, resid_se, resid_z, resid_p))
      cat(sprintf("    χ² difference test: Δχ²=%.2f, df=%d, p=%.4f (nominal, method=%s)\n", 
                  chisq_diff, df_diff, p_diff, lrt_method))
      cat(sprintf("    AIC: ΔAIC=%.2f (%s)\n", 
                  delta_aic, ifelse(delta_aic < -2, "freed model better", "base model adequate")))
      if (!is.na(delta_bic)) {
        cat(sprintf("    BIC: ΔBIC=%.2f (%s)\n", 
                    delta_bic, ifelse(delta_bic < -2, "freed model better", "base model adequate")))
      }
      cat(sprintf("    → %s (nominal)\n\n", interpretation))
    }
  }
  
  cat("========================================\n\n")
} else {
  info("No pairs selected for testing. Single factor appears sufficient.")
}

# Compile and save residual coupling results
if (length(resid_coupling_results) > 0) {
  resid_coupling_table <- rbindlist(resid_coupling_results)
  
  # Apply FDR correction (Benjamini-Hochberg)
  n_tests <- nrow(resid_coupling_table)
  if (n_tests > 1) {
    resid_coupling_table[, p_chisq_fdr := p.adjust(p_chisq_diff, method = "BH")]
    info(sprintf("Applied BH-FDR correction across %d tests", n_tests))
    
    # Final significance uses FDR-corrected p-values
    resid_coupling_table[, significant := !is.na(p_chisq_fdr) & p_chisq_fdr < 0.05]
    resid_coupling_table[, interpretation := ifelse(
      significant & delta_aic < -2,
      "SIGNIFICANT (FDR): Pair shares liability beyond factor",
      ifelse(significant,
             "MARGINAL (FDR): χ² significant but AIC not improved",
             "NOT SIGNIFICANT (FDR): Factor fully explains covariance")
    )]
    
    # Print FDR summary
    cat("\n")
    cat("========================================\n")
    cat("     FDR-CORRECTED RESULTS              \n")
    cat("========================================\n")
    for (i in 1:nrow(resid_coupling_table)) {
      r <- resid_coupling_table[i]
      cat(sprintf("  %s ~~ %s: p_nominal=%.4f, p_FDR=%.4f → %s\n",
                  r$trait1, r$trait2, r$p_chisq_diff, r$p_chisq_fdr,
                  ifelse(r$significant, "SIGNIFICANT", "not significant")))
    }
    cat("========================================\n\n")
    
  } else {
    # Single test: no FDR needed
    resid_coupling_table[, p_chisq_fdr := p_chisq_diff]
    resid_coupling_table[, significant := significant_nominal]
    resid_coupling_table[, interpretation := interpretation_nominal]
  }
  
  fwrite(resid_coupling_table, file.path(outdir, "residual_coupling_tests.tsv"), sep = "\t")
  info("Saved residual_coupling_tests.tsv")
  
  # Count significant pairs
  n_sig <- sum(resid_coupling_table$significant, na.rm = TRUE)
  if (n_sig > 0) {
    sig_pairs <- resid_coupling_table[significant == TRUE, paste0(trait1, "-", trait2)]
    warn(sprintf("%d pair(s) show significant residual coupling (FDR<0.05): %s", n_sig, paste(sig_pairs, collapse = ", ")))
    warn("Consider bifactor model or correlated residuals in final specification.")
  } else {
    info("No pairs significant after FDR correction. Single factor model adequate.")
  }
} else {
  # Create empty file for consistency
  resid_coupling_table <- data.table(
    trait1 = character(),
    trait2 = character(),
    selection_reason = character(),
    approx_residual = numeric(),
    resid_cov_est = numeric(),
    resid_cov_se = numeric(),
    resid_cov_z = numeric(),
    resid_cov_p = numeric(),
    chisq_base = numeric(),
    chisq_freed = numeric(),
    chisq_diff = numeric(),
    df_diff = integer(),
    p_chisq_diff = numeric(),
    aic_base = numeric(),
    aic_freed = numeric(),
    delta_aic = numeric(),
    bic_base = numeric(),
    bic_freed = numeric(),
    delta_bic = numeric(),
    lrt_method = character(),
    significant_nominal = logical(),
    interpretation_nominal = character(),
    p_chisq_fdr = numeric(),
    significant = logical(),
    interpretation = character()
  )
  fwrite(resid_coupling_table, file.path(outdir, "residual_coupling_tests.tsv"), sep = "\t")
  info("Saved residual_coupling_tests.tsv (empty - no pairs tested)")
}

# ------------------------------------------------------------------------------
# Save full model summary
# ------------------------------------------------------------------------------
sink(file.path(outdir, "model_summary.txt"))
cat("=============================================================\n")
cat("         NES SOMATIC CORE FACTOR MODEL SUMMARY               \n")
cat("=============================================================\n\n")
cat("Model specification:\n")
cat(paste0("  ", model, "\n\n"))
cat("Traits:\n")
for (t in core) cat(paste0("  - ", t, "\n"))
cat("\n-------------------------------------------------------------\n")
cat("FIT STATISTICS\n")
cat("-------------------------------------------------------------\n")
print(fit_table, row.names=FALSE)
cat("\n-------------------------------------------------------------\n")
cat("FACTOR LOADINGS\n")
cat("-------------------------------------------------------------\n")
print(loadings_table, row.names=FALSE)
cat("\n-------------------------------------------------------------\n")
cat("RESIDUAL CORRELATIONS (approx)\n")
cat("-------------------------------------------------------------\n")
print(resid_table, row.names=FALSE)
cat("\n-------------------------------------------------------------\n")
cat("RESIDUAL COUPLING TESTS\n")
cat("-------------------------------------------------------------\n")
if (nrow(resid_coupling_table) > 0) {
  cat("\nPairs tested:\n")
  cat(sprintf("  Selection criteria: threshold=%.2f, top_pairs=%d, test_all=%s\n", 
              resid_threshold, top_pairs, test_all_pairs))
  cat(sprintf("  Total pairs tested: %d\n", nrow(resid_coupling_table)))
  cat(sprintf("  FDR correction: %s\n\n", ifelse(nrow(resid_coupling_table) > 1, "BH (applied)", "not needed")))
  
  for (i in 1:nrow(resid_coupling_table)) {
    r <- resid_coupling_table[i]
    cat(sprintf("  %s ~~ %s (selected: %s, method: %s):\n", r$trait1, r$trait2, r$selection_reason, r$lrt_method))
    cat(sprintf("    Approx residual: %.4f\n", r$approx_residual))
    cat(sprintf("    Residual covariance: %.4f (SE=%.4f, p=%.4f)\n", 
                r$resid_cov_est, r$resid_cov_se, r$resid_cov_p))
    cat(sprintf("    χ² difference: %.2f (df=%d, p_nominal=%.4f, p_FDR=%.4f)\n", 
                r$chisq_diff, r$df_diff, r$p_chisq_diff, r$p_chisq_fdr))
    cat(sprintf("    ΔAIC: %.2f\n", r$delta_aic))
    if (!is.na(r$delta_bic)) cat(sprintf("    ΔBIC: %.2f\n", r$delta_bic))
    cat(sprintf("    Interpretation: %s\n\n", r$interpretation))
  }
} else {
  cat("No pairs selected for testing.\n")
  cat("Single factor model appears sufficient.\n")
}
cat("\n-------------------------------------------------------------\n")
cat("FULL MODEL OUTPUT\n")
cat("-------------------------------------------------------------\n")
print(fit$results)
cat("\n")
sink()
info("Saved model_summary.txt")

# ------------------------------------------------------------------------------
# ASCII path diagram
# ------------------------------------------------------------------------------
info("Creating path diagram...")

sink(file.path(outdir, "path_diagram.txt"))
cat("\n")
cat("                    NES SOMATIC CORE FACTOR MODEL\n")
cat("                    =============================\n\n")
cat("                              [ F ]\n")
cat("                           (Latent Factor)\n")
cat("                                |\n")
cat("        +----------+----------+---+----------+----------+\n")
cat("        |          |          |   |          |          |\n")
cat("        v          v          v   v          v          v\n\n")
for (i in 1:nrow(loadings_table)) {
  cat(sprintf("    [%s]\n", loadings_table$trait[i]))
  cat(sprintf("     λ=%.3f\n\n", loadings_table$lambda_std[i]))
}

# Add residual coupling if significant
if (nrow(resid_coupling_table) > 0 && any(resid_coupling_table$significant)) {
  cat("\n")
  cat("Significant residual covariances (FDR<0.05, beyond factor):\n")
  for (i in 1:nrow(resid_coupling_table)) {
    if (resid_coupling_table$significant[i]) {
      cat(sprintf("    %s <~~> %s (r=%.3f, p_FDR=%.4f)\n", 
                  resid_coupling_table$trait1[i],
                  resid_coupling_table$trait2[i],
                  resid_coupling_table$resid_cov_est[i],
                  resid_coupling_table$p_chisq_fdr[i]))
    }
  }
}

cat("\nLegend:\n")
cat("  F = Latent somatic core factor\n")
cat("  λ = Standardised factor loading\n")
if (nrow(resid_coupling_table) > 0 && any(resid_coupling_table$significant)) {
  cat("  <~~> = Significant residual covariance beyond factor\n")
}
sink()

# ------------------------------------------------------------------------------
# Save objects for Script 04
# ------------------------------------------------------------------------------
info("Saving objects for downstream scripts...")

factor_model <- list(
  model = model,
  fit = fit,
  S_core = S_core,
  V_core = V_core,
  core_traits = core,
  loadings = loadings_table,
  fit_stats = fit_table,
  residual_corr = resid_table,
  residual_coupling = resid_coupling_table
)

saveRDS(factor_model, file.path(outdir, "factor_model.rds"))
info("Saved factor_model.rds")

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
cat("\n=============================================================\n")
cat("                    SCRIPT 03 COMPLETE                       \n")
cat("=============================================================\n")
cat(paste0("Output directory: ", outdir, "\n\n"))

# Determine overall recommendation
n_sig_coupling <- if (nrow(resid_coupling_table) > 0) sum(resid_coupling_table$significant, na.rm = TRUE) else 0

if (n_fail == 0 && n_sig_coupling == 0) {
  cat("✓ Model fit acceptable and no significant residual coupling\n")
  cat("  → Single factor model is adequate. Proceed to Script 04 (Factor GWAS)\n")
} else if (n_fail == 0 && n_sig_coupling > 0) {
  cat("⚠ Model fit acceptable BUT significant residual coupling detected\n")
  cat("  Consider:\n")
  cat("  1. Proceeding with base model (conservative: shared liability only)\n")
  cat("  2. Adding freed residual covariances to model specification\n")
  cat("  3. Testing bifactor model for mechanistic substructure\n")
  cat("  See residual_coupling_tests.tsv for details.\n")
} else {
  cat("⚠ Model fit poor - consider:\n")
  cat("  1. Removing low-loading traits\n")
  cat("  2. Adding correlated residuals\n")
  cat("  3. Testing 2-factor model\n")
}
cat("=============================================================\n")
quit(status=0, save="no")