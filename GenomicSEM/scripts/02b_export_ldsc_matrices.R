#!/usr/bin/env Rscript
# =============================================================================
# Script 02b: Export LDSC matrices and visualizations
# =============================================================================
# Run after script 02 creates ldsc_universe.rds
#
# Outputs:
#   results/ldsc_universe/matrices/
#     S_genetic_covariance.tsv         - Genetic covariance matrix
#     rg_genetic_correlation.tsv       - Genetic correlation matrix
#     rg_clustered.tsv                 - Correlation matrix (clustered order)
#     heritabilities.tsv               - h² estimates per trait
#     intercepts.tsv                   - LDSC intercepts
#     pairwise_correlations.tsv        - All pairwise rg (sorted by |rg|)
#     rg_heatmap.pdf / .png            - Clustered heatmap
#     dendrogram.pdf                   - Hierarchical clustering tree
#     summary_report.txt               - Text summary
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# Check for visualization packages
has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)
has_corrplot <- requireNamespace("corrplot", quietly = TRUE)

if (!has_pheatmap && !has_corrplot) {
  message("[INFO] For better heatmaps: install.packages(c('pheatmap', 'corrplot'))")
}

# ------------------------------------------------------------------------------
# Arguments
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i[length(i)] == length(args)) stop(paste0("Missing value after ", flag))
  args[i[length(i)] + 1]
}

ldsc_rds <- get_arg("--ldsc-rds", "results/ldsc_universe/ldsc/ldsc_universe.rds")
outdir   <- get_arg("--outdir", "results/ldsc_universe/matrices")

# ------------------------------------------------------------------------------
# Load LDSC results
# ------------------------------------------------------------------------------
if (!file.exists(ldsc_rds)) {
  stop("[FAIL] LDSC RDS not found: ", ldsc_rds, "\n  Run script 02 first.")
}

cat("[INFO] Loading:", ldsc_rds, "\n")
ldsc <- readRDS(ldsc_rds)

S <- ldsc$S
V <- ldsc$V
I <- ldsc$I

# Get trait names - try multiple locations
traits <- rownames(S)
if (is.null(traits) || length(traits) == 0) {
  traits <- ldsc$trait.names
}
if (is.null(traits) || length(traits) == 0) {
  traits <- colnames(S)
}
if (is.null(traits) || length(traits) == 0) {
  stop("[FAIL] Cannot find trait names in LDSC object. Check ldsc$trait.names or rownames(ldsc$S)")
}

n_traits <- length(traits)

# Apply names to S matrix if missing
if (is.null(rownames(S))) {
  rownames(S) <- traits
  colnames(S) <- traits
}

cat("[INFO] Found", n_traits, "traits:", paste(traits, collapse = ", "), "\n")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Genetic covariance matrix (S)
# ------------------------------------------------------------------------------
cat("[INFO] Exporting genetic covariance matrix...\n")
S_df <- as.data.frame(S)
S_df <- cbind(trait = rownames(S), S_df)
fwrite(S_df, file.path(outdir, "S_genetic_covariance.tsv"), sep = "\t", quote = FALSE)

# ------------------------------------------------------------------------------
# Genetic correlation matrix (rg)
# ------------------------------------------------------------------------------
cat("[INFO] Exporting genetic correlation matrix...\n")
rg <- cov2cor(S)
rg_df <- as.data.frame(rg)
rg_df <- cbind(trait = rownames(rg), rg_df)
fwrite(rg_df, file.path(outdir, "rg_genetic_correlation.tsv"), sep = "\t", quote = FALSE)

# ------------------------------------------------------------------------------
# Heritabilities
# ------------------------------------------------------------------------------
cat("[INFO] Exporting heritabilities...\n")
h2 <- diag(S)
h2_df <- data.table(trait = traits, h2 = h2)
h2_df <- h2_df[order(-h2)]
fwrite(h2_df, file.path(outdir, "heritabilities.tsv"), sep = "\t", quote = FALSE)

cat("[INFO] Heritabilities:\n")
for (i in 1:nrow(h2_df)) {
  cat(sprintf("  %-20s h² = %.4f\n", h2_df$trait[i], h2_df$h2[i]))
}

# ------------------------------------------------------------------------------
# Intercepts
# ------------------------------------------------------------------------------
if (!is.null(I)) {
  cat("[INFO] Exporting LDSC intercepts...\n")
  
  if (is.matrix(I)) {
    I_diag <- diag(I)
  } else if (is.vector(I) && length(I) >= n_traits) {
    I_diag <- I[1:n_traits]
  } else {
    I_diag <- rep(NA_real_, n_traits)
  }
  
  I_df <- data.table(trait = traits, intercept = I_diag)
  fwrite(I_df, file.path(outdir, "intercepts.tsv"), sep = "\t", quote = FALSE)
  
  cat("[INFO] Intercepts (should be ~1.0):\n")
  for (i in 1:nrow(I_df)) {
    flag <- if (!is.na(I_df$intercept[i]) && I_df$intercept[i] > 1.1) " [elevated]" else ""
    cat(sprintf("  %-20s %.4f%s\n", I_df$trait[i], I_df$intercept[i], flag))
  }
}

# ------------------------------------------------------------------------------
# Hierarchical clustering
# ------------------------------------------------------------------------------
cat("[INFO] Clustering traits by genetic correlation...\n")

rg_clean <- rg
rg_clean[!is.finite(rg_clean)] <- 0
diag(rg_clean) <- 1

# Distance = 1 - |rg| so highly correlated traits cluster together
dist_mat <- as.dist(1 - abs(rg_clean))
hc <- hclust(dist_mat, method = "complete")

cluster_order <- hc$order
traits_ordered <- traits[cluster_order]

cat("[INFO] Cluster order:\n  ", paste(traits_ordered, collapse = " → "), "\n")

# Export clustered matrix
rg_clustered <- rg[cluster_order, cluster_order]
rg_clust_df <- as.data.frame(rg_clustered)
rg_clust_df <- cbind(trait = rownames(rg_clustered), rg_clust_df)
fwrite(rg_clust_df, file.path(outdir, "rg_clustered.tsv"), sep = "\t", quote = FALSE)

# ------------------------------------------------------------------------------
# Pairwise correlations
# ------------------------------------------------------------------------------
cat("[INFO] Extracting pairwise correlations...\n")

rg_pairs <- data.table(
  trait1 = character(),
  trait2 = character(),
  rg = numeric()
)

for (i in 1:(n_traits - 1)) {
  for (j in (i + 1):n_traits) {
    rg_pairs <- rbind(rg_pairs, data.table(
      trait1 = traits[i],
      trait2 = traits[j],
      rg = rg[i, j]
    ))
  }
}

rg_pairs <- rg_pairs[order(-abs(rg))]
fwrite(rg_pairs, file.path(outdir, "pairwise_correlations.tsv"), sep = "\t", quote = FALSE)

cat("\n[INFO] Top 10 genetic correlations:\n")
for (i in 1:min(10, nrow(rg_pairs))) {
  cat(sprintf("  %-18s <-> %-18s rg = %+.3f\n", 
              rg_pairs$trait1[i], rg_pairs$trait2[i], rg_pairs$rg[i]))
}

# ------------------------------------------------------------------------------
# Heatmap visualization
# ------------------------------------------------------------------------------
cat("\n[INFO] Creating heatmap...\n")

if (has_pheatmap) {
  library(pheatmap)
  
  # Blue-white-red palette
  colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                               "white", 
                               "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100)
  
  # PDF
  pdf(file.path(outdir, "rg_heatmap.pdf"), width = 12, height = 10)
  pheatmap(
    rg,
    color = colors,
    breaks = seq(-1, 1, length.out = 101),
    clustering_method = "complete",
    clustering_distance_rows = dist_mat,
    clustering_distance_cols = dist_mat,
    display_numbers = TRUE,
    number_format = "%.2f",
    number_color = "black",
    fontsize_number = 7,
    fontsize = 10,
    main = "Genetic Correlations (hierarchically clustered)",
    angle_col = 45
  )
  dev.off()
  
  # PNG
  png(file.path(outdir, "rg_heatmap.png"), width = 1200, height = 1000, res = 100)
  pheatmap(
    rg,
    color = colors,
    breaks = seq(-1, 1, length.out = 101),
    clustering_method = "complete",
    clustering_distance_rows = dist_mat,
    clustering_distance_cols = dist_mat,
    display_numbers = TRUE,
    number_format = "%.2f",
    number_color = "black",
    fontsize_number = 7,
    fontsize = 10,
    main = "Genetic Correlations (hierarchically clustered)",
    angle_col = 45
  )
  dev.off()
  
  cat("[INFO] Saved: rg_heatmap.pdf and rg_heatmap.png\n")
  
} else if (has_corrplot) {
  library(corrplot)
  
  pdf(file.path(outdir, "rg_heatmap.pdf"), width = 12, height = 10)
  corrplot(
    rg_clustered,
    method = "color",
    type = "full",
    order = "original",
    addCoef.col = "black",
    number.cex = 0.7,
    tl.col = "black",
    tl.srt = 45,
    title = "Genetic Correlations (clustered)",
    mar = c(0, 0, 2, 0)
  )
  dev.off()
  
  cat("[INFO] Saved: rg_heatmap.pdf (corrplot)\n")
  
} else {
  # Base R
  pdf(file.path(outdir, "rg_heatmap.pdf"), width = 12, height = 10)
  heatmap(
    rg,
    Rowv = as.dendrogram(hc),
    Colv = as.dendrogram(hc),
    scale = "none",
    col = colorRampPalette(c("blue", "white", "red"))(100),
    margins = c(10, 10),
    main = "Genetic Correlations"
  )
  dev.off()
  
  cat("[INFO] Saved: rg_heatmap.pdf (base R)\n")
}

# ------------------------------------------------------------------------------
# Dendrogram
# ------------------------------------------------------------------------------
pdf(file.path(outdir, "dendrogram.pdf"), width = 10, height = 6)
par(mar = c(5, 4, 4, 2))
plot(hc, 
     main = "Trait Clustering (based on |rg|)",
     xlab = "", sub = "",
     ylab = "Distance (1 - |rg|)")
dev.off()
cat("[INFO] Saved: dendrogram.pdf\n")

# ------------------------------------------------------------------------------
# Summary report
# ------------------------------------------------------------------------------
sink(file.path(outdir, "summary_report.txt"))
cat("=============================================================================\n")
cat("                    LDSC GENETIC CORRELATION SUMMARY                         \n")
cat("=============================================================================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input:", ldsc_rds, "\n")
cat("Traits:", n_traits, "\n\n")

cat("HERITABILITIES\n")
cat("-----------------------------------------------------------------------------\n")
for (i in 1:nrow(h2_df)) {
  cat(sprintf("  %-20s h² = %.4f\n", h2_df$trait[i], h2_df$h2[i]))
}

cat("\nCLUSTER ORDER (most similar grouped)\n")
cat("-----------------------------------------------------------------------------\n")
cat(paste(traits_ordered, collapse = " → "), "\n")

cat("\nSTRONGEST CORRELATIONS\n")
cat("-----------------------------------------------------------------------------\n")
for (i in 1:min(15, nrow(rg_pairs))) {
  cat(sprintf("  %-18s <-> %-18s  rg = %+.3f\n", 
              rg_pairs$trait1[i], rg_pairs$trait2[i], rg_pairs$rg[i]))
}

cat("\nFILES GENERATED\n")
cat("-----------------------------------------------------------------------------\n")
cat("  S_genetic_covariance.tsv    - Genetic covariance matrix\n")
cat("  rg_genetic_correlation.tsv  - Genetic correlation matrix\n")
cat("  rg_clustered.tsv            - Clustered correlation matrix\n")
cat("  heritabilities.tsv          - h² per trait\n")
cat("  intercepts.tsv              - LDSC intercepts\n")
cat("  pairwise_correlations.tsv   - Pairwise rg (sorted)\n")
cat("  rg_heatmap.pdf/png          - Clustered heatmap\n")
cat("  dendrogram.pdf              - Clustering tree\n")
cat("=============================================================================\n")
sink()

cat("[INFO] Saved: summary_report.txt\n")

# ------------------------------------------------------------------------------
# Done
# ------------------------------------------------------------------------------
cat("\n=============================================================\n")
cat("                  SCRIPT 02b COMPLETE                        \n")
cat("=============================================================\n")
cat("Output:", outdir, "\n")
cat("=============================================================\n")