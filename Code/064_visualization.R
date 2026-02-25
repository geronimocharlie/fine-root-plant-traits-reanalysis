# ----- Functional Space Visualization: Compare Imputation Methods -----
# Creates facet grid plots comparing all 4 imputation methods:
#   1. Original missForest (single)
#   2. Mean imputation (baseline)
#   3. 1x OOB Bootstrap (single)
#   4. 2x OOB Bootstrap (single)
#
# OUTPUT:
#   Results/Figures/Functional_Space_PC1_PC2_comparison.png
#   Results/Figures/Functional_Space_PC3_PC4_comparison.png
#
# Prerequisite: scripts 06, 061, 062, 063 have been run

library(ggplot2)
library(MASS)
library(patchwork)  # for combining plots

# -----------------------------------------------------------------------------
# Load PCATotal objects
# -----------------------------------------------------------------------------
PCA_original <- readRDS("Results/imputation/PCATotal_ImputedObs.rds")
PCA_mean     <- readRDS("Results/imputation/PCATotal_mean_imputation.rds")

# For bootstrap methods: load aggregated mean scores/loadings and build pseudo-PCATotal
# compute mean variance across all bootstraps

# Helper function to compute mean and sd of variance across all bootstrap PCATotal objects
compute_variance_stats <- function(outdir) {
  files <- list.files(outdir, pattern = "^PCATotal_ImputedObs_boot_\\d{3}\\.rds$", full.names = TRUE)
  variances <- lapply(files, function(f) {
    obj <- readRDS(f)
    obj$Variance[1:4]
  })
  variance_matrix <- do.call(rbind, variances)  # 50 × 4
  list(
    mean = colMeans(variance_matrix),
    sd   = apply(variance_matrix, 2, sd)
  )
}

# 1x OOB
scores_1x   <- readRDS("Results/imputation/imputed_bootstrap/PCA_scores_boot_summary.rds") # contains mean and sd of the species position over n bootstraps (this is going to be the scatter points & density)
loadings_1x <- readRDS("Results/imputation/imputed_bootstrap/PCA_loadings_boot_summary.rds") # contains mean and sd of the trait loadings over n bootstraps + the single bootstrap loadings (this is going to be the arrows + labels)
variance_1x <- compute_variance_stats("Results/imputation/imputed_bootstrap")  # Mean AND sd of variance across n bootstraps
PCA_1x <- list(
  PCA = list(
    scores = scores_1x$mean,
    loadings = loadings_1x$mean
  ),
  Variance = variance_1x$mean,
  Variance_sd = variance_1x$sd
)

# 2x OOB
scores_2x   <- readRDS("Results/imputation/imputed_bootstrap_2x/PCA_scores_boot_summary.rds")
loadings_2x <- readRDS("Results/imputation/imputed_bootstrap_2x/PCA_loadings_boot_summary.rds")
variance_2x <- compute_variance_stats("Results/imputation/imputed_bootstrap_2x")  # Mean AND sd of variance across n bootstraps
PCA_2x <- list(
  PCA = list(
    scores = scores_2x$mean,
    loadings = loadings_2x$mean
  ),
  Variance = variance_2x$mean,
  Variance_sd = variance_2x$sd
)

cat("Loaded PCA objects:\n")
cat("  - PCA_original:", nrow(PCA_original$PCA$scores), "species\n")
cat("  - PCA_mean:", nrow(PCA_mean$PCA$scores), "species\n")
cat("  - PCA_1x:", nrow(PCA_1x$PCA$scores), "species\n")
cat("  - PCA_2x:", nrow(PCA_2x$PCA$scores), "species\n")

# -----------------------------------------------------------------------------
# Helper: Compute KDE with cumulative probability
# -----------------------------------------------------------------------------
compute_kde_prob <- function(x, y, n = 200) {
  kd <- kde2d(x, y, n = n)
  z <- as.vector(kd$z)
  prob <- z / sum(z)
  ord <- order(prob, decreasing = TRUE)
  cumprob <- cumsum(prob[ord])
  z_prob <- numeric(length(z))
  z_prob[ord] <- cumprob
  list(x = kd$x, y = kd$y, z = matrix(z_prob, nrow = n))
}

# -----------------------------------------------------------------------------
# Plot function for single PCA dataset (returns ggplot object)
# -----------------------------------------------------------------------------
plotFunctionalSpace_prob <- function(PCA_dataset,
                                     pc_x = 1,
                                     pc_y = 2,
                                     probs = c(0.5, 0.75, 0.95, 0.99),
                                     title = "") {
  
  # Extract scores for selected PCs
  scores <- as.data.frame(PCA_dataset$PCA$scores[, c(pc_x, pc_y)])
  colnames(scores) <- c("C1", "C2")

  # Extract trait loadings for selected PCs
  loadings <- as.data.frame(PCA_dataset$PCA$loadings[, c(pc_x, pc_y), drop = FALSE])
  colnames(loadings) <- c("C1", "C2")
  loadings$trait <- rownames(PCA_dataset$PCA$loadings)

  # Above- vs. below-ground traits (order as in traitsSelect)
  above_traits <- c("la", "ln", "ph", "sla", "ssd", "sm")
  below_traits <- c("SRL", "D", "RTD", "N")

  # For PC1/PC2: show only above-ground trait arrows; for PC3/PC4: only below-ground
  if (pc_x %in% c(1, 2) && pc_y %in% c(1, 2)) {
    loadings <- loadings[loadings$trait %in% above_traits, ]
  } else if (pc_x %in% c(3, 4) && pc_y %in% c(3, 4)) {
    loadings <- loadings[loadings$trait %in% below_traits, ]
  }
  
  # Variance explained
  var_exp <- PCA_dataset$Variance[c(pc_x, pc_y)] * 100
  var_sd  <- PCA_dataset$Variance_sd
  has_sd  <- !is.null(var_sd) && length(var_sd) >= max(pc_x, pc_y) && all(!is.na(var_sd[c(pc_x, pc_y)]))

  # Compute KDE with cumulative probability
 kde <- compute_kde_prob(scores$C1, scores$C2)
  kde_df <- expand.grid(C1 = kde$x, C2 = kde$y)
  kde_df$prob <- as.vector(kde$z)

  # Scale factor for loading arrows
  scale_factor <- min(
    diff(range(scores$C1)),
    diff(range(scores$C2))
  ) / max(abs(loadings$C1), abs(loadings$C2)) * 0.4

  # Build axis labels
  xlab <- if (has_sd) {
    paste0("PC", pc_x, " (", round(var_exp[1], 1), " ± ", round(var_sd[pc_x] * 100, 1), "%)")
  } else {
    paste0("PC", pc_x, " (", round(var_exp[1], 1), "%)")
  }
  ylab <- if (has_sd) {
    paste0("PC", pc_y, " (", round(var_exp[2], 1), " ± ", round(var_sd[pc_y] * 100, 1), "%)")
  } else {
    paste0("PC", pc_y, " (", round(var_exp[2], 1), "%)")
  }

  # Build ggplot
  ggplot() +
    geom_raster(
      data = kde_df,
      aes(C1, C2, fill = prob),
      interpolate = TRUE,
      alpha = 0.9
    ) +
    # Probability contours
    geom_contour(
      data = kde_df,
      aes(C1, C2, z = prob),
      breaks = probs,
      color = "black",
      linewidth = 0.5
    ) +
    # Species points
    geom_point(
      data = scores,
      aes(C1, C2),
      size = 0.5,
      alpha = 0.5,
      color = "grey30"
    ) +
    # Trait arrows
    geom_segment(
      data = loadings,
      aes(x = 0, y = 0,
          xend = C1 * scale_factor,
          yend = C2 * scale_factor),
      arrow = arrow(length = unit(0.15, "cm")),
      linewidth = 0.6,
      color = "black"
    ) +
    # Trait labels
    geom_text(
      data = loadings,
      aes(x = C1 * scale_factor * 1.15,
          y = C2 * scale_factor * 1.15,
          label = trait),
      size = 2.5,
      color = "black"
    ) +
    scale_fill_distiller(
      palette = "YlOrRd",
      direction = -1,
      name = "Cumulative\nprobability",
      limits = c(0, 1)
    ) +
    labs(title = title, x = xlab, y = ylab) +
    coord_equal() +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 9),
      legend.position = "none"  # Will add shared legend later
    )
}

# -----------------------------------------------------------------------------
# Create comparison plots
# -----------------------------------------------------------------------------
cat("\nGenerating Functional Space comparison plots...\n")

# List of PCA objects with labels
pca_list <- list(
  "Original (missForest)" = PCA_original,
  "Mean Imputation"       = PCA_mean,
  "1x OOB Bootstrap"      = PCA_1x,
  "2x OOB Bootstrap"      = PCA_2x
)

# -----------------------------------------------------------------------------
# Save plots (each method as separate figure)
# -----------------------------------------------------------------------------
dir.create("Figures", showWarnings = FALSE)

for (nm in names(pca_list)) {
  obj <- pca_list[[nm]]
  file_stub <- gsub("[^A-Za-z0-9]+", "_", nm)

  # PC1 vs PC2
  p12 <- plotFunctionalSpace_prob(obj, pc_x = 1, pc_y = 2,
                                  probs = c(0.50, 0.95, 0.99),
                                  title = nm)
  file_12 <- sprintf("Results/Figures/Functional_Space_%s_PC1_PC2.png", file_stub)
  ggsave(
    filename = file_12,
    plot = p12,
    width = 12, height = 10, units = "cm", dpi = 300
  )
  cat("Saved:", file_12, "\n")

  # PC3 vs PC4
  p34 <- plotFunctionalSpace_prob(obj, pc_x = 3, pc_y = 4,
                                  probs = c(0.50, 0.95, 0.99),
                                  title = nm)
  file_34 <- sprintf("Results/Figures/Functional_Space_%s_PC3_PC4.png", file_stub)
  ggsave(
    filename = file_34,
    plot = p34,
    width = 12, height = 10, units = "cm", dpi = 300
  )
  cat("Saved:", file_34, "\n")
}

cat("\nVisualization complete!\n")
