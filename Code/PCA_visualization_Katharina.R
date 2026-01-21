library(ggplot2)
library(MASS)

setwd("C:/Users/katha/OneDrive/2_documents/4_Studium/EMDS_Semester3/Capstone")
getwd()

# PCA plot function (simple version)
plotPCA <- function(PCA_dataset, title) {
  # Extract scores and loadings
  scores   <- as.data.frame(PCA_dataset$PCA$scores[, 1:2])
  loadings <- as.data.frame(unclass(PCA_dataset$PCA$loadings[, 1:2]))
  
  # Variance explained
  var_exp <- PCA_dataset$Variance[1:2] * 100
  
  scores$species <- rownames(scores)
  loadings$trait <- rownames(loadings)
  
  scale_factor <- max(abs(scores[, 1:2]))
  
  ggplot(scores, aes(RC1, RC2)) +
    geom_point(color = "grey30", size = 2) +
    geom_segment(
      data = loadings,
      aes(
        x = 0, y = 0,
        xend = RC1 * scale_factor,
        yend = RC2 * scale_factor
      ),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "red"
    ) +
    geom_text(
      data = loadings,
      aes(
        x = RC1 * scale_factor * 1.1,
        y = RC2 * scale_factor * 1.1,
        label = trait
      ),
      color = "red",
      size = 3
    ) +
    labs(
      title = title,
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)")
    ) +
    theme_minimal()
}


compute_kde_prob <- function(x, y, n = 200) {
  kd <- kde2d(x, y, n = n)
  
  # Flatten density
  z <- as.vector(kd$z)
  
  # Convert density to cumulative probability
  prob <- z / sum(z)
  ord <- order(prob, decreasing = TRUE)
  cumprob <- cumsum(prob[ord])
  
  z_prob <- numeric(length(z))
  z_prob[ord] <- cumprob
  
  list(
    x = kd$x,
    y = kd$y,
    z = matrix(z_prob, nrow = n)
  )
}

# PCA plot with species density quantiles
plotFunctionalSpace_prob <- function(PCA_dataset,
                                     probs = c(0.5, 0.75, 0.95, 0.99),
                                     title = "") {
  
  scores <- as.data.frame(PCA_dataset$PCA$scores[, 1:2])
  colnames(scores) <- c("C1", "C2")
  
  var_exp <- PCA_dataset$Variance[1:2] * 100
  
  # KDE with probability mass
  kde <- compute_kde_prob(scores$C1, scores$C2)
  
  kde_df <- expand.grid(
    C1 = kde$x,
    C2 = kde$y
  )
  kde_df$prob <- as.vector(kde$z)
  
  # Loadings
  loadings <- as.data.frame(unclass(PCA_dataset$PCA$loadings[, 1:2]))
  colnames(loadings) <- c("C1", "C2")
  loadings$trait <- rownames(loadings)
  
  scale_factor <- max(abs(scores))
  
  ggplot() +
    # --- Filled probability regions ---
    geom_raster(
      data = kde_df,
      aes(C1, C2, fill = prob),
      interpolate = TRUE,
      alpha = 0.9
    ) +
    
    # --- Probability contours ---
    geom_contour(
      data = kde_df,
      aes(C1, C2, z = prob),
      breaks = probs,
      color = "black",
      linewidth = 0.6
    ) +
    
    # --- Species points ---
    geom_point(
      data = scores,
      aes(C1, C2),
      size = 0.7,
      alpha = 0.7,
      color = "grey30"
    ) +
    
    # --- Trait arrows ---
    geom_segment(
      data = loadings,
      aes(
        x = 0, y = 0,
        xend = C1 * scale_factor,
        yend = C2 * scale_factor
      ),
      arrow = arrow(length = unit(0.2, "cm")),
      linewidth = 0.7
    ) +
    geom_text(
      data = loadings,
      aes(
        x = C1 * scale_factor * 1.1,
        y = C2 * scale_factor * 1.1,
        label = trait
      ),
      size = 3
    ) +
    scale_fill_distiller(
      palette = "YlOrRd",
      direction = -1,
      name = "Cumulative\nprobability"
    ) +
    
    labs(
      title = title,
      x = paste0("C1 (", round(var_exp[1], 1), "%)"),
      y = paste0("C2 (", round(var_exp[2], 1), "%)")
    ) +
    
    coord_equal() +
    theme_minimal()
}


# export table function
exportPCA_table <- function(PCA_dataset,
                            n_axes = 4,
                            rotation = c("none", "varimax"),
                            digits = 2) {
  
  rotation <- match.arg(rotation)
  
  # --- Loadings ---
  loadings <- PCA_dataset$PCA$loadings[, 1:n_axes, drop = FALSE]
  
  if (rotation == "varimax") {
    vm <- varimax(loadings)
    loadings <- vm$loadings
  }
  
  loadings <- round(unclass(loadings), digits)
  
  # --- Variance ---
  variance <- PCA_dataset$Variance[1:n_axes] * 100
  
  # --- Eigenvalues (derived, robust) ---
  n_vars <- nrow(PCA_dataset$PCA$loadings)
  eigenvalues <- PCA_dataset$Variance[1:n_axes] * n_vars
  
  eigen_table <- rbind(
    eigenvalue   = round(eigenvalues, digits),
    `% variance` = round(variance, digits)
  )
  
  # --- Final table ---
  table_out <- rbind(
    eigen_table,
    loadings
  )
  
  return(as.data.frame(table_out))
}

# Pipeline
fun_pipeline <- function(file_name, plot_name, table_file_name, rotation="none"){
  PCA_file <- readRDS(paste0("./PCA_results/", file_name, ".rds"))
  plot <- plotFunctionalSpace_prob(PCA_dataset = PCA_file,
                                   probs = c(0.50, 0.95, 0.99),
                                   title = plot_name)
  ggsave(
    filename = paste0("./PCA_plots/", plot_name, ".png"),
    plot = plot,
    width = 18,
    height = 15,
    units = "cm",
    dpi = 300
  )
  pca_table <- exportPCA_table(PCA_dataset = PCA_file, n_axes = 2, rotation = rotation)
  write.csv(pca_table, paste0("./PCA_plots/", table_file_name, ".csv"))}


# make plots
fun_pipeline("PCA_imputed_combined_full", "imputed_combined_plot", "imputed_combined")
fun_pipeline("PCA_imputed_combined_half" , "imputed_combined_half_plot", "imputed_combined_half")
fun_pipeline("PCA_imputed_combined_quarter", "imputed_combined_quarter_plot", "imputed_combined_quarter")
fun_pipeline("PCA_imputed_combined_tenth", "imputed_combined_tenth_plot", "imputed_combined_tenth")

fun_pipeline("PCA_imputed_below", "imputed_below_plot", "imputed_below")
fun_pipeline("PCA_imputed_above", "imputed_above_plot", "imputed_above")

fun_pipeline("PCA_below", "below_plot", "below")
fun_pipeline("PCA_above", "above_plot", "above")
fun_pipeline("PCA_combined", "combined_plot", "combined")

fun_pipeline("PCA_imputed_continental_above", "imputed_continental_above_plot", "imputed_continental_above")
fun_pipeline("PCA_imputed_continental_below", "imputed_continental_below_plot", "imputed_continental_below")
fun_pipeline("PCA_imputed_temperate_above", "imputed_temperate_above_plot", "imputed_temperate_above")
fun_pipeline("PCA_imputed_temperate_below", "imputed_temperate_below_plot", "imputed_temperate_below")
fun_pipeline("PCA_imputed_tropical_above", "imputed_tropical_above_plot", "imputed_tropical_above")
fun_pipeline("PCA_imputed_tropical_below", "imputed_tropical_below_plot", "imputed_tropical_below")

fun_pipeline("PCA_continental_above", "continental_above_plot", "continental_above")
fun_pipeline("PCA_continental_below", "continental_below_plot", "continental_below")
fun_pipeline("PCA_temperate_above", "temperate_above_plot", "temperate_above")
fun_pipeline("PCA_temperate_below", "temperate_below_plot", "temperate_below")
fun_pipeline("PCA_tropical_above", "tropical_above_plot", "tropical_above")
fun_pipeline("PCA_tropical_below", "tropical_below_plot", "tropical_below")
