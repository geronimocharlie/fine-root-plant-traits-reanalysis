library(ggplot2)
library(MASS)

setwd("C:/Users/katha/OneDrive/2_documents/4_Studium/EMDS_Semester3/Capstone")
getwd()

# PCA plot function (simple version)
plotPCA <- function(PCA_dataset, 
                    axes = c(1,2), 
                    probs = c(0.5, 0.75, 0.95, 0.99), 
                    title,
                    show_traits = NULL) {
  
  # Extract scores and loadings
  scores   <- as.data.frame(PCA_dataset$PCA$scores[, axes])
  colnames(scores) <- c("C1", "C2")
  
  # Variance explained
  var_exp <- PCA_dataset$Variance[axes] * 100
  
  # KDE with probability mass
  kde <- compute_kde_prob(scores$C1, scores$C2)
  kde_df <- expand.grid(
    C1 = kde$x,
    C2 = kde$y
  )
  kde_df$prob <- as.vector(kde$z)
  
  # Loadings
  loadings <- as.data.frame(unclass(PCA_dataset$PCA$loadings[, axes]))
  colnames(loadings) <- c("C1", "C2")
  loadings$trait <- rownames(loadings)
  
  # Filter traits if requested
  if(!is.null(show_traits)) {
    loadings <- loadings[loadings$trait %in% show_traits, ]
  }
  
  scale_factor <- max(abs(scores))
  
  ggplot() +
    geom_raster(
      data = kde_df,
      aes(C1, C2, fill = prob),
      interpolate = TRUE,
      alpha = 0.9
    ) +
    geom_contour(
      data = kde_df,
      aes(C1, C2, z = prob),
      breaks = probs,
      color = "black",
      linewidth = 0.6
    ) +
    geom_point(
      data = scores,
      aes(C1, C2),
      size = 0.7,
      alpha = 0.7,
      color = "grey30"
    ) +
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
      x = paste0("PC", axes[1], " (", round(var_exp[1],1), "%)"),
      y = paste0("PC", axes[2], " (", round(var_exp[2],1), "%)")
    ) +
    coord_equal() +
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

# export table function
exportPCA_table <- function(PCA_dataset,
                            rotation = "varimax",
                            digits = 2) {
  
  n_axes <- min(2, ncol(PCA_dataset$PCA$loadings))
  
  # --- Loadings ---
  loadings <- PCA_dataset$PCA$loadings[, 1:n_axes, drop = FALSE]
  
  vm <- varimax(loadings)
    loadings <- vm$loadings
  
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
fun_pipeline <- function(file_name, plot_name, table_file_name, rotation="varimax"){
  PCA_file <- readRDS(paste0("./PCA_results/", file_name, ".rds"))
  plot_12 <- plotPCA(PCA_dataset = PCA_file,
                                           axes = c(1,2),
                                           probs = c(0.50, 0.95, 0.99),
                                           title = paste0(plot_name, " (PC1 vs PC2)"),
                     c("ph", "sla", "ln", "sm", "la", "ssd"))
  ggsave(
    filename = paste0("./PCA_plots/", plot_name, "_PC1_PC2.png"),
    plot = plot_12,
    width = 18,
    height = 15,
    units = "cm",
    dpi = 300
  )
  
  if(PCA_file$dimensions >= 4) {
    plot_34 <- plotPCA(PCA_dataset = PCA_file,
                                             axes = c(3,4),
                                             probs = c(0.50, 0.95, 0.99),
                                             title = paste0(plot_name, " (PC3 vs PC4)"),
                       c("SRL", "D", "RTD", "N"))
    ggsave(
      filename = paste0("./PCA_plots/", plot_name, "_PC3_PC4.png"),
      plot = plot_34,
      width = 18,
      height = 15,
      units = "cm",
      dpi = 300
    )
  }
  
  pca_table <- exportPCA_table(PCA_dataset = PCA_file, rotation = rotation)
  write.csv(pca_table, paste0("./PCA_plots/", table_file_name, ".csv"))}


# make plots
fun_pipeline("PCA_imputed_combined", "imputed_combined_plot", "imputed_combined")
fun_pipeline("PCA_imputed_combined_half" , "imputed_combined_half_plot", "imputed_combined_half")
fun_pipeline("PCA_imputed_combined_quarter", "imputed_combined_quarter_plot", "imputed_combined_quarter")
fun_pipeline("PCA_imputed_combined_tenth", "imputed_combined_tenth_plot", "imputed_combined_tenth")

fun_pipeline("PCA_imputed_continental", "imputed_continental_plot", "imputed_continental")
fun_pipeline("PCA_imputed_temperate", "imputed_temperate_plot", "imputed_temperate")
fun_pipeline("PCA_imputed_tropical", "imputed_tropical_plot", "imputed_tropical")
