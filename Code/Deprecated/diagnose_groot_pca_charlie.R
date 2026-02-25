diagnose_pca_dims <- function(pca_path) {
  cat("\n==============================\n")
  cat("Diagnosing:", pca_path, "\n")
  cat("==============================\n")

  obj <- readRDS(pca_path)

  # --- What your script uses
  if (!is.null(obj$traitsUse)) {
    cat("traitsUse columns (your k_biome): ", ncol(obj$traitsUse), "\n", sep = "")
    cat("traitsUse colnames: ", paste(colnames(obj$traitsUse), collapse = ", "), "\n", sep = "")
  } else {
    cat("WARNING: obj$traitsUse is NULL\n")
  }

  # --- What paran retained (if stored)
  if (!is.null(obj$dimensions)) {
    cat("obj$dimensions (paran retained): ", obj$dimensions, "\n", sep = "")
  } else {
    cat("NOTE: obj$dimensions not stored\n")
  }

  # --- Loadings matrix (what colleagues often look at)
  if (!is.null(obj$PCA) && !is.null(obj$PCA$loadings)) {
    L <- as.matrix(obj$PCA$loadings)
    cat("loadings matrix dims: ", nrow(L), " x ", ncol(L), "\n", sep = "")
    cat("loading colnames: ", paste(colnames(L), collapse = ", "), "\n", sep = "")
  } else {
    cat("WARNING: obj$PCA$loadings not found\n")
  }

  # --- Variance explained (if present)
  if (!is.null(obj$Variance)) {
    cat("Variance vector length: ", length(obj$Variance), "\n", sep = "")
    cat("Variance values: ", paste(round(obj$Variance, 4), collapse = ", "), "\n", sep = "")
  } else if (!is.null(obj$PCA) && !is.null(obj$PCA$Vaccounted)) {
    cat("PCA$Vaccounted available. (printing proportions row if possible)\n")
    print(obj$PCA$Vaccounted)
  } else {
    cat("NOTE: variance info not found\n")
  }

  # --- “Effective dimension” check: are later axes basically zero?
  # For rotated factors, a quick heuristic is the column SD of scores:
  if (!is.null(obj$traitsUse)) {
    sds <- apply(obj$traitsUse, 2, sd, na.rm = TRUE)
    cat("SD of scores per axis: ", paste(round(sds, 4), collapse = ", "), "\n", sep = "")
    cat("Heuristic: if last axis SD is ~0, it’s effectively collapsed.\n")
  }

  invisible(TRUE)
}

# Example:
diagnose_pca_dims("Results/PCA_imputed_continental_groot.rds")



biome_files <- c(
  temperate   = "Results/PCA_imputed_temperate_groot.rds",
  tropical    = "Results/PCA_imputed_tropical_groot.rds",
  continental = "Results/PCA_imputed_continental_groot.rds"
)

dims_list <- lapply(names(biome_files), function(b) {
  obj <- readRDS(biome_files[[b]])

  k_traitsUse <- if (!is.null(obj$traitsUse)) ncol(obj$traitsUse) else NA_integer_
  k_dim       <- if (!is.null(obj$dimensions)) obj$dimensions else NA_integer_
  k_loadings  <- if (!is.null(obj$PCA) && !is.null(obj$PCA$loadings)) ncol(as.matrix(obj$PCA$loadings)) else NA_integer_

  sds <- if (!is.null(obj$traitsUse)) apply(obj$traitsUse, 2, sd, na.rm = TRUE) else NA
  last_sd <- if (!all(is.na(sds))) tail(sds, 1) else NA_real_

  data.frame(
    biome = b,
    k_dimensions_paran = k_dim,
    k_traitsUse = k_traitsUse,
    k_loadings = k_loadings,
    sd_last_axis = last_sd,
    stringsAsFactors = FALSE
  )
})

dims_tbl <- do.call(rbind, dims_list)

print(dims_tbl)
write.csv(dims_tbl,
          "Results/Procrustes_results/biome_pca_dimension_diagnostics.csv",
          row.names = FALSE)


check_object_identity <- function(pca_path) {
  obj <- readRDS(pca_path)
  cat("\nFile:", pca_path, "\n")
  if (!is.null(obj$traits)) cat("traits dim:", dim(obj$traits), "\n")
  if (!is.null(obj$PCA) && !is.null(obj$PCA$n.obs)) cat("n.obs:", obj$PCA$n.obs, "\n")
  if (!is.null(obj$PCA) && !is.null(obj$PCA$call)) {
    cat("PCA call:\n")
    print(obj$PCA$call)
  }
}
check_object_identity("Results/PCA_imputed_tropical_groot.rds")


plot_variance_profile <- function(pca_path) {
  obj <- readRDS(pca_path)
  if (is.null(obj$Variance)) stop("No obj$Variance stored in: ", pca_path)

  v <- as.numeric(obj$Variance)
  df <- data.frame(axis = seq_along(v), variance = v)

  ggplot(df, aes(x = axis, y = variance)) +
    geom_col() +
    labs(
      title = basename(pca_path),
      x = "Axis",
      y = "Variance explained (as stored)"
    ) +
    theme_minimal()
}

# Example:
plot_variance_profile("Results/PCA_imputed_tropical_groot.rds")

obj <- readRDS("Results/PCA_imputed_continental_groot.rds")
cat(
  dimensions = obj$dimensions,
  k_traitsUse = ncol(obj$traitsUse),
  k_loadings  = ncol(as.matrix(obj$PCA$loadings))
)
