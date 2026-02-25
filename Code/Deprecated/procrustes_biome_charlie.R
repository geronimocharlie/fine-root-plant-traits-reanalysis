############################################################

# What this script does:
#   For each Biome × Domain (Above/Below):
#     1) Load the global TOTAL reference PCA (imputed).
#     2) Extract the correct reference axes:
#          above -> PCs 1:2
#          below -> PCs 3:4
#     3) Fit a biome-only PCA on the biome trait dataset
#        (same pipeline as perform_PCA: scale + psych::principal + varimax
#         + sqrtEigen rescaling).
#     4) Match shared species IDs and run Procrustes.
#     5) Save per-test .rds + one summary .csv
#
# Folder convention:
#   - scripts: Code/
#   - outputs: Results/
############################################################

## ----------------------------
## 0) Packages
## ----------------------------
pkgs <- c("ade4", "psych", "paran")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(ade4)
library(psych)
library(paran)

set.seed(1)

## ----------------------------
## 1) Paths
## ----------------------------
# IMPORTANT: This must be the TOTAL (combined) imputed PCA object
# whose traitsUse contains at least 4 columns (PC1–PC4).
ref_total_file <- file.path("Results", "PCA_imputed_combined_full.rds")

data_files <- list(
  temperate   = list(above = file.path("Data", "temperate_above.rds"),
                     below = file.path("Data", "temperate_below.rds")),
  tropical    = list(above = file.path("Data", "tropical_above.rds"),
                     below = file.path("Data", "tropical_below.rds")),
  continental = list(above = file.path("Data", "continental_above.rds"),
                     below = file.path("Data", "continental_below.rds"))
)

out_dir <- file.path("Results", "Procrustes_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ----------------------------
## 2) Helper: Fit PCA on biome-only dataset (same logic as perform_PCA)
## ----------------------------
fit_biome_pca_like <- function(traits_df) {
  traits_df <- as.data.frame(traits_df)

  if (is.null(rownames(traits_df)) || any(rownames(traits_df) == "")) {
    stop("Biome traits data must have rownames = species IDs.")
  }
  if (is.null(colnames(traits_df)) || any(colnames(traits_df) == "")) {
    stop("Biome traits data must have colnames = trait names.")
  }

  # Drop common non-trait columns if present
  drop_cols <- intersect(colnames(traits_df), c("genus", "family", "order", "biomesKoeppenGroup", "biome"))
  if (length(drop_cols) > 0) {
    traits_df <- traits_df[, setdiff(colnames(traits_df), drop_cols), drop = FALSE]
  }

  non_num <- colnames(traits_df)[!sapply(traits_df, is.numeric)]
  if (length(non_num) > 0) {
    stop("Non-numeric trait columns found in biome dataset: ", paste(non_num, collapse = ", "))
  }

  k <- paran(traits_df, graph = FALSE)$Retained
  if (is.na(k) || k < 1) stop("paran failed to retain a valid number of dimensions.")

  p <- psych::principal(scale(traits_df),
                        nfactors = k,
                        rotate   = "varimax",
                        covar    = FALSE)

  # Same rescaling as perform_PCA
  sqrtEigen <- sqrt(colSums(p$loadings^2))
  for (i in 1:k) {
    p$scores[, i] <- p$scores[, i] * sqrtEigen[i]
  }
  sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(p$loadings)),
                         byrow = TRUE,
                         nrow  = nrow(p$loadings))
  p$loadings <- p$loadings / sqrtEigenMat

  traitsUse <- as.data.frame(p$scores[, 1:k, drop = FALSE])
  rownames(traitsUse) <- rownames(traits_df)

  list(dimensions = k, traitsUse = traitsUse)
}

## ----------------------------
## 3) Load TOTAL reference PCA (global imputed)
## ----------------------------
if (!file.exists(ref_total_file)) stop("Missing TOTAL reference PCA file: ", ref_total_file)

ref_total <- readRDS(ref_total_file)

if (!("traitsUse" %in% names(ref_total))) stop("Reference PCA missing traitsUse: ", ref_total_file)
ref_total_scores <- as.data.frame(ref_total$traitsUse)

if (ncol(ref_total_scores) < 4) {
  stop("Reference TOTAL PCA has <4 axes in traitsUse. Need PC1–PC4 for above/below mapping.")
}
if (is.null(rownames(ref_total_scores)) || any(rownames(ref_total_scores) == "")) {
  stop("Reference TOTAL traitsUse must have rownames (species IDs).")
}

# Domain-to-axis mapping (THE FIX)
axis_map <- list(
  above = 1:2,  # PC1–PC2
  below = 3:4   # PC3–PC4
)

## ----------------------------
## 4) Main loop: 3 biomes × 2 domains = 6 tests
## ----------------------------
summary_rows <- list()

for (biome in names(data_files)) {
  for (domain in c("above", "below")) {

    cat("\n====================================================\n")
    cat("Biome test:", biome, "| Domain:", domain, "\n")
    cat("====================================================\n")

    biome_path <- data_files[[biome]][[domain]]
    if (!file.exists(biome_path)) stop("Missing biome dataset: ", biome_path)

    # ---- Fit biome-only PCA
    biome_traits <- readRDS(biome_path)
    biome_pca <- fit_biome_pca_like(biome_traits)
    biome_scores <- as.data.frame(biome_pca$traitsUse)

    # ---- Choose correct reference axes
    axes_ref <- axis_map[[domain]]
    ref_dom_scores <- ref_total_scores[, axes_ref, drop = FALSE]

    # ---- Biome PCA must provide at least 2 axes for Procrustes
    if (ncol(biome_scores) < 2) {
      stop("Biome PCA retained <2 axes for ", biome, " / ", domain, ". Cannot run Procrustes.")
    }

    # Use first two biome axes (PC1–PC2 of the biome-only PCA)
    axes_biome <- 1:2
    biome_scores_2d <- biome_scores[, axes_biome, drop = FALSE]

    # ---- Match species IDs
    common_ids <- intersect(rownames(ref_dom_scores), rownames(biome_scores_2d))
    if (length(common_ids) < 10) {
      warning("Very small overlap of species for ", biome, " / ", domain, ": n = ", length(common_ids),
              "\nInterpretation may be unstable.")
    }

    X_ref   <- ref_dom_scores[common_ids, , drop = FALSE]
    X_biome <- biome_scores_2d[common_ids, , drop = FALSE]

    # Ensure identical ordering
    X_biome <- X_biome[rownames(X_ref), , drop = FALSE]
    stopifnot(identical(rownames(X_ref), rownames(X_biome)))

    cat("Shared species used:", nrow(X_ref), "\n")
    cat("Reference axes:      PC", min(axes_ref), "–PC", max(axes_ref), "\n", sep = "")
    cat("Biome axes:          PC1–PC2 (biome-only PCA)\n")
    cat("Biome PCA dims (paran):", biome_pca$dimensions, "\n")

    # ---- Procrustes test
    set.seed(1)
    proc_res <- ade4::procuste.rtest(as.data.frame(X_ref),
                                    as.data.frame(X_biome),
                                    nrepet = 9999)

    print(proc_res)

    # ---- Save per-test output
    out <- list(
      biome = biome,
      domain = domain,
      biome_data_file = biome_path,
      reference_pca_file = ref_total_file,
      reference_axes = paste0("PC", axes_ref),
      biome_axes = c("PC1", "PC2"),
      n_shared_species = nrow(X_ref),
      biome_pca_dimensions = biome_pca$dimensions,
      reference_scores = X_ref,
      biome_scores = X_biome,
      procrustes = proc_res
    )

    out_file <- file.path(out_dir, paste0("biome_procrustes_", biome, "_", domain, "_FIXEDAXES.rds"))
    saveRDS(out, out_file)
    cat("Saved:", out_file, "\n")

    # ---- Summary row
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      biome = biome,
      domain = domain,
      n_shared_species = nrow(X_ref),
      reference_axes = paste0("PC", min(axes_ref), "-PC", max(axes_ref)),
      biome_axes = "PC1-PC2",
      obs = unname(proc_res$obs),
      p_value = unname(proc_res$pvalue),
      stringsAsFactors = FALSE
    )
  }
}

## ----------------------------
## 5) Write summary table
## ----------------------------
summary_df <- do.call(rbind, summary_rows)
summary_csv <- file.path(out_dir, "biome_procrustes_summary_FIXEDAXES.csv")
write.csv(summary_df, summary_csv, row.names = FALSE)

cat("\n====================================================\n")
cat("ALL DONE.\n")
cat("Summary written to:\n  ", summary_csv, "\n", sep = "")
cat("Per-test .rds files are in:\n  ", out_dir, "\n", sep = "")
cat("====================================================\n")
