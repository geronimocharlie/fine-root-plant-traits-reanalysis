############################################################
# ----------------------------------------------------------
# Purpose:
#   Run the FIRST (corrected) Procrustes test:
#     "Does the not-imputed configuration of species occupy
#      the same trait-space 'landscape' as the imputed one?"
#
# Correct approach:
#   1) Load the *reference* PCA object computed on the FULL
#      imputed dataset (already saved by your colleague).
#   2) Load the not-imputed combined trait dataset.
#   3) Standardize not-imputed traits using the *reference*
#      means/sds (critical).
#   4) Project not-imputed data into the *same* PCA space
#      using the reference rotated loadings.
#   5) Match shared species IDs and run Procrustes on PC1–PC2.
#
# Folder convention:
#   - scripts: Code/
#   - outputs: Results/
#     - reference PCA: Results/PCA_results/PCA_imputed_combined_full.rds
#     - not-imputed data: Results/PCA_data/not_imputed_combined.rds
#     - save procrustes output: Results/Procrustes_results/
############################################################

## ----------------------------
## 0) Packages
## ----------------------------
pkgs <- c("ade4")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
library(ade4)

## ----------------------------
## 1) Paths
## ----------------------------
ref_pca_file <- file.path("Results", "PCA_imputed_combined_full.rds")
notimp_file  <- file.path("Data", "not_imputed_combined.rds")

out_dir <- file.path("Results", "Procrustes_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ----------------------------
## 2) Load reference PCA + data
## ----------------------------
if (!file.exists(ref_pca_file)) stop("Missing reference PCA file: ", ref_pca_file)
if (!file.exists(notimp_file))  stop("Missing not-imputed data file: ", notimp_file)

ref <- readRDS(ref_pca_file)
notimp <- readRDS(notimp_file)

# Ensure data.frames for safe handling
ref_traits   <- as.data.frame(ref$traits)
ref_scores   <- as.data.frame(ref$traitsUse)           # reference scores (imputed data)
ref_loadings <- as.matrix(ref$PCA$loadings)            # rotated loadings (traits x components)
ref_means    <- ref$means
ref_sds      <- ref$sds

notimp <- as.data.frame(notimp)

## ----------------------------
## 3) Align trait columns between not-imputed data and reference PCA
## ----------------------------
# We require that not-imputed contains (at least) the traits used in the reference.
# Best practice: align by trait names (column names).

ref_trait_names <- colnames(ref_traits)
if (is.null(ref_trait_names) || any(ref_trait_names == "")) {
  stop("Reference traits have no valid column names. Cannot align traits safely.")
}

if (is.null(colnames(notimp)) || any(colnames(notimp) == "")) {
  stop("Not-imputed data has no valid column names. Cannot align traits safely.")
}

missing_in_notimp <- setdiff(ref_trait_names, colnames(notimp))
if (length(missing_in_notimp) > 0) {
  stop("Not-imputed data is missing trait columns required by reference PCA:\n  - ",
       paste(missing_in_notimp, collapse = "\n  - "))
}

# Keep only reference traits, in the reference order
notimp_aligned <- notimp[, ref_trait_names, drop = FALSE]

# Ensure numeric
non_num <- ref_trait_names[!sapply(notimp_aligned, is.numeric)]
if (length(non_num) > 0) {
  stop("Non-numeric trait columns in not-imputed data (must be numeric for projection): ",
       paste(non_num, collapse = ", "))
}

# Align loadings rows to trait names if possible
if (!is.null(rownames(ref_loadings)) && all(ref_trait_names %in% rownames(ref_loadings))) {
  ref_loadings <- ref_loadings[ref_trait_names, , drop = FALSE]
} else {
  # If loadings have no rownames, we rely on consistent ordering
  # (this is less safe, but your earlier check indicated alignment is OK).
  warning("Reference loadings do not have usable rownames; relying on column order.")
}

## ----------------------------
## 4) Standardize not-imputed data using REFERENCE means/sds
## ----------------------------
# Critical: Do NOT standardize using not-imputed mean/sd.
# We want coordinates in the reference map.

Z_not <- scale(notimp_aligned, center = ref_means, scale = ref_sds)

## ----------------------------
## 5) Project not-imputed species into the reference PCA space
## ----------------------------
# Projection rule for this workflow:
#   scores_new = Z_new %*% Loadings_ref
#
# Here, Loadings_ref are the rotated loadings saved in ref$PCA$loadings.

scores_not <- Z_not %*% ref_loadings
scores_not <- as.data.frame(scores_not)
colnames(scores_not) <- colnames(ref_loadings)

# Keep species IDs
if (is.null(rownames(notimp_aligned)) || any(rownames(notimp_aligned) == "")) {
  stop("Not-imputed data has missing/empty rownames. Need species IDs as rownames.")
}
rownames(scores_not) <- rownames(notimp_aligned)

## ----------------------------
## 6) Match shared species IDs (required for Procrustes)
## ----------------------------
if (is.null(rownames(ref_scores)) || any(rownames(ref_scores) == "")) {
  stop("Reference scores (traitsUse) have missing/empty rownames.")
}

common_ids <- intersect(rownames(ref_scores), rownames(scores_not))
if (length(common_ids) < 10) {
  warning("Very small overlap of species between reference and not-imputed: n = ",
          length(common_ids), "\nInterpretation may be unstable.")
}

# Use PC1–PC4 (first four dimensions) for the standard Procrustes comparison

if (ncol(ref_scores) < 4 || ncol(scores_not) < 4) {
  stop("Cannot use PC1–PC4: reference or projected scores have < 4 dimensions.")
}

axes_to_use <- 1:4
X_ref <- ref_scores[common_ids, axes_to_use, drop = FALSE]
X_not <- scores_not[common_ids, axes_to_use, drop = FALSE]

# Ensure same ordering
X_not <- X_not[rownames(X_ref), , drop = FALSE]
stopifnot(identical(rownames(X_ref), rownames(X_not)))

cat("\n==============================\n")
cat("Procrustes test setup\n")
cat("==============================\n")
cat("Reference PCA file: ", ref_pca_file, "\n", sep = "")
cat("Not-imputed data:   ", notimp_file, "\n", sep = "")
cat("Shared species used:", nrow(X_ref), "\n")
cat("Axes compared:      PC1-PC4\n\n")

## ----------------------------
## 7) Run Procrustes test (permutation)
## ----------------------------
set.seed(1)  # for reproducible permutation test
proc_res <- ade4::procuste.rtest(as.data.frame(X_ref),
                                as.data.frame(X_not),
                                nrepet = 9999)

print(proc_res)

## ----------------------------
## 8) Save outputs
## ----------------------------
out <- list(
  reference_pca_file = ref_pca_file,
  notimputed_file    = notimp_file,
  n_shared_species   = nrow(X_ref),
  axes_used          = c("PC1", "PC2"),
  scores_reference   = X_ref,
  scores_projected   = X_not,
  procrustes         = proc_res
)

saveRDS(out, file.path(out_dir, "procrustes_imputed_ref_vs_notimputed_projected_PC1PC4.rds"))

cat("\nSaved Procrustes output to:\n  ",
    file.path(out_dir, "procrustes_imputed_ref_vs_notimputed_projected_PC1PC4.rds"),
    "\n", sep = "")
