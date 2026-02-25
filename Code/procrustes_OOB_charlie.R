############################################################

# ----------------------------------------------------------
# Goal:
#   Quantify whether imputation uncertainty (OOB bootstrap) changes
#   the underlying trait "landscape" (ordination geometry).
#
# What we do (correct concept):
#   1) Fix ONE reference PCA space using:
#        OOB_imputation/PCATotal_mean_Imputation.rds
#      (assumed to be a PCA object that contains traitsUse scores).
#   2) For each bootstrap-imputed TRAIT dataset (50 files):
#        - Standardize using the reference PCA's means/sds
#        - Project into the fixed reference PCA space using reference loadings
#        - Run Procrustes vs the reference configuration (same species)
#   3) Repeat for the 2x OOB-error bootstrap set (50 files).
#   4) Save:
#        - per-iteration results (r, p, n species, axes used)
#        - summary csv
#        - full RDS with details
#
# Assumptions you stated:
#   - Everything for this test is inside folder: OOB_imputation/
#   - Subfolders:
#       OOB_imputation/imputed_bootstrap/      (50 .rds files)
#       OOB_imputation/imputed_bootstrap_2x/   (50 .rds files)
#   - Reference PCA objects:
#       OOB_imputation/PCATotal_mean_Imputation.rds
#       (Optional: OOB_imputation/PCATotal_ImputedObs.rds not required)
#   - Bootstrap files are TRAIT dataframes:
#       rows = species IDs, cols = traits (numeric), same trait names/order
#
# Folder convention:
#   - scripts: Code/
#   - outputs: Results/
############################################################

## ----------------------------
## 0) Packages
## ----------------------------
pkgs <- c("ade4")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
library(ade4)

set.seed(1)

## ----------------------------
## 1) Paths
## ----------------------------
oob_dir <- "OOB_imputation"

ref_file <- file.path(oob_dir, "PCATotal_mean_imputation.rds")

boot_dirs <- list(
  oob_1x = file.path(oob_dir, "imputed_bootstrap"),
  oob_2x = file.path(oob_dir, "imputed_bootstrap_2x")
)

# Where to save outputs
out_dir <- file.path("Results", "Procrustes_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ----------------------------
## 2) One-liner to test a bootstrap file (you asked for this)
## ----------------------------
# Run this ONCE interactively to confirm file type:
# tmp <- readRDS(list.files(file.path(oob_dir,"imputed_bootstrap"), full.names=TRUE)[1]); str(tmp)

## ----------------------------
## 3) Load reference PCA object (fixed map)
## ----------------------------
if (!file.exists(ref_file)) stop("Missing reference file: ", ref_file)
ref <- readRDS(ref_file)

# We expect the same structure you validated earlier:
required_top <- c("traitsUse", "means", "sds", "PCA")
missing_top <- setdiff(required_top, names(ref))
if (length(missing_top) > 0) {
  stop("Reference object is missing required elements:\n  - ",
       paste(missing_top, collapse = "\n  - "))
}

ref_scores_all <- as.data.frame(ref$traitsUse)
ref_means <- ref$means
ref_sds   <- ref$sds
ref_loadings <- as.matrix(ref$PCA$loadings)

if (is.null(rownames(ref_scores_all)) || any(rownames(ref_scores_all) == "")) {
  stop("Reference traitsUse must have rownames (species IDs).")
}

# Decide axes to use: PC1–PC4 if available, otherwise as many as possible (min 2).
max_axes <- 4
k_ref <- ncol(ref_scores_all)
axes_to_use <- 1:min(max_axes, k_ref)
if (length(axes_to_use) < 2) stop("Reference has <2 axes; Procrustes needs at least 2.")

# Subset reference scores to the axes we will compare
ref_scores <- ref_scores_all[, axes_to_use, drop = FALSE]

cat("\n==============================\n")
cat("OOB Procrustes setup\n")
cat("==============================\n")
cat("Reference PCA:        ", ref_file, "\n", sep = "")
cat("Axes used:            PC", min(axes_to_use), "–PC", max(axes_to_use), "\n", sep = "")
cat("Reference species (n):", nrow(ref_scores), "\n\n")

## ----------------------------
## 4) Helper: standardize + project + procrustes for one bootstrap file
## ----------------------------
run_one_bootstrap <- function(traits_path, ref_scores, ref_means, ref_sds, ref_loadings, axes_to_use,
                              nrepet = 9999, seed = 1) {

  traits <- readRDS(traits_path)
  traits <- as.data.frame(traits)


  # Drop common non-trait columns if present (taxonomy / grouping)
  drop_cols <- intersect(colnames(traits), c("genus", "family", "order", "biomesKoeppenGroup", "biome"))
  if (length(drop_cols) > 0) {
    traits <- traits[, setdiff(colnames(traits), drop_cols), drop = FALSE]
  }


  # Basic checks
  if (is.null(rownames(traits)) || any(rownames(traits) == "")) {
    stop("Bootstrap traits must have rownames (species IDs): ", traits_path)
  }
  if (is.null(colnames(traits)) || any(colnames(traits) == "")) {
    stop("Bootstrap traits must have colnames (trait names): ", traits_path)
  }
  non_num <- colnames(traits)[!sapply(traits, is.numeric)]
  if (length(non_num) > 0) {
    stop("Non-numeric columns in bootstrap traits (", basename(traits_path), "): ",
         paste(non_num, collapse = ", "))
  }

  # Align trait columns to reference scaling vectors (by name if possible)
  # Reference means/sds may or may not be named; handle both cases.
  if (!is.null(names(ref_means)) && all(names(ref_means) %in% colnames(traits))) {
    traits <- traits[, names(ref_means), drop = FALSE]
  } else {
    # If unnamed, assume same order/length as the original reference trait matrix.
    if (length(ref_means) != ncol(traits) || length(ref_sds) != ncol(traits)) {
      stop("Trait column count mismatch vs reference means/sds in file: ", basename(traits_path))
    }
  }

  # Align loadings rows if they have names
  L <- ref_loadings
  if (!is.null(rownames(L)) && all(colnames(traits) %in% rownames(L))) {
    L <- L[colnames(traits), , drop = FALSE]
  }

  # Standardize using REFERENCE parameters
  Z <- scale(traits, center = ref_means, scale = ref_sds)

  # Project into reference space
  scores_proj_all <- as.data.frame(as.matrix(Z) %*% as.matrix(L))
  # Keep only axes_to_use
  scores_proj <- scores_proj_all[, axes_to_use, drop = FALSE]
  rownames(scores_proj) <- rownames(traits)

  # Match species with reference scores
  common_ids <- intersect(rownames(ref_scores), rownames(scores_proj))
  if (length(common_ids) < 10) {
    warning("Very small overlap (n=", length(common_ids), ") for file: ", basename(traits_path))
  }

  X_ref <- ref_scores[common_ids, , drop = FALSE]
  X_prj <- scores_proj[common_ids, , drop = FALSE]
  X_prj <- X_prj[rownames(X_ref), , drop = FALSE]
  stopifnot(identical(rownames(X_ref), rownames(X_prj)))

  # Procrustes permutation test
  set.seed(seed)
  proc <- ade4::procuste.rtest(as.data.frame(X_ref), as.data.frame(X_prj), nrepet = nrepet)

  data.frame(
    file = basename(traits_path),
    n_shared_species = nrow(X_ref),
    axes = paste0("PC", min(axes_to_use), "-PC", max(axes_to_use)),
    obs = unname(proc$obs),
    p_value = unname(proc$pvalue),
    stringsAsFactors = FALSE
  ) -> summary_row

  list(
    summary = summary_row,
    procrustes = proc,
    scores_ref = X_ref,
    scores_proj = X_prj
  )
}

## ----------------------------
## 5) Loop over both bootstrap sets
## ----------------------------
all_summaries <- list()
all_details <- list()

for (set_name in names(boot_dirs)) {

  boot_dir <- boot_dirs[[set_name]]
  if (!dir.exists(boot_dir)) stop("Missing bootstrap directory: ", boot_dir)

  files <- list.files(boot_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No .rds files found in: ", boot_dir)

  # Sort to ensure _001 ... _050 order
  files <- sort(files)

  cat("Running set:", set_name, "\n")
  cat("Bootstrap dir:", boot_dir, "\n")
  cat("Files found:", length(files), "\n\n")

  set_summaries <- list()
  set_details <- list()

  for (fp in files) {
    cat("  ->", basename(fp), "\n")
    res <- run_one_bootstrap(
      traits_path  = fp,
      ref_scores   = ref_scores,
      ref_means    = ref_means,
      ref_sds      = ref_sds,
      ref_loadings = ref_loadings,
      axes_to_use  = axes_to_use,
      nrepet       = 9999,
      seed         = 1
    )

    res$summary$set <- set_name
    set_summaries[[length(set_summaries) + 1]] <- res$summary

    # Keep full detail (can be large, but only 100 iters total)
    set_details[[basename(fp)]] <- res
  }

  all_summaries[[set_name]] <- do.call(rbind, set_summaries)
  all_details[[set_name]] <- set_details

  cat("\nFinished set:", set_name, "\n\n")
}

summary_df <- do.call(rbind, all_summaries)

## ----------------------------
## 6) Save outputs
## ----------------------------
# CSV summary
csv_file <- file.path(out_dir, "oob_procrustes_bootstrap_summary.csv")
write.csv(summary_df, csv_file, row.names = FALSE)

# RDS with all details
rds_file <- file.path(out_dir, "oob_procrustes_bootstrap_full_results.rds")
saveRDS(
  list(
    reference_file = ref_file,
    axes_to_use = axes_to_use,
    reference_n_species = nrow(ref_scores),
    summary = summary_df,
    details = all_details
  ),
  rds_file
)

cat("\n==============================\n")
cat("ALL DONE\n")
cat("==============================\n")
cat("Saved summary CSV:\n  ", csv_file, "\n", sep = "")
cat("Saved full results RDS:\n  ", rds_file, "\n", sep = "")
cat("==============================\n\n")

############################################################
# Optional quick summary prints
############################################################
cat("\nQuick summary (mean obs by set):\n")
print(aggregate(obs ~ set, data = summary_df, FUN = mean))

cat("\nQuick summary (min obs by set):\n")
print(aggregate(obs ~ set, data = summary_df, FUN = min))
