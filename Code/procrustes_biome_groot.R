############################################################
# ----------------------------------------------------------
# Biome test (NEW setup):
# - One PCA per biome on combined/groot traits
# - Compare each biome PCA configuration to the GLOBAL imputed PCA space
# - Handle underrepresented biomes where PCA retains <4 dimensions:
#     k_compare = min(4, k_ref, k_biome), must be >= 2
#
# Inputs (edit if needed):
#   Results/PCA_imputed_combined_full.rds              (GLOBAL reference)
#   Results/PCA_imputed_temperate_groot.rds
#   Results/PCA_imputed_tropical_groot.rds
#   Results/PCA_imputed_continental_groot.rds
#
# Outputs:
#   Results/Procrustes_results/biome_groot_procrustes_<biome>.rds
#   Results/Procrustes_results/biome_groot_procrustes_summary.csv
############################################################

## ----------------------------
## 0) Packages
## ----------------------------
pkgs <- c("ade4", "dplyr", "readr")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(ade4)
library(dplyr)
library(readr)

set.seed(1)

## ----------------------------
## 1) Paths
## ----------------------------
ref_file <- file.path("Results", "PCA_imputed_combined_full.rds")  # <-- change if needed

biome_pca_files <- list(
  temperate   = file.path("Results", "PCA_imputed_temperate_groot.rds"),
  tropical    = file.path("Results", "PCA_imputed_tropical_groot.rds"),
  continental = file.path("Results", "PCA_imputed_continental_groot.rds")
)

out_dir <- file.path("Results", "Procrustes_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ----------------------------
## 2) Helper checks
## ----------------------------
get_scores <- function(pca_obj, file_label = "PCA object") {
  if (is.null(pca_obj$traitsUse)) stop(file_label, " is missing $traitsUse.")
  X <- as.data.frame(pca_obj$traitsUse)
  if (is.null(rownames(X)) || any(rownames(X) == "")) stop(file_label, " traitsUse has no valid rownames (species IDs).")
  if (ncol(X) < 1) stop(file_label, " traitsUse has 0 columns.")
  X
}

## ----------------------------
## 3) Load reference PCA
## ----------------------------
if (!file.exists(ref_file)) stop("Missing global reference PCA: ", ref_file)
ref <- readRDS(ref_file)
X_ref_all <- get_scores(ref, "Global reference PCA")

k_ref <- ncol(X_ref_all)
cat("\n==============================\n")
cat("Biome Procrustes (groot) setup\n")
cat("==============================\n")
cat("Reference PCA file: ", ref_file, "\n", sep = "")
cat("Reference dimensions available: ", k_ref, "\n", sep = "")
cat("==============================\n\n")

## ----------------------------
## 4) Run per-biome Procrustes
## ----------------------------
summary_rows <- list()

for (biome in names(biome_pca_files)) {

  biome_file <- biome_pca_files[[biome]]
  if (!file.exists(biome_file)) stop("Missing biome PCA file: ", biome_file)

  cat("--------------------------------------------------\n")
  cat("Biome:", biome, "\n")
  cat("Biome PCA file:", biome_file, "\n")

  b <- readRDS(biome_file)
  X_biome_all <- get_scores(b, paste0("Biome PCA (", biome, ")"))

  k_biome <- ncol(X_biome_all)

  # choose shared dimensionality for Procrustes
  k_compare <- min(4, k_ref, k_biome)

  if (k_compare < 2) {
    warning("Biome ", biome, " has k_compare < 2 (k_biome=", k_biome, "). Skipping Procrustes.")
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      biome = biome,
      n_shared_species = NA_integer_,
      k_ref = k_ref,
      k_biome = k_biome,
      k_compare = k_compare,
      axes_compared = NA_character_,
      obs = NA_real_,
      p_value = NA_real_,
      notes = "Skipped: insufficient dimensionality (<2)",
      stringsAsFactors = FALSE
    )
    next
  }

  axes <- 1:k_compare
  X_ref <- X_ref_all[, axes, drop = FALSE]
  X_bio <- X_biome_all[, axes, drop = FALSE]

  # match shared species
  common_ids <- intersect(rownames(X_ref), rownames(X_bio))
  X_ref2 <- X_ref[common_ids, , drop = FALSE]
  X_bio2 <- X_bio[common_ids, , drop = FALSE]
  X_bio2 <- X_bio2[rownames(X_ref2), , drop = FALSE]
  stopifnot(identical(rownames(X_ref2), rownames(X_bio2)))

  cat("Shared species used:", nrow(X_ref2), "\n")
  cat("k_ref:", k_ref, "| k_biome:", k_biome, "| k_compare:", k_compare, "\n", sep = " ")
  cat("Axes compared: PC1–PC", k_compare, "\n", sep = "")

  # Procrustes permutation test
  set.seed(1)
  proc <- ade4::procuste.rtest(as.data.frame(X_ref2), as.data.frame(X_bio2), nrepet = 9999)
  print(proc)

  # save per-biome output
  out <- list(
    biome = biome,
    reference_file = ref_file,
    biome_file = biome_file,
    n_shared_species = nrow(X_ref2),
    k_ref = k_ref,
    k_biome = k_biome,
    k_compare = k_compare,
    axes_compared = paste0("PC1–PC", k_compare),
    reference_scores = X_ref2,
    biome_scores = X_bio2,
    procrustes = proc
  )

  out_file <- file.path(out_dir, paste0("biome_groot_procrustes_", biome, ".rds"))
  saveRDS(out, out_file)
  cat("Saved:", out_file, "\n\n")

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    biome = biome,
    n_shared_species = nrow(X_ref2),
    k_ref = k_ref,
    k_biome = k_biome,
    k_compare = k_compare,
    axes_compared = paste0("PC1–PC", k_compare),
    obs = unname(proc$obs),
    p_value = unname(proc$pvalue),
    notes = "",
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(summary_rows)

summary_csv <- file.path(out_dir, "biome_groot_procrustes_summary.csv")
write.csv(summary_df, summary_csv, row.names = FALSE)

cat("\n==============================\n")
cat("ALL DONE\n")
cat("==============================\n")
cat("Summary written to:\n  ", summary_csv, "\n", sep = "")
cat("==============================\n\n")



############################################################
# Code/06_plot_biome_groot_procrustes.R
# ----------------------------------------------------------
# Reads:
#   Results/Procrustes_results/biome_groot_procrustes_summary.csv
# Creates:
#   - congruence plot with 0.9/0.8 lines and labels (n, k)
#   - retained dimension barplot (k_biome)
############################################################

library(dplyr)
library(readr)
library(ggplot2)

in_csv <- file.path("Results", "Procrustes_results", "biome_groot_procrustes_summary.csv")
if (!file.exists(in_csv)) stop("Missing summary csv: ", in_csv)

df <- read_csv(in_csv, show_col_types = FALSE)

# Drop skipped (k_compare < 2) from plotting
df_plot <- df %>% filter(!is.na(obs), !is.na(n_shared_species))

# Decide what "obs" means:
# In your previous biome tests obs behaved like similarity (~0.7–1.0).
# We'll label it as "Procrustes similarity" consistently.
# (If you ever see obs ~0.05 here, then it's residual-like; we'd adjust label.)

df_plot <- df_plot %>%
  mutate(
    biome = factor(biome, levels = c("temperate", "tropical", "continental")),
    congruence_level = case_when(
      obs < 0.8 ~ "below 0.8",
      obs < 0.9 ~ "0.8–0.9",
      TRUE ~ "≥ 0.9"
    ),
    label = paste0("n=", n_shared_species, "\n",
                   "k_biome=", k_biome, ", k_used=", k_compare, "\n",
                   axes_compared)
  )

dir_fig <- file.path("Results", "Figures")
if (!dir.exists(dir_fig)) dir.create(dir_fig, recursive = TRUE)

# Plot 1: congruence
y_min <- max(0.65, min(df_plot$obs, na.rm = TRUE) - 0.05)
y_max <- min(1.01, max(df_plot$obs, na.rm = TRUE) + 0.02)

p1 <- ggplot(df_plot, aes(x = biome, y = obs)) +
  geom_hline(yintercept = 0.9, linetype = 2, linewidth = 0.8) +
  geom_hline(yintercept = 0.8, linetype = 3, linewidth = 0.8) +
  geom_segment(aes(xend = biome, y = y_min, yend = obs),
               linewidth = 0.8, alpha = 0.35) +
  geom_point(aes(color = congruence_level), size = 3.2) +
  geom_text(aes(label = sprintf("%.2f", obs)), vjust = -0.8, size = 3.4) +
  geom_text(aes(label = label), vjust = 1.7, size = 3) +
  scale_color_manual(
    values = c("≥ 0.9" = "black", "0.8–0.9" = "darkred", "below 0.8" = "firebrick4"),
    name = "Congruence level"
  ) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(
    title = "Biome test: biome PCA (groot) vs global imputed reference space",
    x = "Biome",
    y = "Procrustes similarity (higher = more similar)",
    caption = "Reference lines: 0.9 = high congruence; 0.8 = moderate congruence. Labels show n and retained dimensions."
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(dir_fig, "figure_biome_groot_procrustes.png"), p1, width = 9, height = 6, dpi = 300)
ggsave(file.path(dir_fig, "figure_biome_groot_procrustes.pdf"), p1, width = 9, height = 6)

# Plot 2: retained dimensionality (collapse view)
p2 <- ggplot(df_plot, aes(x = biome, y = k_biome)) +
  geom_col() +
  labs(
    title = "Retained PCA dimensionality per biome (paran)",
    x = "Biome",
    y = "Retained dimensions (k_biome)"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(dir_fig, "figure_biome_groot_retained_dimensions.png"), p2, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(dir_fig, "figure_biome_groot_retained_dimensions.pdf"), p2, width = 7, height = 4.5)

p1
p2
