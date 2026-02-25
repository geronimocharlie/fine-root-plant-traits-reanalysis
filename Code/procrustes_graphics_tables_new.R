############################################################
# ----------------------------------------------------------
# Creates ONE final table + plots from:
#  (A) imputation test: global imputed reference vs not-imputed projected
#      -> reads Procrustes .rds outputs in Results/Procrustes_results/
#  (B) biome test (NEW groot setup): biome_groot_procrustes_summary.csv
#  (C) OOB test: oob_procrustes_bootstrap_summary.csv (m2 residuals)
#
# Outputs:
#  Results/Procrustes_results/final_procrustes_summary_table.csv
#  Results/Figures/
#     figure_final_biome_groot.png/pdf
#     figure_final_oob_R_equiv.png/pdf
#     figure_final_imputation.png/pdf
############################################################

## ----------------------------
## Packages
## ----------------------------
pkgs <- c("dplyr", "readr", "ggplot2", "stringr", "purrr", "tibble")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(purrr)
library(tibble)

## ----------------------------
## Paths
## ----------------------------
res_dir <- "Results"
proc_dir <- file.path(res_dir, "Procrustes_results")
fig_dir  <- file.path(res_dir, "Figures")

if (!dir.exists(proc_dir)) stop("Missing directory: ", proc_dir)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

## ==========================================================
## (A) IMPUTATION TEST(S): read saved Procrustes .rds outputs
## ==========================================================
# We assume files like:
#   procrustes_imputed_ref_vs_notimputed_projected_PC1PC2.rds
#   procrustes_imputed_ref_vs_notimputed_projected_PC1PC4.rds (if you have it)
#
# Each .rds was saved as a list with fields like:
#   $axes_compared, $n_shared_species, $procrustes (ade4 procuste.rtest object)
#
# If your object schema differs slightly, we robustly extract what we can.

imp_files <- list.files(proc_dir,
                        pattern = "^procrustes_imputed_ref_vs_notimputed.*\\.rds$",
                        full.names = TRUE)

read_imputation_rds <- function(fp) {
  x <- readRDS(fp)

  # try to extract procuste.rtest object
  proc_obj <- NULL
  if (is.list(x) && !is.null(x$procrustes)) proc_obj <- x$procrustes
  if (is.null(proc_obj) && inherits(x, "procustetest")) proc_obj <- x  # unlikely, but safe

  obs <- if (!is.null(proc_obj) && !is.null(proc_obj$obs)) unname(proc_obj$obs) else NA_real_
  pval <- if (!is.null(proc_obj) && !is.null(proc_obj$pvalue)) unname(proc_obj$pvalue) else NA_real_

  # axes label
  axes <- NA_character_
  if (is.list(x) && !is.null(x$axes_compared)) axes <- x$axes_compared
  if (is.na(axes)) {
    # infer from filename
    axes <- str_extract(basename(fp), "PC\\d+PC\\d+")
    if (!is.na(axes)) {
      axes <- str_replace_all(axes, "PC", "PC")
      axes <- str_replace(axes, "PC(\\d+)PC(\\d+)", "PC\\1–PC\\2")
    }
  }

  n_shared <- NA_integer_
  if (is.list(x) && !is.null(x$n_shared_species)) n_shared <- as.integer(x$n_shared_species)

  tibble(
    analysis = "imputation_vs_notimputed",
    comparison = "global imputed reference vs not-imputed (projected)",
    biome = NA_character_,
    set = NA_character_,
    axes_compared = axes,
    k_compare = suppressWarnings(as.integer(str_match(axes, "PC1–PC(\\d+)")[,2])),
    metric = "procrustes_similarity",
    obs = obs,
    p_value = pval,
    n_shared_species = n_shared,
    source_file = basename(fp)
  )
}

imputation_tbl <- if (length(imp_files) > 0) {
  map_dfr(imp_files, read_imputation_rds)
} else {
  tibble()
}

## ==========================================================
## (B) BIOME TEST (NEW): read biome_groot_procrustes_summary.csv
## ==========================================================
biome_csv <- file.path(proc_dir, "biome_groot_procrustes_summary.csv")
biome_tbl <- if (file.exists(biome_csv)) {
  read_csv(biome_csv, show_col_types = FALSE) %>%
    mutate(
      analysis = "biome_groot_vs_global",
      comparison = "biome PCA (groot) vs global imputed reference",
      set = NA_character_,
      metric = "procrustes_similarity",
      # standardize column name for similarity statistic
      obs = obs,
      p_value = p_value,
      axes_compared = axes_compared,
      biome = biome
    ) %>%
    select(
      analysis, comparison, biome, set,
      axes_compared, k_compare,
      metric, obs, p_value, n_shared_species,
      k_ref, k_biome, notes
    )
} else {
  tibble()
}

# ==========================================================
# (C) OOB: read oob_procrustes_bootstrap_summary.csv (m²)
#         and summarize to one row per regime (mean)
# ==========================================================
oob_csv <- file.path(proc_dir, "oob_procrustes_bootstrap_summary.csv")
oob_tbl <- tibble()

if (file.exists(oob_csv)) {
  o <- readr::read_csv(oob_csv, show_col_types = FALSE)

  # Standardize metric column name
  if ("statistic_value" %in% names(o)) {
    o <- dplyr::rename(o, m2 = statistic_value)
  } else if ("obs" %in% names(o)) {
    o <- dplyr::rename(o, m2 = obs)
  } else {
    stop("OOB CSV has neither 'statistic_value' nor 'obs'. Columns: ",
         paste(names(o), collapse = ", "))
  }

  # Standardize axes column name (optional)
  if ("axes" %in% names(o) && !("axes_compared" %in% names(o))) {
    o <- dplyr::rename(o, axes_compared = axes)
  }

  # Ensure required columns exist
  required <- c("set", "m2", "n_shared_species")
  missing <- setdiff(required, names(o))
  if (length(missing) > 0) {
    stop("OOB CSV missing required columns: ", paste(missing, collapse = ", "),
         "\nAvailable columns: ", paste(names(o), collapse = ", "))
  }

  # Add R-equivalent for interpretability (optional)
  o <- o %>%
    dplyr::mutate(r_equiv = sqrt(pmax(0, 1 - m2)))

  # Summarize to one row per uncertainty regime
  oob_tbl <- o %>%
    dplyr::group_by(set) %>%
    dplyr::summarise(
      analysis = "oob_imputation_uncertainty",
      comparison = "bootstrap imputations vs mean-imputation reference",
      biome = NA_character_,
      axes_compared = dplyr::first(if ("axes_compared" %in% names(o)) axes_compared else NA_character_),
      k_compare = NA_integer_,
      metric = "procrustes_residual_m2",
      obs = mean(m2, na.rm = TRUE),
      # Optional extra columns (highly recommended for context)
      obs_sd  = sd(m2, na.rm = TRUE),
      obs_min = min(m2, na.rm = TRUE),
      obs_max = max(m2, na.rm = TRUE),
      r_equiv_mean = mean(r_equiv, na.rm = TRUE),
      n_shared_species = dplyr::first(n_shared_species),
      p_value = NA_real_,
      .groups = "drop"
    )
} else {
  message("WARNING: OOB summary CSV not found at: ", oob_csv)
}
## ==========================================================
## Combine final summary table
## ==========================================================
final_tbl <- bind_rows(
  imputation_tbl,
  biome_tbl %>% mutate(set = NA_character_) %>%
    select(any_of(names(imputation_tbl)), everything()),
  oob_tbl
) %>%
  # keep a stable column order
  relocate(analysis, comparison, biome, set, axes_compared, k_compare, metric, obs, p_value, n_shared_species) %>%
  arrange(analysis, biome, set, axes_compared)

out_csv <- file.path("Results/Tables/final_procrustes_summary_table.csv")
write_csv(final_tbl, out_csv)
cat("Saved final summary table:\n  ", out_csv, "\n", sep = "")

## ==========================================================
## PLOTS
## ==========================================================

## ---- Plot 1: Biome groot (if available)
if (nrow(biome_tbl) > 0) {
  biome_plot <- biome_tbl %>%
    filter(!is.na(obs), !is.na(n_shared_species)) %>%
    mutate(
      biome = factor(biome, levels = c("temperate", "tropical", "continental")),
      congruence_level = case_when(
        obs < 0.8 ~ "below 0.8",
        obs < 0.9 ~ "0.8–0.9",
        TRUE ~ "≥ 0.9"
      ),
      label = paste0("n=", n_shared_species, "\n",
                     "k_used=", k_compare, " (k_biome=", k_biome, ")\n",
                     axes_compared)
    )

  y_min <- max(0.65, min(biome_plot$obs, na.rm = TRUE) - 0.05)
  y_max <- min(1.01, max(biome_plot$obs, na.rm = TRUE) + 0.02)

  p_biome <- ggplot(biome_plot, aes(x = biome, y = obs)) +
    geom_hline(yintercept = 0.9, linetype = 2, linewidth = 0.8) +
    geom_hline(yintercept = 0.8, linetype = 3, linewidth = 0.8) +
    geom_segment(aes(xend = biome, y = y_min, yend = obs), linewidth = 0.8, alpha = 0.35) +
    geom_point(aes(color = congruence_level), size = 3.2) +
    geom_text(aes(label = sprintf("%.2f", obs)), vjust = -0.8, size = 3.4) +
    geom_text(aes(label = label), vjust = 1.7, size = 3) +
    scale_color_manual(
      values = c("≥ 0.9" = "black", "0.8–0.9" = "darkred", "below 0.8" = "firebrick4"),
      name = "Congruence level"
    ) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(
      title = "Biome test (groot): biome PCA vs global imputed reference space",
      x = "Biome",
      y = "Procrustes similarity (higher = more similar)",
      caption = "Lines: 0.9 = high congruence; 0.8 = moderate congruence. Labels show n and dimensionality used."
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(fig_dir, "figure_final_biome_groot.png"), p_biome, width = 9, height = 6, dpi = 300)
  ggsave(file.path(fig_dir, "figure_final_biome_groot.pdf"), p_biome, width = 9, height = 6)
  cat("Saved biome figure to Results/Figures.\n")
}

## ---- Plot 2: OOB (R-equivalent) (if available)
if (nrow(oob_tbl) > 0) {
  oob_for_plot <- oob_tbl %>%
    dplyr::mutate(set = factor(set, levels = c("oob_1x", "oob_2x")))

  p_oob <- ggplot(oob_for_plot, aes(x = set, y = r_equiv)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, alpha = 0.6) +
    labs(
      title = "OOB test: stability under imputation uncertainty",
      x = "Uncertainty regime",
      y = "Procrustes similarity (R-equivalent; higher = more similar)",
      caption = "R-equivalent computed as sqrt(1 - m²), where m² is the Procrustes residual."
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(fig_dir, "figure_final_oob_R_equiv.png"), p_oob, width = 7, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, "figure_final_oob_R_equiv.pdf"), p_oob, width = 7, height = 5)
}
## ---- Plot 3: Imputation test(s) as points (if available)
if (nrow(imputation_tbl) > 0) {
  imp_plot <- imputation_tbl %>%
    mutate(label = paste0(axes_compared, "\n", "n=", n_shared_species))

  p_imp <- ggplot(imp_plot, aes(x = axes_compared, y = obs)) +
    geom_point(size = 3.2) +
    geom_text(aes(label = label), vjust = 1.5, size = 3) +
    geom_hline(yintercept = 0.9, linetype = 2, linewidth = 0.8) +
    geom_hline(yintercept = 0.8, linetype = 3, linewidth = 0.8) +
    coord_cartesian(ylim = c(max(0.65, min(imp_plot$obs) - 0.05), 1.01)) +
    labs(
      title = "Imputation test: not-imputed projected into global imputed reference",
      x = "Axes compared",
      y = "Procrustes similarity (higher = more similar)",
      caption = "Lines: 0.9 = high congruence; 0.8 = moderate congruence."
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(fig_dir, "figure_final_imputation.png"), p_imp, width = 7, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, "figure_final_imputation.pdf"), p_imp, width = 7, height = 5)
  cat("Saved imputation figure to Results/Figures.\n")
}

cat("\n==============================\nDONE\n==============================\n")