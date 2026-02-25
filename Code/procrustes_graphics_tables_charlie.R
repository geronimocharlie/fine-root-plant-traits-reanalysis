############################################################

# ----------------------------------------------------------
# Loads results for:
#   (1) Imputation test (imputed vs not-imputed)
#   (2) Biome tests (3 biomes × 2 domains)
#   (3) OOB uncertainty tests (bootstrap 1x vs 2x)
#
# Produces:
#   - Results/Figures/figure_biome_procrustes.(png/pdf)
#   - Results/Figures/figure_oob_procrustes.(png/pdf)
#   - Results/Tables/procrustes_tests_summary.csv
#
# Ensures table includes:
#   - test category
#   - comparison
#   - biome / domain / set (as applicable)
#   - n_shared_species
#   - axes compared (PCs)
#   - statistic value + p-value
#
# IMPORTANT:
#   ade4::procuste.rtest can output different statistics across setups
#   (sometimes "correlation-like" close to 1, sometimes residual-like close to 0).
#   In our pipeline:
#     - Imputation + Biome scripts produced values near 1 (interpret: higher = more similar)
#     - OOB script produced values ~0.05 (interpret: lower = more similar; residual-like)
#   We keep them as separate "statistic_type" and do not mix them on the same y-axis.
############################################################

## ----------------------------
## 0) Packages
## ----------------------------
pkgs <- c("ggplot2", "dplyr", "readr", "stringr", "tibble")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

## ----------------------------
## 1) Paths
## ----------------------------
dir_fig <- file.path("Results", "Figures")
dir_tab <- file.path("Results", "Tables")
if (!dir.exists(dir_fig)) dir.create(dir_fig, recursive = TRUE)
if (!dir.exists(dir_tab)) dir.create(dir_tab, recursive = TRUE)

# Imputation test output (adjust name if you used PC1PC2 instead of PC1PC4)
imp_rds_candidates <- c(
  file.path("Results", "Procrustes_results", "procrustes_imputed_ref_vs_notimputed_projected_PC1PC4.rds"),
  file.path("Results", "Procrustes_results", "procrustes_imputed_ref_vs_notimputed_projected_PC1PC2.rds")
)

# Biome summary CSV (prefer FIXEDAXES)
biome_csv_candidates <- c(
  file.path("Results", "Procrustes_results", "biome_procrustes_summary_FIXEDAXES.csv"),
  file.path("Results", "Procrustes_results", "biome_procrustes_summary.csv")
)

# OOB outputs
oob_csv <- file.path("Results", "Procrustes_results", "oob_procrustes_bootstrap_summary.csv")

## ----------------------------
## 2) Helpers
## ----------------------------
pick_first_existing <- function(paths) {
  for (p in paths) if (file.exists(p)) return(p)
  return(NA_character_)
}

safe_read_rds <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  readRDS(path)
}

## ----------------------------
## 3) Load Imputation test (RDS)
## ----------------------------
imp_rds <- pick_first_existing(imp_rds_candidates)
imp_obj <- safe_read_rds(imp_rds)

imp_df <- tibble()
if (!is.null(imp_obj)) {
  # Our saved list contains:
  #  - n_shared_species
  #  - axes_used (character vector) OR "axes_used" field
  #  - procrustes object with $obs and $pvalue
  # We make minimal assumptions and fail loudly if missing.
  stopifnot(!is.null(imp_obj$n_shared_species))
  stopifnot(!is.null(imp_obj$procrustes))
  stopifnot(!is.null(imp_obj$procrustes$obs))
  stopifnot(!is.null(imp_obj$procrustes$pvalue))

  axes_used <- imp_obj$axes_used
  if (is.null(axes_used)) axes_used <- NA_character_

  # normalize axes display
  axes_label <- if (length(axes_used) > 1) {
    paste0(min(axes_used), "–", max(axes_used))  # if it's like c("PC1","PC2"...), this isn't numeric
  } else {
    as.character(axes_used)
  }
  # Better: if it's "PC1","PC2" etc:
  if (all(str_detect(axes_used, "^PC\\d+$"))) {
    nums <- as.integer(str_remove(axes_used, "PC"))
    axes_label <- paste0("PC", min(nums), "–PC", max(nums))
  } else if (is.character(axes_used) && length(axes_used) == 1 && !is.na(axes_used)) {
    axes_label <- axes_used
  } else if (is.na(axes_used)[1]) {
    axes_label <- "NA"
  }

  imp_df <- tibble(
    test_category  = "imputation",
    comparison     = "imputed reference vs not-imputed projected",
    biome          = NA_character_,
    domain         = "combined",
    set            = NA_character_,
    n_shared_species = as.integer(imp_obj$n_shared_species),
    axes_compared  = axes_label,
    statistic_type = "procrustes_similarity (higher = more similar)",
    statistic_value = as.numeric(imp_obj$procrustes$obs),
    p_value        = as.numeric(imp_obj$procrustes$pvalue),
    source_file    = imp_rds
  )
} else {
  message("WARNING: Imputation RDS not found. Looked for:\n  - ",
          paste(imp_rds_candidates, collapse = "\n  - "))
}

## ----------------------------
## 4) Load Biome tests (CSV summary)
## ----------------------------
biome_csv <- pick_first_existing(biome_csv_candidates)
biome_df <- tibble()

if (!is.na(biome_csv) && file.exists(biome_csv)) {
  b <- read_csv(biome_csv, show_col_types = FALSE)

  # Expected columns from our scripts:
  # biome, domain, n_shared_species, (reference_axes or axes_used), obs, p_value
  # We'll support both variants:
  if (!("biome" %in% names(b))) stop("Biome summary is missing column 'biome': ", biome_csv)
  if (!("domain" %in% names(b))) stop("Biome summary is missing column 'domain': ", biome_csv)
  if (!("n_shared_species" %in% names(b))) stop("Biome summary missing 'n_shared_species': ", biome_csv)
  if (!("obs" %in% names(b))) stop("Biome summary missing 'obs': ", biome_csv)
  if (!("p_value" %in% names(b))) stop("Biome summary missing 'p_value': ", biome_csv)

  # axes column name may differ
  axes_col <- if ("reference_axes" %in% names(b)) "reference_axes" else if ("axes_used" %in% names(b)) "axes_used" else NA_character_
  if (is.na(axes_col)) {
    # fallback
    b$reference_axes <- NA_character_
    axes_col <- "reference_axes"
  }

  biome_df <- b %>%
    transmute(
      test_category    = "biome",
      comparison       = "biome-only PCA vs global imputed reference space",
      biome            = as.character(biome),
      domain           = as.character(domain),
      set              = NA_character_,
      n_shared_species = as.integer(n_shared_species),
      axes_compared    = as.character(.data[[axes_col]]),
      statistic_type   = "procrustes_similarity (higher = more similar)",
      statistic_value  = as.numeric(obs),
      p_value          = as.numeric(p_value),
      source_file      = biome_csv
    )
} else {
  message("WARNING: Biome summary CSV not found. Looked for:\n  - ",
          paste(biome_csv_candidates, collapse = "\n  - "))
}

## ----------------------------
## 5) Load OOB tests (CSV summary)
## ----------------------------
oob_df <- tibble()

if (file.exists(oob_csv)) {
  o <- read_csv(oob_csv, show_col_types = FALSE)

  # expected columns from our OOB script:
  # file, n_shared_species, axes, obs, p_value, set
  required <- c("set", "n_shared_species", "axes", "obs")
  missing <- setdiff(required, names(o))
  if (length(missing) > 0) stop("OOB summary missing columns: ", paste(missing, collapse = ", "))

  # p_value may exist; if not, keep NA
  if (!("p_value" %in% names(o))) o$p_value <- NA_real_

  oob_df <- o %>%
    transmute(
      test_category    = "oob",
      comparison       = "bootstrap-imputed projection vs mean-imputation reference",
      biome            = NA_character_,
      domain           = "combined",
      set              = as.character(set),
      n_shared_species = as.integer(n_shared_species),
      axes_compared    = as.character(axes),
      statistic_type   = "procrustes_residual_m2 (lower = more similar)",
      statistic_value  = as.numeric(obs),
      p_value          = as.numeric(p_value),
      source_file      = oob_csv
    )
} else {
  message("WARNING: OOB summary CSV not found at: ", oob_csv)
}

## ----------------------------
## 6) Master table (all tests)
## ----------------------------
summary_all <- bind_rows(imp_df, biome_df, oob_df) %>%
  mutate(
    p_value = ifelse(is.na(p_value), NA_real_, p_value)
  )

out_table <- file.path(dir_tab, "procrustes_tests_summary.csv")
write.csv(summary_all, out_table, row.names = FALSE)

cat("\nSaved master summary table:\n  ", out_table, "\n", sep = "")

## ----------------------------
## 7) Figure A: Biome results (similarity-style, higher = better)
## ----------------------------

# ============================================================
# Final Biome plot:
# - Facet (above / below)
# - Lollipop
# - n once in x-axis labels
# - Reference lines at 0.9 (high) and 0.8 (moderate)
# - Color highlighting below thresholds
# - Tighter zoom
# - Title unchanged
# ============================================================


# ---- Prepare data
biome_plot <- biome_df %>%
  mutate(
    biome  = factor(biome, levels = c("temperate", "tropical", "continental")),
    domain = factor(domain, levels = c("above", "below")),
    congruence_level = case_when(
      statistic_value < 0.8 ~ "below 0.8",
      statistic_value < 0.9 ~ "0.8–0.9",
      TRUE ~ "≥ 0.9"
    )
  )

# ---- n per biome (only once in axis label)
n_df <- biome_plot %>%
  group_by(biome) %>%
  summarise(n = first(n_shared_species), .groups = "drop")

biome_labels <- setNames(
  paste0(levels(biome_plot$biome), "\n", "n=", 
         n_df$n[match(levels(biome_plot$biome), n_df$biome)]),
  levels(biome_plot$biome)
)

# ---- Domain facet labels with axis mapping
domain_labs <- c(
  above = "Aboveground (PC1–PC2)",
  below = "Belowground (PC3–PC4)"
)

# ---- Zoom range
y_min <- 0.72
y_max <- 1.01

# ---- Plot
p_biome_final <- ggplot(biome_plot, aes(x = biome, y = statistic_value)) +

  # Reference lines
  geom_hline(yintercept = 0.9, linetype = 2, linewidth = 0.8) +
  geom_hline(yintercept = 0.8, linetype = 3, linewidth = 0.8) +

  # Lollipop stems
  geom_segment(aes(xend = biome, y = y_min, yend = statistic_value),
               linewidth = 0.8, alpha = 0.35) +

  # Points with threshold-based coloring
  geom_point(aes(color = congruence_level), size = 3.2) +

  # Value labels
  geom_text(aes(label = sprintf("%.2f", statistic_value)),
            vjust = -0.8, size = 3.4) +

  # Facet by domain
  facet_wrap(~domain, ncol = 1, labeller = labeller(domain = domain_labs)) +

  scale_x_discrete(labels = biome_labels) +

  scale_color_manual(
    values = c(
      "≥ 0.9"   = "black",
      "0.8–0.9" = "darkred",
      "below 0.8" = "firebrick4"
    ),
    name = "Congruence level"
  ) +

  coord_cartesian(ylim = c(y_min, y_max)) +

  labs(
    title = "Biome test: biome-only PCA vs global imputed reference space",
    x = "Biome",
    y = "Procrustes similarity (higher = more similar)",
    #caption = paste(
    #  "Reference lines:",
    #  "0.9 = high congruence;",
    #  "0.8 = moderate congruence.",
    #  "Aboveground uses PC1–PC2; Belowground uses PC3–PC4."
    #)
  ) +

  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

# ---- Save
ggsave("Results/Figures/figure_biome_procrustes_final.png",
       p_biome_final, width = 9, height = 6, dpi = 300)
ggsave("Results/Figures/figure_biome_procrustes_final.pdf",
       p_biome_final, width = 9, height = 6)

p_biome_final


## ----------------------------
## 8) Figure B: OOB results (residual-style, lower = better)
## ----------------------------
if (nrow(oob_df) > 0) {

  oob_plot <- oob_df %>%
    mutate(
      set = factor(set, levels = c("oob_1x", "oob_2x"))
    )

  p2 <- ggplot(oob_plot, aes(x = set, y = statistic_value)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, alpha = 0.6) +
    labs(
      title = "OOB test: projection stability under imputation uncertainty",
      x = "Imputation uncertainty set",
      y = "Procrustes residual m² (lower = more similar)"
    ) +
    theme_minimal()

  ggsave(file.path(dir_fig, "figure_oob_procrustes.png"), p2, width = 7, height = 5, dpi = 300)
  ggsave(file.path(dir_fig, "figure_oob_procrustes.pdf"), p2, width = 7, height = 5)

  cat("Saved OOB figure:\n  ",
      file.path(dir_fig, "figure_oob_procrustes.png"),
      "\n", sep = "")
}

## ----------------------------
## 9) Optional: A tiny one-line console report
## ----------------------------
cat("\n--- Quick report ---\n")
if (nrow(imp_df) > 0) {
  cat("Imputation test:", round(imp_df$statistic_value[1], 4),
      "| axes:", imp_df$axes_compared[1],
      "| n:", imp_df$n_shared_species[1], "\n")
}
if (nrow(biome_df) > 0) {
  cat("Biome tests (min/max similarity):",
      round(min(biome_df$statistic_value, na.rm = TRUE), 4), "/",
      round(max(biome_df$statistic_value, na.rm = TRUE), 4), "\n")
}
if (nrow(oob_df) > 0) {
  cat("OOB tests (mean residual by set):\n")
  print(oob_df %>% group_by(set) %>% summarise(mean_m2 = mean(statistic_value, na.rm = TRUE)))
}
cat("--------------------\n")
 
 # ============================================================
# 8(a) OOB plot using R-equivalent derived from m^2
# ============================================================

# oob_df expected columns:
# set, statistic_value (= m^2), n_shared_species, axes_compared
stopifnot(nrow(oob_df) > 0)

oob_plot <- oob_df %>%
  mutate(
    m2 = statistic_value,
    # Convert to an intuitive similarity scale (higher = more similar)
    # r_equiv = sqrt(1 - m2)
    r_equiv = sqrt(pmax(0, 1 - m2)),
    set = factor(set, levels = c("oob_1x", "oob_2x"))
  )

p_oob_r <- ggplot(oob_plot, aes(x = set, y = r_equiv)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(
    title = "OOB test: projection stability under imputation uncertainty",
    x = "Imputation uncertainty set",
    y = "Procrustes similarity (R-equivalent; higher = more similar)",
    caption = "R-equivalent computed as sqrt(1 - m²), where m² is the Procrustes residual from ade4::procuste.rtest."
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path("Results", "Figures", "figure_oob_procrustes_R_equiv.png"),
       p_oob_r, width = 7, height = 5, dpi = 300)
ggsave(file.path("Results", "Figures", "figure_oob_procrustes_R_equiv.pdf"),
       p_oob_r, width = 7, height = 5)

p_oob_r

