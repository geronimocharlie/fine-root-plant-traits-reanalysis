library(ggplot2)
library(MASS)
library(dplyr)
library(tidyr)
library(stringr)

setwd("C:/Users/katha/OneDrive/2_documents/4_Studium/EMDS_Semester3/Capstone")

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
  
  scale_factor <- min(
    max(abs(scores$C1)) / max(abs(loadings$C1)),
    max(abs(scores$C2)) / max(abs(loadings$C2))
  )
  
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



# ==== Pipeline =====
fun_pipeline <- function(file_name, plot_name, rotation="varimax"){
  PCA_file <- readRDS(paste0("./PCA_results/", file_name, ".rds"))
  plot_12 <- plotPCA(PCA_dataset = PCA_file,
                                           axes = c(1,2),
                                           probs = c(0.50, 0.95, 0.99),
                                           title = paste0(plot_name, " (PC1 vs PC2)")) #c("ph", "sla", "ln", "sm", "la", "ssd"))
                     
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
                                             title = paste0(plot_name, " (PC3 vs PC4)")) #c("SRL", "D", "RTD", "N"))
                       
    ggsave(
      filename = paste0("./PCA_plots/", plot_name, "_PC3_PC4.png"),
      plot = plot_34,
      width = 18,
      height = 15,
      units = "cm",
      dpi = 300
    )
  }}


# # ==== Apply Pipeline ====
# fun_pipeline("PCA_imputed_combined", "Imputed Full")
# fun_pipeline("PCA_imputed_combined_half" , "Imputed Half")
# fun_pipeline("PCA_imputed_combined_quarter", "Imputed Quarter")
# fun_pipeline("PCA_imputed_combined_tenth", "Imputed Tenth")
# # 
# # # biomes GRooT
# fun_pipeline("PCA_imputed_continental_groot", "Continental GRooT")
# fun_pipeline("PCA_imputed_temperate_groot", "Temperate GRooT")
# fun_pipeline("PCA_imputed_tropical_groot", "Tropical GRooT")
# # 
# # # 9 biomes
# fun_pipeline("PCA_imputed_bf", "Boreal Forest")
# fun_pipeline("PCA_imputed_tgd", "Temperate Grassland Desert")
# fun_pipeline("PCA_imputed_trf", "Temperate Rain Forest")
# fun_pipeline("PCA_imputed_troprf", "Tropical Rain Forest")
# fun_pipeline("PCA_imputed_tropss", "Tropical Seasonal Forest Savanna")
# fun_pipeline("PCA_imputed_tsf", "Temperate Seasonal Forest")
# fun_pipeline("PCA_imputed_ws", "Woodland Shrubland")
# fun_pipeline("PCA_imputed_sd", "Subtropical Desert")
# # 
# # # 9 biomes categorized
# fun_pipeline("PCA_imputed_continental", "Continental")
# fun_pipeline("PCA_imputed_temperate", "Temperate")
# fun_pipeline("PCA_imputed_tropical", "Tropical")
# fun_pipeline("PCA_imputed_multiple", "Multiple")
# fun_pipeline("PCA_imputed_other", "Other")

# ==== PART 2: LOADING TABLES ====
# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
setwd("C:/Users/katha/OneDrive/2_documents/4_Studium/EMDS_Semester3/Capstone")

# Trait order (above lowercase, below uppercase)
trait_order <- c("ph", "ssd", "sm", "la", "ln", "sla", "SRL", "D", "RTD", "N")

# Helper to extract loadings matrix
extract_load <- function(mat, n) {
  m <- as.matrix(mat[trait_order, 1:n, drop = FALSE])
  colnames(m) <- paste0("PC", 1:n)
  m
}

# Build long data frame with optional SD column
make_df <- function(mat, method, n, sd_mat = NULL) {
  df <- as.data.frame(mat, optional = TRUE, stringsAsFactors = FALSE) %>%
    setNames(paste0("PC", 1:n)) %>%
    mutate(Trait = rownames(mat)) %>%
    pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "value") %>%
    mutate(Method = method)
  
  if (!is.null(sd_mat)) {
    sd_df <- as.data.frame(sd_mat, optional = TRUE, stringsAsFactors = FALSE) %>%
      setNames(paste0("PC", 1:n)) %>%
      mutate(Trait = rownames(sd_mat)) %>%
      pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "sd")
    df <- left_join(df, sd_df, by = c("Trait", "PC"))
  } else {
    df$sd <- NA_real_
  }
  df
}


# PCA stability
PCA_combined <- readRDS("PCA_results/PCA_imputed_combined.rds")
PCA_half     <- readRDS("PCA_results/PCA_imputed_combined_half.rds")
PCA_quarter  <- readRDS("PCA_results/PCA_imputed_combined_quarter.rds")
PCA_tenth  <- readRDS("PCA_results/PCA_imputed_combined_tenth.rds")

load_original <- extract_load(PCA_combined$PCA$loadings, 4)
load_half     <- extract_load(PCA_half$PCA$loadings, 4)
load_quarter       <- extract_load(PCA_quarter$PCA$loadings, 3)
load_tenth       <- extract_load(PCA_tenth$PCA$loadings, 3)

long_df1 <- bind_rows(
  make_df(load_original, "Full", 4),
  make_df(load_half,     "Half", 4),
  make_df(load_quarter,   "Quarter", 3),
  make_df(load_tenth,     "Tenth", 3)
)


# Flag the largest absolute value per Trait & Method
long_df1 <- long_df1 %>%
  group_by(Method, Trait) %>%
  mutate(is_max = abs(value) == max(abs(value))) %>%
  ungroup()

# biomes groot
PCA_temperate_groot <- readRDS("PCA_results/PCA_imputed_temperate_groot.rds")
PCA_continental_groot <- readRDS("PCA_results/PCA_imputed_continental_groot.rds")
PCA_tropical_groot <- readRDS("PCA_results/PCA_imputed_tropical_groot.rds")

load_temperate_groot <- extract_load(PCA_temperate_groot$PCA$loadings, 3)
load_continental_groot <- extract_load(PCA_continental_groot$PCA$loadings, 3)
load_tropical_groot <- extract_load(PCA_tropical_groot$PCA$loadings, 3)

long_df_groot <- bind_rows(
  make_df(load_original, "Original", 4),
  make_df(load_temperate_groot,     "Temperate GRooT", 3),
  make_df(load_continental_groot,   "Continental GRooT", 3),
  make_df(load_tropical_groot,     "Tropical GRooT", 3)
)

long_df_groot <- long_df_groot %>%
  group_by(Method, Trait) %>%
  mutate(is_max = abs(value) == max(abs(value))) %>%
  ungroup()


# biomes categorized
PCA_temperate <- readRDS("PCA_results/PCA_imputed_temperate.rds")
PCA_continental <- readRDS("PCA_results/PCA_imputed_continental.rds")
PCA_tropical <- readRDS("PCA_results/PCA_imputed_tropical.rds")
PCA_multiple <- readRDS("PCA_results/PCA_imputed_multiple.rds")
PCA_other <- readRDS("PCA_results/PCA_imputed_other.rds")

load_temperate <- extract_load(PCA_temperate$PCA$loadings, 4)
load_continental <- extract_load(PCA_continental$PCA$loadings, 3)
load_tropical <- extract_load(PCA_tropical$PCA$loadings, 3)
load_multiple <- extract_load(PCA_multiple$PCA$loadings, 4)
load_other <- extract_load(PCA_other$PCA$loadings, 3)

long_df2 <- bind_rows(
  make_df(load_original, "Original", 4),
  make_df(load_temperate,     "Temperate", 4),
  make_df(load_continental,   "Continental", 3),
  make_df(load_tropical,     "Tropical", 3),
  make_df(load_multiple,     "Multiple", 4),
  make_df(load_other,     "Other", 3)
)

long_df2 <- long_df2 %>%
  group_by(Method, Trait) %>%
  mutate(is_max = abs(value) == max(abs(value))) %>%
  ungroup()


# all 9 biomes
PCA_bf <- readRDS("PCA_results/PCA_imputed_bf.rds")
PCA_tgd <- readRDS("PCA_results/PCA_imputed_tgd.rds")
PCA_trf <- readRDS("PCA_results/PCA_imputed_trf.rds")
PCA_troprf <- readRDS("PCA_results/PCA_imputed_troprf.rds")
PCA_tropss <- readRDS("PCA_results/PCA_imputed_tropss.rds")
PCA_tsf <- readRDS("PCA_results/PCA_imputed_tsf.rds")
PCA_ws <- readRDS("PCA_results/PCA_imputed_ws.rds")
PCA_sd <- readRDS("PCA_results/PCA_imputed_sd.rds")

load_bf <- extract_load(PCA_bf$PCA$loadings, 3)
load_tgd <- extract_load(PCA_tgd$PCA$loadings, 4)
load_trf <- extract_load(PCA_trf$PCA$loadings, 3)
load_troprf <- extract_load(PCA_troprf$PCA$loadings, 3)
load_tropss <- extract_load(PCA_tropss$PCA$loadings, 4)
load_tsf <- extract_load(PCA_tsf$PCA$loadings, 4)
load_ws <- extract_load(PCA_ws$PCA$loadings, 4)
load_sd <- extract_load(PCA_sd$PCA$loadings, 3)

long_df_9b <- bind_rows(
  make_df(load_bf, "Boreal Forest", 3),
  make_df(load_tgd,     "Temperate Grassland Desert", 4),
  make_df(load_trf,   "Temperate Rain Forest", 3),
  make_df(load_tsf,     "Temperate Seasonal Forest", 4),
  make_df(load_troprf, "Tropical Rain Forest", 3),
  make_df(load_tropss, "Tropical Seasonal Forest Savanna", 4),
  make_df(load_ws, "Woodland Shrubland", 4),
  make_df(load_sd, "Subtropical Desert", 3)
)

long_df_9b <- long_df_9b %>%
  group_by(Method, Trait) %>%
  mutate(is_max = abs(value) == max(abs(value))) %>%
  ungroup()

# Order factors for layout
long_df1$Trait  <- factor(long_df1$Trait, levels = trait_order)
long_df1$Method <- factor(long_df1$Method, levels = c("Full", "Half", "Quarter", "Tenth"))
long_df1$PC     <- factor(long_df1$PC, levels = paste0("PC", 1:4))

long_df2$Trait  <- factor(long_df2$Trait, levels = trait_order)
long_df2$Method <- factor(long_df2$Method, levels = c("Original", "Temperate", "Continental", "Tropical", "Multiple", "Other"))
long_df2$PC     <- factor(long_df2$PC, levels = paste0("PC", 1:4))

long_df_groot$Trait  <- factor(long_df_groot$Trait, levels = trait_order)
long_df_groot$Method <- factor(long_df_groot$Method, levels = c("Original", "Temperate GRooT", "Continental GRooT", "Tropical GRooT"))
long_df_groot$PC     <- factor(long_df_groot$PC, levels = paste0("PC", 1:4))

long_df_9b$Trait  <- factor(long_df_9b$Trait, levels = trait_order)
long_df_9b$Method <- factor(long_df_9b$Method, levels = c("Temperate Grassland Desert", "Temperate Rain Forest", "Temperate Seasonal Forest", 
                                                          "Tropical Rain Forest", "Tropical Seasonal Forest Savanna", "Subtropical Desert",
                                                          "Boreal Forest", "Woodland Shrubland"))
long_df_9b$PC     <- factor(long_df_9b$PC, levels = paste0("PC", 1:4))

# Labels
long_df1$label <- sprintf("%.2f", long_df1$value)
long_df_groot$label <- sprintf("%.2f", long_df_groot$value)
long_df2$label <- sprintf("%.2f", long_df2$value)
long_df_9b$label <- sprintf("%.2f", long_df_9b$value)

# Color scale: based on uncertainty (SD), only for bootstrap methods
# Use a sequential scale for SD (higher SD = more uncertainty)
col_scale <- scale_fill_distiller(
  palette = "YlOrRd",
  direction = 1,
  na.value = "grey90",
  name = "SD"
)

# Plot
p1 <- ggplot(long_df1, aes(PC, Trait, fill = sd)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label, fontface = ifelse(is_max, "bold", "plain")), size = 3) +
  col_scale +
  scale_y_discrete(limits = rev(trait_order)) +
  facet_wrap(~ Method, ncol = 2) +
  labs(title = "PCA Loadings by Dataset Size",
       subtitle = "Bold = largest |loading| per trait",
       x = "Principal Component",
       y = "Trait") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0)
  )

p2 <- ggplot(long_df_groot, aes(PC, Trait, fill = sd)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label, fontface = ifelse(is_max, "bold", "plain")), size = 3) +
  col_scale +
  scale_y_discrete(limits = rev(trait_order)) +
  facet_wrap(~ Method, ncol = 2) +
  labs(title = "PCA Loadings by GRooT Biome",
       subtitle = "Bold = largest |loading| per trait",
       x = "Principal Component",
       y = "Trait") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0)
  )

p3 <- ggplot(long_df2, aes(PC, Trait, fill = sd)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label, fontface = ifelse(is_max, "bold", "plain")), size = 3) +
  col_scale +
  scale_y_discrete(limits = rev(trait_order)) +
  facet_wrap(~ Method, ncol = 2) +
  labs(title = "PCA Loadings by Biome",
       subtitle = "Bold = largest |loading| per trait",
       x = "Principal Component",
       y = "Trait") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0)
  )

p4 <- ggplot(long_df_9b, aes(PC, Trait, fill = sd)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label, fontface = ifelse(is_max, "bold", "plain")), size = 3) +
  col_scale +
  scale_y_discrete(limits = rev(trait_order)) +
  facet_wrap(~ Method, ncol = 2) +
  labs(title = "PCA Loadings by Biome",
       subtitle = "Bold = largest |loading| per trait",
       x = "Principal Component",
       y = "Trait") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0)
  )

# # Save
dir.create("Figures", showWarnings = FALSE)
ggsave("PCA_plots/Loadings_heatmap_bysize.png", p1, width = 18, height = 16, units = "cm", dpi = 300)
ggsave("PCA_plots/Loadings_heatmap_bybiome_groot.png", p2, width = 18, height = 16, units = "cm", dpi = 300)
ggsave("PCA_plots/Loadings_heatmap_bybiome.png", p3, width = 18, height = 16, units = "cm", dpi = 300)
ggsave("PCA_plots/Loadings_heatmap_bybiome_9b.png", p4, width = 18, height = 16, units = "cm", dpi = 300)
