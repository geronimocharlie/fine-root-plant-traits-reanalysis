setwd("C:/Users/katha/OneDrive/2_documents/4_Studium/EMDS_Semester3/Capstone")

# Import data
rootTraits <- read.table("./Data/Root_traits.txt")[, c("SRL", "D", "RTD", "N", "biomesKoeppenGroup")] # already log transformed
aboveTraits <- read.table("./Data/Above_traits.txt")
aboveTraits <- log10(aboveTraits) # log transform 
aboveTaxonomy <- read.table("./Data/Above_taxonomy.txt")
rootBiomes <- read.table("./Data/Root_traits.txt")[, c("biomesKoeppenGroup")]

# Clean data/ omit NAs
rootTraitsComplete <- na.omit(rootTraits) # Species with complete information below ground (748) # with known biomes (380)
rootTraitsComplete_nobiomes <- subset(rootTraitsComplete, select = -(biomesKoeppenGroup))
rootTraitsCompleteBiomes <- rootTraits[complete.cases(rootTraits[, c("SRL", "D", "RTD", "N")]), ]
aboveTraitsCompleteRows <- complete.cases(aboveTraits)
aboveTraitsComplete <- aboveTraits[aboveTraitsCompleteRows, ] # Species with complete information above-ground (2630)

# ===================
# WITHOUT IMPUTATION
# ===================

# Combine data sets
taxonomy <- aboveTaxonomy[aboveTraitsCompleteRows, ]
aboveInRoots <- intersect(rownames(aboveTraitsComplete), rownames(rootTraitsCompleteBiomes))
plantTraitsAbove <- aboveTraitsComplete[aboveInRoots, ]
plantTraitsRoots <- rootTraitsCompleteBiomes[aboveInRoots, ][rownames(plantTraitsAbove), ]
identical(rownames(plantTraitsRoots), rownames(plantTraitsAbove))
AllTraits <- cbind(plantTraitsAbove, plantTraitsRoots)
AllTraits_nobiomes <- subset(AllTraits, select = -(biomesKoeppenGroup))
# add taxonomy info
AllTraitsAllInfoTax <- taxonomy[rownames(AllTraits),]
AllTraitsAllInfo <- cbind(AllTraits, AllTraitsAllInfoTax) # Species with complete information both above- and belowground (301)

# Group by biome 
subset_temperate <- subset(AllTraitsAllInfo, biomesKoeppenGroup == "Temperate") # 71
subset_temperate_clean <- subset(subset_temperate, select = -c(biomesKoeppenGroup, genus, family, order))

subset_tropical <- subset(AllTraitsAllInfo, biomesKoeppenGroup == "Tropical") # 27
subset_tropical_clean <- subset(subset_tropical, select = -c(biomesKoeppenGroup, genus, family, order))

subset_continental <- subset(AllTraitsAllInfo, biomesKoeppenGroup == "Continental") # 40
subset_continental_clean <- subset(subset_continental, select = -c(biomesKoeppenGroup, genus, family, order))

# 139 data points compared to 301 loss of 160 data points due to unknown biomes

# ================
# WITH IMPUTATION
# ================
imputed_data <- readRDS("./data/imputedTraits.rds")
imputed_data_clean <- subset(imputed_data, select = -c(genus, family, order))

# sorting by root biomes
imputed_data_rootbiomes <- imputed_data_clean
imputed_data_rootbiomes$biomesKoeppenGroup <- NA

common_idx <- intersect(rownames(imputed_data_rootbiomes), rownames(rootTraits))
imputed_data_rootbiomes[common_idx, "biomesKoeppenGroup"] <- rootTraits[common_idx, "biomesKoeppenGroup"]

# sorting by full biome groups
biomes <- read.table("./data/biomes_1218sp.txt")
imputed_data_clean$biomes <- NA

common_idx <- intersect(rownames(imputed_data_clean), rownames(biomes))
imputed_data_clean[common_idx, "biomes"] <- apply(
  biomes[common_idx, ],
  1,
  function(x) {
    biome_names <- colnames(biomes)[which(x == 1)]
    paste(biome_names, collapse = ", ")
  }
)

temperate <- c("Temperate.grassland.desert",
               "Temperate.rain.forest",
               "Temperate.seasonal.forest")

tropical <- c("Tropical.rain.forest",
              "Tropical.seasonal.forest.savanna")

continental <- c("Boreal.forest",
                 "Tundra")

# Neue Spalte
imputed_data_clean$biome_short <- sapply(
  imputed_data_clean$biomes,
  function(b) {
    
    if (is.na(b)) return(NA)
    
    # einzelne Biome extrahieren
    biome_vec <- trimws(unlist(strsplit(b, ",")))
    
    # Zugehörigkeit prüfen
    groups <- c(
      if (any(biome_vec %in% temperate)) "Temperate",
      if (any(biome_vec %in% tropical)) "Tropical",
      if (any(biome_vec %in% continental)) "Continental"
    )
    
    if (length(groups) > 1) return("Multiple")
    if (length(groups) == 1) return(groups)
    return("Other")
  }
)

# Group by root biome 
imputed_subset_temperate <- subset(imputed_data_rootbiomes, biomesKoeppenGroup == "Temperate") # 348
imputed_subset_temperate_clean <- subset(imputed_subset_temperate, select = -c(biomesKoeppenGroup))

imputed_subset_tropical <- subset(imputed_data_rootbiomes, biomesKoeppenGroup == "Tropical") # 82
imputed_subset_tropical_clean <- subset(imputed_subset_tropical, select = -c(biomesKoeppenGroup))

imputed_subset_continental <- subset(imputed_data_rootbiomes, biomesKoeppenGroup == "Continental") # 120
imputed_subset_continental_clean <- subset(imputed_subset_continental, select = -c(biomesKoeppenGroup))

# 550 data points compared to 1218, loss of 668 data points due to unknown biomes


# Group by full biomes but categorized
imputed_fullb_subset_temperate <- subset(imputed_data_clean, biome_short == "Temperate") # 707
imputed_fullb_subset_temperate_clean <- subset(imputed_fullb_subset_temperate, select = -c(biomes, biome_short))

imputed_fullb_subset_tropical <- subset(imputed_data_clean, biome_short == "Tropical") # 107
imputed_fullb_subset_tropical_clean <- subset(imputed_fullb_subset_tropical, select = -c(biomes, biome_short))

imputed_fullb_subset_continental <- subset(imputed_data_clean, biome_short == "Continental") # 11
imputed_fullb_subset_continental_clean <- subset(imputed_fullb_subset_continental, select = -c(biomes, biome_short))

imputed_fullb_subset_multiple <- subset(imputed_data_clean, biome_short == "Multiple") # 347
imputed_fullb_subset_multiple_clean <- subset(imputed_fullb_subset_multiple, select = -c(biomes, biome_short))

imputed_fullb_subset_other <- subset(imputed_data_clean, biome_short == "Other") # 46
imputed_fullb_subset_other_clean <- subset(imputed_fullb_subset_other, select = -c(biomes, biome_short))

# Group by all 9 biomes
# Group by full biomes but categorized
imputed_subset_bf <-subset(imputed_data_clean,grepl("Boreal\\.forest", biomes)) # 169
imputed_subset_bf_clean <- subset(imputed_subset_bf, select = -c(biomes, biome_short))

imputed_subset_tun <-subset(imputed_data_clean,grepl("Tundra", biomes)) # 14
imputed_subset_tun_clean <- subset(imputed_subset_tun, select = -c(biomes, biome_short))

imputed_subset_tgd <-subset(imputed_data_clean,grepl("Temperate\\.grassland\\.desert", biomes)) # 244
imputed_subset_tgd_clean <- subset(imputed_subset_tgd, select = -c(biomes, biome_short))

imputed_subset_trf <-subset(imputed_data_clean,grepl("Temperate\\.rain\\.forest", biomes)) # 139
imputed_subset_trf_clean <- subset(imputed_subset_trf, select = -c(biomes, biome_short))

imputed_subset_tsf <-subset(imputed_data_clean,grepl("Temperate\\.seasonal\\.forest", biomes)) # 885
imputed_subset_tsf_clean <- subset(imputed_subset_tsf, select = -c(biomes, biome_short))

imputed_subset_troprf <-subset(imputed_data_clean,grepl("Tropical\\.rain\\.forest", biomes)) # 94
imputed_subset_troprf_clean <- subset(imputed_subset_troprf, select = -c(biomes, biome_short))

imputed_subset_tropss <-subset(imputed_data_clean,grepl("Tropical\\.seasonal\\.forest\\.savanna", biomes)) # 294
imputed_subset_tropss_clean <- subset(imputed_subset_tropss, select = -c(biomes, biome_short))

imputed_subset_ws <-subset(imputed_data_clean,grepl("Woodland\\.shrubland", biomes)) # 861
imputed_subset_ws_clean <- subset(imputed_subset_ws, select = -c(biomes, biome_short))

imputed_subset_sd <-subset(imputed_data_clean,grepl("Subtropical\\.desert", biomes)) # 95
imputed_subset_sd_clean <- subset(imputed_subset_sd, select = -c(biomes, biome_short))

# =============
# Save as RDS
# =============
fun_save_RDS <- function(file, name){
            output_file <- paste0("Data/", name, ".rds")
            saveRDS(file, output_file)}

# save files root biomes
fun_save_RDS(imputed_subset_temperate_clean, "imputed_temperate") # temperate
fun_save_RDS(imputed_subset_tropical_clean, "imputed_tropical") # tropical
fun_save_RDS(imputed_subset_continental_clean, "imputed_continental") # continental

# save files full biomes categorized
fun_save_RDS(imputed_fullb_subset_temperate_clean, "imputed_temperate_fullb") # temperate
fun_save_RDS(imputed_fullb_subset_tropical_clean, "imputed_tropical_fullb") # tropical
fun_save_RDS(imputed_fullb_subset_continental_clean, "imputed_continental_fullb") # continental
fun_save_RDS(imputed_fullb_subset_multiple_clean, "imputed_multiple_fullb") # continental
fun_save_RDS(imputed_fullb_subset_other_clean, "imputed_other_fullb") # continental

# save files all 9 biomes
fun_save_RDS(imputed_subset_bf_clean, "imputed_bf")
fun_save_RDS(imputed_subset_tun_clean, "imputed_tun")
fun_save_RDS(imputed_subset_tgd_clean, "imputed_tgd")
fun_save_RDS(imputed_subset_trf_clean, "imputed_trf")
fun_save_RDS(imputed_subset_tsf_clean, "imputed_tsf")
fun_save_RDS(imputed_subset_troprf_clean, "imputed_troprf")
fun_save_RDS(imputed_subset_tropss_clean, "imputed_tropss")
fun_save_RDS(imputed_subset_ws_clean, "imputed_ws")
fun_save_RDS(imputed_subset_sd_clean, "imputed_sd")

# for PCA stability check
# PCA combined für 1218 imputed data points, dann 1/2 (609) der Daten, 1/4 (304) der Daten, 1/10 (122) der Daten
n <- nrow(imputed_data_clean)
idx_half <- sample(seq_len(n), size = floor(n / 2), replace = FALSE)
idx_quarter <- sample(seq_len(n), size = floor(n / 4), replace = FALSE)
idx_tenth <- sample(seq_len(n), size = floor(n / 10), replace = FALSE)
imputed_data_half <- imputed_data_clean[idx_half, ] # select random half of data
imputed_data_quarter <- imputed_data_clean[idx_quarter, ] # quarter
imputed_data_tenth <- imputed_data_clean[idx_tenth, ] # tenth

# save files
fun_save_RDS(imputed_data_clean, "imputed_combined") # recreating paper
fun_save_RDS(imputed_data_half, "imputed_combined_half") # for stability test
fun_save_RDS(imputed_data_quarter, "imputed_combined_quarter")
fun_save_RDS(imputed_data_tenth, "imputed_combined_tenth")




