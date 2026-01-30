#setwd("C:/Users/katha/OneDrive/2_documents/4_Studium/EMDS_Semester3/Capstone")

# Import data
rootTraits <- read.table("./data/Root_traits.txt")[, c("SRL", "D", "RTD", "N", "biomesKoeppenGroup")] # already log transformed
aboveTraits <- read.table("./data/Above_traits.txt")
aboveTraits <- log10(aboveTraits) # log transform 
aboveTaxonomy <- read.table("./data/Above_taxonomy.txt")
rootBiomes <- read.table("./data/Root_traits.txt")[, c("biomesKoeppenGroup")]

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
subset_temperate_biome <- subset(AllTraitsAllInfo, biomesKoeppenGroup == "Temperate") # 71
subset_temperate_nobiome <- subset(subset_temperate_biome, select = -c(biomesKoeppenGroup, genus, family, order))
subset_temperate_below <- subset(subset_temperate_nobiome, select = c("SRL", "D", "RTD", "N"))
subset_temperate_above <- subset(subset_temperate_nobiome, select = c("la", "ln", "ph", "sla", "ssd", "sm"))

subset_tropical_biome <- subset(AllTraitsAllInfo, biomesKoeppenGroup == "Tropical") # 27
subset_tropical_nobiome <- subset(subset_tropical_biome, select = -c(biomesKoeppenGroup, genus, family, order))
subset_tropical_below <- subset(subset_tropical_nobiome, select = c("SRL", "D", "RTD", "N"))
subset_tropical_above <- subset(subset_tropical_nobiome, select = c("la", "ln", "ph", "sla", "ssd", "sm"))

subset_continental_biome <- subset(AllTraitsAllInfo, biomesKoeppenGroup == "Continental") # 40
subset_continental_nobiome <- subset(subset_continental_biome, select = -c(biomesKoeppenGroup, genus, family, order))
subset_continental_below <- subset(subset_continental_nobiome, select = c("SRL", "D", "RTD", "N"))
subset_continental_above <- subset(subset_continental_nobiome, select = c("la", "ln", "ph", "sla", "ssd", "sm"))

# 139 data points compared to 301 loss of 160 data points due to unknown biomes

# ================
# WITH IMPUTATION
# ================
imputed_data <- readRDS("./data/imputedTraits.rds")
imputed_data_clean <- subset(imputed_data, select = -c(genus, family, order))
imputed_data_biomes <- imputed_data
imputed_data_biomes$biomesKoeppenGroup <- NA

common_idx <- intersect(rownames(imputed_data_biomes), rownames(rootTraits))
imputed_data_biomes[common_idx, "biomesKoeppenGroup"] <- rootTraits[common_idx, "biomesKoeppenGroup"]

# Group by biome 
# mit imputed, wie viele Datenpunkte pro biom?
imputed_subset_temperate_biome <- subset(imputed_data_biomes, biomesKoeppenGroup == "Temperate") # 348
imputed_subset_temperate_nobiome <- subset(imputed_subset_temperate_biome, select = -c(biomesKoeppenGroup, genus, family, order))
imputed_subset_temperate_below <- subset(imputed_subset_temperate_nobiome, select = c("SRL", "D", "RTD", "N"))
imputed_subset_temperate_above <- subset(imputed_subset_temperate_nobiome, select = c("la", "ln", "ph", "sla", "ssd", "sm"))


imputed_subset_tropical_biome <- subset(imputed_data_biomes, biomesKoeppenGroup == "Tropical") # 82
imputed_subset_tropical_nobiome <- subset(imputed_subset_tropical_biome, select = -c(biomesKoeppenGroup, genus, family, order))
imputed_subset_tropical_below <- subset(imputed_subset_tropical_nobiome, select = c("SRL", "D", "RTD", "N"))
imputed_subset_tropical_above <- subset(imputed_subset_tropical_nobiome, select = c("la", "ln", "ph", "sla", "ssd", "sm"))


imputed_subset_continental_biome <- subset(imputed_data_biomes, biomesKoeppenGroup == "Continental") # 120
imputed_subset_continental_nobiome <- subset(imputed_subset_continental_biome, select = -c(biomesKoeppenGroup, genus, family, order))
imputed_subset_continental_below <- subset(imputed_subset_continental_nobiome, select = c("SRL", "D", "RTD", "N"))
imputed_subset_continental_above <- subset(imputed_subset_continental_nobiome, select = c("la", "ln", "ph", "sla", "ssd", "sm"))

# 550 data points compared to 1218, loss of 668 data points due to unknown biomes

# =============
# Save as RDS
# =============
fun_save_RDS <- function(file, name){
            output_file <- paste0("PCA_data/", name, ".rds")
            saveRDS(file, output_file)}

# temperate
fun_save_RDS(subset_temperate_above, "temperate_above")
fun_save_RDS(subset_temperate_below, "temperate_below")
fun_save_RDS(imputed_subset_temperate_above, "imputed_temperate_above")
fun_save_RDS(imputed_subset_temperate_below, "imputed_temperate_below")

# tropical
fun_save_RDS(subset_tropical_above, "tropical_above")
fun_save_RDS(subset_tropical_below, "tropical_below")
fun_save_RDS(imputed_subset_tropical_above, "imputed_tropical_above")
fun_save_RDS(imputed_subset_tropical_below, "imputed_tropical_below")

# continental
fun_save_RDS(subset_continental_above, "continental_above")
fun_save_RDS(subset_continental_below, "continental_below")
fun_save_RDS(imputed_subset_continental_above, "imputed_continental_above")
fun_save_RDS(imputed_subset_continental_below, "imputed_continental_below")

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
fun_save_RDS(imputed_data_clean, "imputed_combined")
fun_save_RDS(imputed_data_half, "imputed_combined_half")
fun_save_RDS(imputed_data_quarter, "imputed_combined_quarter")
fun_save_RDS(imputed_data_tenth, "imputed_combined_tenth")

# repetitions paper
imputed_below <- subset(imputed_data_clean, select = -c(la, ln, ph, sla, ssd, sm))
imputed_above <- subset(imputed_data_clean, select = -c(SRL, D, RTD, N)) # imputed data here less than not imputed data 
fun_save_RDS(imputed_below, "imputed_below")
fun_save_RDS(imputed_above, "imputed_above")

# these 3 could serve as benchline comparison
fun_save_RDS(AllTraits_nobiomes, "not_imputed_combined")
fun_save_RDS(rootTraitsComplete_nobiomes, "not_imputed_below")
fun_save_RDS(aboveTraitsComplete, "not_imputed_above")

