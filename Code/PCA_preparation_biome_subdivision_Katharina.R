setwd("C:/Users/katha/OneDrive/2_documents/4_Studium/EMDS_Semester3/Capstone")

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
imputed_data_biomes <- imputed_data
imputed_data_biomes$biomesKoeppenGroup <- NA

common_idx <- intersect(rownames(imputed_data_biomes), rownames(rootTraits))
imputed_data_biomes[common_idx, "biomesKoeppenGroup"] <- rootTraits[common_idx, "biomesKoeppenGroup"]

# Group by biome 
imputed_subset_temperate <- subset(imputed_data_biomes, biomesKoeppenGroup == "Temperate") # 348
imputed_subset_temperate_clean <- subset(imputed_subset_temperate, select = -c(biomesKoeppenGroup, genus, family, order))

imputed_subset_tropical <- subset(imputed_data_biomes, biomesKoeppenGroup == "Tropical") # 82
imputed_subset_tropical_clean <- subset(imputed_subset_tropical, select = -c(biomesKoeppenGroup, genus, family, order))

imputed_subset_continental <- subset(imputed_data_biomes, biomesKoeppenGroup == "Continental") # 120
imputed_subset_continental_clean <- subset(imputed_subset_continental, select = -c(biomesKoeppenGroup, genus, family, order))

# 550 data points compared to 1218, loss of 668 data points due to unknown biomes

# =============
# Save as RDS
# =============
fun_save_RDS <- function(file, name){
            output_file <- paste0("PCA_data/", name, ".rds")
            saveRDS(file, output_file)}

# save files
fun_save_RDS(imputed_subset_temperate_clean, "imputed_temperate") # temperate
fun_save_RDS(imputed_subset_tropical_clean, "imputed_tropical") # tropical
fun_save_RDS(imputed_subset_continental_clean, "imputed_continental") # continental

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



