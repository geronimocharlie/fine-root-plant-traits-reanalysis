# Import packages
library(dplyr)
library(psych)
library(ks)
library(TPD)
library(paran)

#setwd("C:/Users/katha/OneDrive/2_documents/4_Studium/EMDS_Semester3/Capstone")

# =====================================
# CHECK FILES AND CREATE OUTPUT DIRS
# =====================================
if (!file.exists("data/Root_traits.txt")) {
  stop("File not found: data/Root_traits.txt")
}
if (!file.exists("data/Above_traits.txt")) {
  stop("File not found: data/Above_traits.txt")
}
if (!file.exists("data/Above_taxonomy.txt")) {
  stop("File not found: data/Above_taxonomy.txt")
}
if (!file.exists("data/imputedTraits.rds")) {
  stop("File not found: data/imputedTraits.rds")
}

if (!dir.exists("PCA_results")) {
  dir.create("PCA_results", recursive = TRUE)
}


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

subset_tropical_biome <- subset(AllTraitsAllInfo, biomesKoeppenGroup == "Tropical") # 27
subset_tropical_nobiome <- subset(subset_tropical_biome, select = -c(biomesKoeppenGroup, genus, family, order))

subset_continental_biome <- subset(AllTraitsAllInfo, biomesKoeppenGroup == "Continental") # 40
subset_continental_nobiome <- subset(subset_continental_biome, select = -c(biomesKoeppenGroup, genus, family, order))
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

imputed_subset_tropical_biome <- subset(imputed_data_biomes, biomesKoeppenGroup == "Tropical") # 82
imputed_subset_tropical_nobiome <- subset(imputed_subset_tropical_biome, select = -c(biomesKoeppenGroup, genus, family, order))

imputed_subset_continental_biome <- subset(imputed_data_biomes, biomesKoeppenGroup == "Continental") # 120
imputed_subset_continental_nobiome <- subset(imputed_subset_continental_biome, select = -c(biomesKoeppenGroup, genus, family, order))
# 550 data points compared to 1218, loss of 668 data points due to unknown biomes

# ===============
# PCAs & Plots
# ===============

# ==== IMPUTED ====
# 1. PCA combined für 1218 imputed data points, dann 1/2 (609) der Daten, 1/4 (304) der Daten, 1/10 (122) der Daten
# -> wie sehr ändern sich die PCAs sind sie bei wenig Daten stabil?
# wenn immer andere PCAs rauskommen für verschiedene 1/4 dann kann man sich biome PCAs auch sparen
n <- nrow(imputed_data_clean)
idx_half <- sample(seq_len(n), size = floor(n / 2), replace = FALSE)
idx_quarter <- sample(seq_len(n), size = floor(n / 4), replace = FALSE)
idx_tenth <- sample(seq_len(n), size = floor(n / 10), replace = FALSE)

# imputed_data_clean
imputed_data_half <- imputed_data_clean[idx_half, ]
imputed_data_quarter <- imputed_data_clean[idx_quarter, ]
imputed_data_tenth <- imputed_data_clean[idx_tenth, ]

# === PCA workflow to check stability ===
# Helper function to perform PCA with error handling
perform_PCA <- function(traits_df, pca_name, gridSize = 100) {
  tryCatch({
    cat("\n=== Computing", pca_name, "===\n")
    flush.console()
    
    pca_list <- list()
    pca_list$traits <- traits_df
    pca_list$dimensions <- paran(pca_list$traits)$Retained
    pca_list$means <- apply(pca_list$traits, 2, mean)
    pca_list$sds <- apply(pca_list$traits, 2, sd)
    pca_list$PCA <- psych::principal(scale(pca_list$traits), nfactors=pca_list$dimensions, 
                                     rotate="varimax", covar = FALSE)  # covar = FALSE to save RAM
    pca_list$Variance <- pca_list$PCA$Vaccounted[2,]
    sqrtEigen <- sqrt(colSums(pca_list$PCA$loadings**2))
    for(i in 1:pca_list$dimensions){
      pca_list$PCA$scores[, i] <- pca_list$PCA$scores[, i] * sqrtEigen[i] 
    }
    sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(pca_list$PCA$loadings)), byrow=T, 
                           nrow = nrow(pca_list$PCA$loadings))
    pca_list$PCA$loadings <- pca_list$PCA$loadings / sqrtEigenMat
    pca_list$traitsUse <- data.frame(pca_list$PCA$scores[, 1:pca_list$dimensions])
    sdTraits <- sqrt(diag(Hpi.diag(pca_list$traitsUse)))
    pca_list$TPDs <- TPDsMean(species = rownames(pca_list$traitsUse), 
                              means = pca_list$traitsUse, 
                              sds = matrix(rep(sdTraits, nrow(pca_list$traitsUse)), byrow=T, 
                                           ncol=pca_list$dimensions),
                              n_divisions = gridSize)
    
    output_file <- paste0("PCA_results/", pca_name, ".rds")
    saveRDS(pca_list, output_file)
    cat("Successfully saved:", output_file, "\n")
    return(pca_list)
    
  }, error = function(e) {
    cat("ERROR in", pca_name, ":", e$message, "\n")
    cat("Skipping this analysis and continuing...\n")
    return(NULL)
  })
}

# 1. PCA imputed combined full
gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
PCA_combined_full <- perform_PCA(imputed_data_clean, "PCA_imputed_combined_full", gridSize)
rm(PCA_combined_full); gc()  # Free memory

# 2. PCA imputed combined half
gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
PCA_combined_full <- perform_PCA(imputed_data_half, "PCA_imputed_combined_half", gridSize)
rm(PCA_combined_full); gc()  # Free memory

# 3. PCA imputed combined quarter
gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
PCA_combined_full <- perform_PCA(imputed_data_quarter, "PCA_imputed_combined_quarter", gridSize)
rm(PCA_combined_full); gc()  # Free memory

# 4. PCA imputed combined tenth
gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
PCA_combined_full <- perform_PCA(imputed_data_tenth, "PCA_imputed_combined_tenth", gridSize)
rm(PCA_combined_full); gc()  # Free memory

#  now compare if PCA is stable, if not, no need to look at biome specific PCAs

# 5. PCA below
#imputed_below <- subset(imputed_data_clean, select = -c(la, ln, ph, sla, ssd, sm))

#gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
#PCA_combined_full <- perform_PCA(imputed_below, "PCA_imputed_below", gridSize)
#rm(PCA_combined_full); gc()  # Free memory

# 6. PCA above
#imputed_above <- subset(imputed_data_clean, select = -c(SRL, D, RTD, N))# imputed data here less than not imputed data 

#gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
#PCA_combined_full <- perform_PCA(imputed_above, "PCA_imputed_above", gridSize)
#rm(PCA_combined_full); gc()  # Free memory


# PCA für imputed biomes (seperated), falls PCA stabil
# 7. continental
#gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
#PCA_combined_full <- perform_PCA(imputed_subset_continental_nobiome, "PCA_imputed_continental", gridSize)
#rm(PCA_combined_full); gc()  # Free memory

# 8. temperate
#gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
#PCA_combined_full <- perform_PCA(imputed_subset_temperate_nobiome, "PCA_imputed_temperate", gridSize)
#rm(PCA_combined_full); gc()  # Free memory

# 9. tropical
#gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
#PCA_combined_full <- perform_PCA(imputed_subset_tropical_nobiome, "PCA_imputed_tropical", gridSize)
#rm(PCA_combined_full); gc()  # Free memory


# ==== NOT IMPUTED ====
# 1. PCA combined
gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
PCA_combined_full <- perform_PCA(AllTraits_nobiomes, "PCA_combined", gridSize)
rm(PCA_combined_full); gc()  # Free memory

# 2. PCA below
#gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
#PCA_combined_full <- perform_PCA(rootTraitsComplete_nobiomes, "PCA_below", gridSize)
#rm(PCA_combined_full); gc()  # Free memory

# 3. PCA above
#gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
#PCA_combined_full <- perform_PCA(aboveTraitscomplete, "PCA_above", gridSize)
#rm(PCA_combined_full); gc()  # Free memory

# 4. PCA per biome
# temperate
#gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
#PCA_combined_full <- perform_PCA(subset_temperate_nobiome, "PCA_temperate", gridSize)
#rm(PCA_combined_full); gc()  # Free memory

# tropical
#gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
#PCA_combined_full <- perform_PCA(subset_tropical_nobiome, "PCA_tropical", gridSize)
#rm(PCA_combined_full); gc()  # Free memory

# continental
gridSize <- 100  # Reduced from 200 to save memory (40% reduction)
PCA_combined_full <- perform_PCA(subset_continental_nobiome, "PCA_continental", gridSize)
rm(PCA_combined_full); gc()  # Free memory




