# Import packages
library(dplyr)
library(psych)
library(ks)
library(TPD)
library(paran)

setwd("C:/Users/katha/OneDrive/2_documents/4_Studium/EMDS_Semester3/Capstone")

# === conduct PCA ===
# function to perform PCA with error handling
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

gridSize <- 50  # Reduced from 200 to save memory (40% reduction)

# 1. PCA imputed combined full
imputed_combined <- readRDS("./PCA_data/imputed_combined.rds")
PCA_imputed_combined_full <- perform_PCA(imputed_combined, "PCA_imputed_combined", gridSize)
rm(PCA_imputed_combined_full); gc()  # Free memory

# 2. PCA imputed combined half
imputed_combined_half <- readRDS("./PCA_data/imputed_combined_half.rds")
PCA_imputed_combined_half <- perform_PCA(imputed_combined_half, "PCA_imputed_combined_half", gridSize)
rm(PCA_imputed_combined_half); gc()  # Free memory

# 3. PCA imputed combined quarter
imputed_combined_quarter <- readRDS("./PCA_data/imputed_combined_quarter.rds")
PCA_imputed_combined_quarter <- perform_PCA(imputed_combined_quarter, "PCA_imputed_combined_quarter", gridSize)
rm(PCA_imputed_combined_quarter); gc()  # Free memory

# 4. PCA imputed combined tenth
imputed_combined_tenth <- readRDS("./PCA_data/imputed_combined_tenth.rds")
PCA_imputed_combined_tenth <- perform_PCA(imputed_combined_tenth, "PCA_imputed_combined_tenth", gridSize)
rm(PCA_imputed_combined_tenth); gc()  # Free memory

# 5. PCA imputed below
imputed_below <- readRDS("./PCA_data/imputed_below.rds")
PCA_imputed_below <- perform_PCA(imputed_below, "PCA_imputed_below", gridSize)
rm(PCA_imputed_below); gc()  # Free memory

# 6. PCA imputed above
imputed_above <- readRDS("./PCA_data/imputed_above.rds")
PCA_imputed_above <- perform_PCA(imputed_above, "PCA_imputed_above", gridSize)
rm(PCA_imputed_above); gc()  # Free memory

# 7. continental below and above
imputed_continental_below <- readRDS("./PCA_data/imputed_continental_below.rds")
PCA_imputed_continental_below <- perform_PCA(imputed_continental_below, "PCA_imputed_continental_below", gridSize)
rm(PCA_imputed_continental_below); gc()  # Free memory

imputed_continental_above <- readRDS("./PCA_data/imputed_continental_above.rds")
PCA_imputed_continental_above <- perform_PCA(imputed_continental_above, "PCA_imputed_continental_above", gridSize)
rm(PCA_imputed_continental_above); gc()  # Free memory

# 8. temperate below and above
imputed_temperate_below <- readRDS("./PCA_data/imputed_temperate_below.rds")
PCA_imputed_temperate_below <- perform_PCA(imputed_temperate_below, "PCA_imputed_temperate_below", gridSize)
rm(PCA_imputed_temperate_below); gc()  # Free memory

imputed_temperate_above <- readRDS("./PCA_data/imputed_temperate_above.rds")
PCA_imputed_temperate_above <- perform_PCA(imputed_temperate_above, "PCA_imputed_temperate_above", gridSize)
rm(PCA_imputed_temperate_above); gc()  # Free memory

# 9. tropical below and above
imputed_tropical_below <- readRDS("./PCA_data/imputed_tropical_below.rds")
PCA_imputed_tropical_below <- perform_PCA(imputed_tropical_below, "PCA_imputed_tropical_below", gridSize)
rm(PCA_imputed_tropical_below); gc()  # Free memory

imputed_tropical_above <- readRDS("./PCA_data/imputed_tropical_above.rds")
PCA_imputed_tropical_above <- perform_PCA(imputed_tropical_above, "PCA_imputed_tropical_above", gridSize)
rm(PCA_imputed_tropical_above); gc()  # Free memory


# ==== NOT IMPUTED ====
# 1. PCA combined
not_imputed_combined <- readRDS("./PCA_data/not_imputed_combined.rds")
PCA_combined <- perform_PCA(not_imputed_combined, "PCA_combined", gridSize)
rm(PCA_combined); gc()  # Free memory

# 2. PCA below
not_imputed_below <- readRDS("./PCA_data/not_imputed_below.rds")
PCA_below <- perform_PCA(not_imputed_below, "PCA_below", gridSize)
rm(PCA_below); gc()  # Free memory

# 3. PCA above
not_imputed_above <- readRDS("./PCA_data/not_imputed_above.rds")
PCA_above <- perform_PCA(not_imputed_above, "PCA_above", gridSize)
rm(PCA_above); gc()  # Free memory

# 4. PCA per biome
# temperate above and below
temperate_above <- readRDS("./PCA_data/temperate_above.rds")
PCA_temperate_above <- perform_PCA(temperate_above, "PCA_temperate_above", gridSize)
rm(PCA_temperate_above); gc()  # Free memory

temperate_below <- readRDS("./PCA_data/temperate_below.rds")
PCA_temperate_below <- perform_PCA(temperate_below, "PCA_temperate_below", gridSize)
rm(PCA_temperate_below); gc()  # Free memory

# tropical
tropical_above <- readRDS("./PCA_data/tropical_above.rds")
PCA_tropical_above <- perform_PCA(tropical_above, "PCA_tropical_above", gridSize)
rm(PCA_tropical_above); gc()  # Free memory

tropical_below <- readRDS("./PCA_data/tropical_below.rds")
PCA_tropical_below <- perform_PCA(tropical_below, "PCA_tropical_below", gridSize)
rm(PCA_tropical_below); gc()  # Free memory


# continental
continental_above <- readRDS("./PCA_data/continental_above.rds")
PCA_continental_above <- perform_PCA(continental, "PCA_continental_above", gridSize)
rm(PCA_continental_above); gc()  # Free memory

continental_below <- readRDS("./PCA_data/continental_below.rds")
PCA_continental_below <- perform_PCA(continental_below, "PCA_continental_below", gridSize)
rm(PCA_continental_below); gc()  # Free memory





