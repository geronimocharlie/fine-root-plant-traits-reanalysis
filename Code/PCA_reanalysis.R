# Import packages
library(dplyr)
library(psych)
library(ks)
library(TPD)
library(paran)

# === conduct PCA ===
# function to perform PCA with error handling
perform_PCA <- function(traits_df, pca_name, 
                        gridSize_full = 30,
                        gridSize_2D = 200) {
  tryCatch({
    cat("\n=== Computing", pca_name, "===\n")
    flush.console()
    
    pca_list <- list()
    pca_list$traits <- traits_df
    pca_list$dimensions <- paran(pca_list$traits)$Retained
    pca_list$dimensions <- min(paran(pca_list$traits)$Retained, 4) 
    pca_list$means <- apply(pca_list$traits, 2, mean)
    pca_list$sds <- apply(pca_list$traits, 2, sd)
    
    # pca computation
    pca_list$PCA <- psych::principal(scale(pca_list$traits), nfactors=pca_list$dimensions, 
                                     rotate="varimax", covar = TRUE) 
    pca_list$Variance <- pca_list$PCA$Vaccounted[2,]
    
    # eigenvalue extraction
    sqrtEigen <- sqrt(colSums(pca_list$PCA$loadings**2))
    
    # score rescaling
    for(i in 1:pca_list$dimensions){
      pca_list$PCA$scores[, i] <- pca_list$PCA$scores[, i] * sqrtEigen[i] 
    }
    
    # Convert loadings to eigenvectors
    sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(pca_list$PCA$loadings)),
                           byrow = TRUE,
                           nrow = nrow(pca_list$PCA$loadings))
    
    # convert loadings to eigenvectors
    pca_list$PCA$loadings <- pca_list$PCA$loadings / sqrtEigenMat
    
    # axis orientation flipping
    for(i in 1:pca_list$dimensions){
      if(i == 1 & pca_list$PCA$loadings["ph", i] < 0){ #ph is positive
        pca_list$PCA$loadings[, i] <- -1 * pca_list$PCA$loadings[, i]
        pca_list$PCA$scores[, i] <- -1 *pca_list$PCA$scores[, i]
      }
      if(i == 2 & pca_list$PCA$loadings["sla", i] < 0){ #sla is positive
        pca_list$PCA$loadings[, i] <- -1 * pca_list$PCA$loadings[, i]
        pca_list$PCA$scores[, i] <- -1 *pca_list$PCA$scores[, i]
      }
      if(i == 3 & pca_list$PCA$loadings["SRL", i] > 0){ #SRL is negative
        pca_list$PCA$loadings[, i] <- -1 * pca_list$PCA$loadings[, i]
        pca_list$PCA$scores[, i] <- -1 *pca_list$PCA$scores[, i]
      }
      if(i == 4 & pca_list$PCA$loadings["N", i] < 0){ #N is ppositive
        pca_list$PCA$loadings[, i] <- -1 * pca_list$PCA$loadings[, i]
        pca_list$PCA$scores[, i] <- -1 *pca_list$PCA$scores[, i]
      }
    }
    
    # convert discrete PCA scores into probability density functions in trait space
    pca_list$traitsUse <- data.frame(pca_list$PCA$scores[, 1:pca_list$dimensions])
    sdTraits <- sqrt(diag(Hpi.diag(pca_list$traitsUse)))
    
    pca_list$TPDs <- TPDsMean(
      species = rownames(pca_list$traitsUse), 
      means = pca_list$traitsUse, 
      sds = matrix(rep(sdTraits, nrow(pca_list$traitsUse)), byrow=TRUE, 
                   ncol=pca_list$dimensions),
      n_divisions = gridSize_full
    )
    
    # pairwise tpd (trait probability density)
    if (pca_list$dimensions >= 2) {
      
      pca_list$traitsUse2D <- list()
      pca_list$TPDs2D <- list()
      
      combs <- combn(pca_list$dimensions, 2)
      
      for (j in seq_len(ncol(combs))) {
        
        pair <- combs[, j]
        name <- paste0("Comp", pair[1], "_Comp", pair[2])
        
        pca_list$traitsUse2D[[name]] <-
          data.frame(pca_list$PCA$scores[, pair])
        
        sdTraits <- sqrt(diag(
          Hpi.diag(pca_list$traitsUse2D[[name]])
        ))
        
        pca_list$TPDs2D[[name]] <- TPDsMean(
          species = rownames(pca_list$traitsUse2D[[name]]),
          means   = pca_list$traitsUse2D[[name]],
          sds     = matrix(rep(sdTraits,
                               nrow(pca_list$traitsUse2D[[name]])),
                           byrow = TRUE,
                           ncol  = 2),
          n_divisions = gridSize_2D
        )
      }
    }
    
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

gridSize_full <- 30
gridSize_2D <- 200 # decrease to save RAM


# 1. PCA imputed combined full
imputed_combined <- readRDS("./Data/imputed_combined.rds")
PCA_imputed_combined_full <- perform_PCA(imputed_combined, "PCA_imputed_combined", gridSize_full, gridSize_2D)
rm(PCA_imputed_combined_full); gc()  # Free memory

# 2. PCA imputed combined half
imputed_combined_half <- readRDS("./Data/imputed_combined_half.rds")
PCA_imputed_combined_half <- perform_PCA(imputed_combined_half, "PCA_imputed_combined_half", gridSize_full, gridSize_2D)
rm(PCA_imputed_combined_half); gc()  # Free memory

# 3. PCA imputed combined quarter
imputed_combined_quarter <- readRDS("./Data/imputed_combined_quarter.rds")
PCA_imputed_combined_quarter <- perform_PCA(imputed_combined_quarter, "PCA_imputed_combined_quarter", gridSize_full, gridSize_2D)
rm(PCA_imputed_combined_quarter); gc()  # Free memory
 
# 4. PCA imputed combined tenth
imputed_combined_tenth <- readRDS("./Data/imputed_combined_tenth.rds")
PCA_imputed_combined_tenth <- perform_PCA(imputed_combined_tenth, "PCA_imputed_combined_tenth", gridSize_full, gridSize_2D)
rm(PCA_imputed_combined_tenth); gc()  # Free memory

# 7. continental 
imputed_continental <- readRDS("./Data/imputed_continental.rds")
PCA_imputed_continental <- perform_PCA(imputed_continental, "PCA_imputed_continental", gridSize_full, gridSize_2D)
rm(PCA_imputed_continental); gc()  # Free memory

# 8. temperate 
imputed_temperate <- readRDS("./PCA_data/imputed_temperate.rds")
PCA_imputed_temperate <- perform_PCA(imputed_temperate, "PCA_imputed_temperate", gridSize_full, gridSize_2D)
rm(PCA_imputed_temperate); gc()  # Free memory

# 9. tropical 
imputed_tropical <- readRDS("./Data/imputed_tropical.rds")
PCA_imputed_tropical <- perform_PCA(imputed_tropical, "PCA_imputed_tropical", gridSize_full, gridSize_2D)
rm(PCA_imputed_tropical); gc()  # Free memory


# Full biome separation
# 1. Boreal forest
imputed_bf <- readRDS("./Data/imputed_bf.rds")
PCA_imputed_bf <- perform_PCA(imputed_bf, "PCA_imputed_bf", gridSize_full, gridSize_2D)
rm(PCA_imputed_bf); gc()  # Free memory

# 2. Tundra
imputed_tun <- readRDS("./Data/imputed_tun.rds")
PCA_imputed_tun <- perform_PCA(imputed_tun, "PCA_imputed_tun", gridSize_full, gridSize_2D)
rm(PCA_imputed_tun); gc()  # Free memory

# 3. Temperate Grassland Desert
imputed_tgd <- readRDS("./Data/imputed_tgd.rds")
PCA_imputed_tgd <- perform_PCA(imputed_tgd, "PCA_imputed_tgd", gridSize_full, gridSize_2D)
rm(PCA_imputed_tgd); gc()  # Free memory

# 4. Temperate Rain Forest
imputed_trf <- readRDS("./Data/imputed_trf.rds")
PCA_imputed_trf <- perform_PCA(imputed_trf, "PCA_imputed_trf", gridSize_full, gridSize_2D)
rm(PCA_imputed_trf); gc()  # Free memory

# 5. Temperate seasonal forest
imputed_tsf <- readRDS("./Data/imputed_tsf.rds")
PCA_imputed_tsf <- perform_PCA(imputed_tsf, "PCA_imputed_tsf", gridSize_full, gridSize_2D)
rm(PCA_imputed_tsf); gc()  # Free memory

# 6. Tropical Rain Forest
imputed_troprf <- readRDS("./Data/imputed_troprf.rds")
PCA_imputed_troprf <- perform_PCA(imputed_troprf, "PCA_imputed_troprf", gridSize_full, gridSize_2D)
rm(PCA_imputed_troprf); gc()  # Free memory

# 7. Tropical seasonal Forest Savanna
imputed_tropss <- readRDS("./Data/imputed_tropss.rds")
PCA_imputed_tropss <- perform_PCA(imputed_tropss, "PCA_imputed_tropss", gridSize_full, gridSize_2D)
rm(PCA_imputed_tropss); gc()  # Free memory

# 8. Woodland Shrubland
imputed_ws <- readRDS("./Data/imputed_ws.rds")
PCA_imputed_ws <- perform_PCA(imputed_ws, "PCA_imputed_ws", gridSize_full, gridSize_2D)
rm(PCA_imputed_ws); gc()  # Free memory

# 9. Subtropical desert
imputed_sd <- readRDS("./Data/imputed_sd.rds")
PCA_imputed_sd <- perform_PCA(imputed_sd, "PCA_imputed_sd", gridSize_full, gridSize_2D)
rm(PCA_imputed_sd); gc()  # Free memory


# Full biomes categorized
# 1. PCA imputed multiple biomes
imputed_multiple <- readRDS("./Data/imputed_multiple_fullb.rds")
PCA_imputed_multiple <- perform_PCA(imputed_multiple, "PCA_imputed_multiple", gridSize_full, gridSize_2D)
rm(PCA_imputed_mulitple); gc()  # Free memory

# 2. PCA imputed other biomes
imputed_other <- readRDS("./Data/imputed_other_fullb.rds")
PCA_imputed_other <- perform_PCA(imputed_other, "PCA_imputed_other", gridSize_full, gridSize_2D)
rm(PCA_imputed_other); gc()  # Free memory

# 3. continental 
imputed_continental_fullb <- readRDS("./Data/imputed_continental_fullb.rds")
PCA_imputed_continental_fullb <- perform_PCA(imputed_continental_fullb, "PCA_imputed_continental_fullb", gridSize_full, gridSize_2D)
rm(PCA_imputed_continental_fullb); gc()  # Free memory

# 4. temperate
imputed_temperate_fullb <- readRDS("./Data/imputed_temperate_fullb.rds")
PCA_imputed_temperate_fullb <- perform_PCA(imputed_temperate_fullb, "PCA_imputed_temperate_fullb", gridSize_full, gridSize_2D)
rm(PCA_imputed_temperate_fullb); gc()  # Free memory

# 5. tropical 
imputed_tropical_fullb <- readRDS("./Data/imputed_tropical_fullb.rds")
PCA_imputed_tropical_fullb <- perform_PCA(imputed_tropical_fullb, "PCA_imputed_tropical_fullb", gridSize_full, gridSize_2D)
rm(PCA_imputed_tropical_fullb); gc()  # Free memory





