### In this script we prepare plant trait data for analyses: 
# 0. Loading trait data (average of trait values at the species level from GROOT (fine root traits), which are already log-transformed, and loading of aboveground trait data from TRY (needs log-transformation). Taxonomical information for species with aboveground traits comes from The Plant List (package 'Taxonstand').
# 1. Selecting species with complete trait information for aboveground and fine-root traits separately, and for all traits combined.
# 2. Estimation of functional space (including estimation of dimensionality, estimation of PCA and varimax rotation, and estimation of TPD functions) based on complete information for aboveground traits, fine-root traits, and all traits combined.
cat(paste0("\n\n Starting script #1 \n\n"))
########################################################################
### 0. Loading trait and taxonomy data at species level
########################################################################
rootTraits <- read.table("data/Root_traits.txt")[, c("SRL", "D", "RTD", "N")] # Load fine-root traits: Specific Root Length, Diameter, Root Tissue Density, Nitrogen content
aboveTraits <- read.table("data/Above_traits.txt") # Load aboveground traits from TRY database
aboveTraits <- log10(aboveTraits) # Log-transform aboveground traits for normality
aboveTaxonomy <- read.table("data/Above_taxonomy.txt") # Load taxonomic info (genus, family) for species with aboveground traits

########################################################################
### 1. Selection of species with complete empirical information
########################################################################
rootTraitsComplete <- na.omit(rootTraits) # Remove species with any missing fine-root traits → 748 species
aboveTraitsCompleteRows <- complete.cases(aboveTraits) # Logical vector: which species have complete aboveground traits
aboveTraitsComplete <- aboveTraits[aboveTraitsCompleteRows, ] # Extract species with complete aboveground data → 2630 species

taxonomy <- aboveTaxonomy[aboveTraitsCompleteRows, ] # Apply same filter to taxonomy data to keep alignment
aboveInRoots <- intersect(rownames(aboveTraitsComplete), rownames(rootTraitsComplete)) # Find species present in BOTH datasets
plantTraitsAbove <- aboveTraitsComplete[aboveInRoots, ] # Subset aboveground traits to species also having root data
plantTraitsRoots <- rootTraitsComplete[aboveInRoots, ][rownames(plantTraitsAbove), ] # Subset root traits and ensure same order
identical(rownames(plantTraitsRoots), rownames(plantTraitsAbove)) # Verify row names match perfectly
AllTraitsAllInfo <- cbind(plantTraitsAbove, plantTraitsRoots) # Combine above + below traits for species with both
traitsUse <- AllTraitsAllInfo # Copy for processing
AllTraitsAllInfoTax <- taxonomy[rownames(AllTraitsAllInfo),] # Extract taxonomy for the combined dataset
AllTraitsAllInfo <- cbind(traitsUse, AllTraitsAllInfoTax) # Final dataset: traits + taxonomy → 301 species

########################################################################
### 2. Estimation of functional spaces for abovegroud, belowground, and combined
########################################################################

# Define alpha parameter for TPD calculations ##NEW
alphaUse <- 0.95

####################################
#### A. Only aboveground traits ####
gridSize <- 200 # Grid resolution for TPD calculations (200x200 = 40,000 cells)
PCAAbove <- list() # Container for all aboveground PCA results
PCAAbove$traits <- aboveTraitsComplete # Store the trait data used
PCAAbove$dimensions <- paran(PCAAbove$traits)$Retained # Determine optimal number of PCA dimensions using parallel analysis
PCAAbove$means <- apply(PCAAbove$traits, 2, mean) # Calculate trait means for each column
PCAAbove$sds <- apply(PCAAbove$traits, 2, sd) # Calculate trait standard deviations
PCAAbove$PCA <- psych::principal(scale(PCAAbove$traits), nfactors=PCAAbove$dimensions, 
                                 rotate="varimax", covar = T) # Run PCA with varimax rotation on standardized traits
PCAAbove$Variance <- PCAAbove$PCA$Vaccounted[2,] # Extract proportion of variance explained by each component
sqrtEigen <- sqrt(colSums(PCAAbove$PCA$loadings**2)) # Calculate scaling factors from loadings
for(i in 1:PCAAbove$dimensions){
  PCAAbove$PCA$scores[, i] <- PCAAbove$PCA$scores[, i] * sqrtEigen[i] # Scale PCA scores by sqrt(eigenvalues)
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCAAbove$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCAAbove$PCA$loadings)) # Create matrix of scaling factors
PCAAbove$PCA$loadings <- PCAAbove$PCA$loadings / sqrtEigenMat # Convert loadings to unit eigenvectors
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCAAbove$traitsUse <- data.frame(PCAAbove$PCA$scores[, 1:PCAAbove$dimensions]) # Extract species scores in PCA space
sdTraits <- sqrt(diag(Hpi.diag(PCAAbove$traitsUse))) # Calculate bandwidth for kernel density estimation
PCAAbove$TPDs <- TPDsMean(species = rownames(PCAAbove$traitsUse), 
                    means = PCAAbove$traitsUse, 
                    sds = matrix(rep(sdTraits, nrow(PCAAbove$traitsUse)), byrow=T, 
                                 ncol=PCAAbove$dimensions),
                    n_divisions = gridSize) # Calculate Trait Probability Densities for each species
saveRDS(PCAAbove, paste0("data/PCAAboveONLY_CompleteObs.rds")) # Save aboveground PCA results
####################################
#### B. Only fine-root traits   ####
gridSize <- 200 # Same grid resolution as aboveground analysis
PCABelow <- list() # Container for belowground PCA results
plantTraitsBelow <- c("SRL", "D", "RTD", "N") # Define the 4 root trait names
PCABelow$traits <- rootTraitsComplete # Use complete root trait dataset (748 species)
PCABelow$dimensions <- paran(PCABelow$traits)$Retained # Determine optimal dimensions for root traits
PCABelow$means <- apply(PCABelow$traits, 2, mean)
PCABelow$sds <- apply(PCABelow$traits, 2, sd)
PCABelow$PCA <- psych::principal(scale(PCABelow$traits), nfactors=PCABelow$dimensions, 
                                 rotate="varimax", covar = T) # PCA with varimax rotation on root traits
PCABelow$Variance <- PCABelow$PCA$Vaccounted[2,] # Variance explained by each component
sqrtEigen <- sqrt(colSums(PCABelow$PCA$loadings**2)) # Calculate scaling factors
for(i in 1:PCABelow$dimensions){
  PCABelow$PCA$scores[, i] <- PCABelow$PCA$scores[, i] * sqrtEigen[i] # Scale scores
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCABelow$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCABelow$PCA$loadings))
PCABelow$PCA$loadings <- PCABelow$PCA$loadings / sqrtEigenMat # Convert to unit eigenvectors
# Check all are 1: colSums(PCABelow$PCA$loadings**2)
PCABelow$traitsUse <- data.frame(PCABelow$PCA$scores[, 1:PCABelow$dimensions]) # Species positions in root trait space
sdTraits <- sqrt(diag(Hpi.diag(PCABelow$traitsUse))) # Bandwidth for kernel density estimation
PCABelow$TPDs <- TPDsMean(species = rownames(PCABelow$traitsUse), 
                          means = PCABelow$traitsUse, 
                          sds = matrix(rep(sdTraits, nrow(PCABelow$traitsUse)), 
                                       byrow=T, ncol=PCABelow$dimensions),
                          n_divisions = gridSize) # Calculate TPDs for root functional space
saveRDS(PCABelow, paste0("data/PCABelowONLY_CompleteObs.rds")) # Save belowground results 


##### PROCRUSTES ABOVE-BELOW:
commonSP <- intersect(rownames(PCABelow$traitsUse), rownames(PCAAbove$traitsUse)) # Find species common to both trait spaces
vegan::protest(PCABelow$traitsUse[commonSP, ], PCAAbove$traitsUse[commonSP, ]) # Test similarity between above/below functional spaces

####################################
####     C. All traits         ##### 
AllTraits <- AllTraitsAllInfo[, c("la", "ln", "ph", "sla", "ssd", "sm",
                                  "SRL", "D", "RTD", "N")] # Extract 10 traits: 6 above + 4 below
gridSize <- 30 # Smaller grid for 4D space (30^4 = 810,000 cells vs 200^2 = 40,000 for 2D)
PCATotal <- list() # Container for combined trait analysis
PCATotal$traits <- AllTraits
PCATotal$dimensions <- paran(PCATotal$traits)$Retained # Determine dimensions via parallel analysis
PCATotal$means <- apply(PCATotal$traits, 2, mean)
PCATotal$sds <- apply(PCATotal$traits, 2, sd)
PCATotal$AllInfo <- AllTraitsAllInfo # Store complete dataset including taxonomy
PCATotal$PCA <- psych::principal(scale(PCATotal$traits), nfactors=PCATotal$dimensions, 
                                 rotate="varimax", covar = T) # PCA on all 10 traits with varimax rotation
PCATotal$Variance <- PCATotal$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCATotal$PCA$loadings**2))
for(i in 1:PCATotal$dimensions){
  PCATotal$PCA$scores[, i] <- PCATotal$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCATotal$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCATotal$PCA$loadings))
PCATotal$PCA$loadings <- PCATotal$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCATotal$PCA$loadings**2)
#### Make sure orientation of all vectors is always consistent:
for(i in 1:PCATotal$dimensions){
  if(i == 1 & PCATotal$PCA$loadings["ph", i] < 0){ # Force C1 to have positive plant height loading
    PCATotal$PCA$loadings[, i] <- -1 * PCATotal$PCA$loadings[, i]
    PCATotal$PCA$scores[, i] <- -1 *PCATotal$PCA$scores[, i]
  }
  if(i == 2 & PCATotal$PCA$loadings["sla", i] < 0){ # Force C2 to have positive SLA loading
    PCATotal$PCA$loadings[, i] <- -1 * PCATotal$PCA$loadings[, i]
    PCATotal$PCA$scores[, i] <- -1 *PCATotal$PCA$scores[, i]
  }
  if(i == 3 & PCATotal$PCA$loadings["SRL", i] > 0){ # Force C3 to have negative SRL loading
    PCATotal$PCA$loadings[, i] <- -1 * PCATotal$PCA$loadings[, i]
    PCATotal$PCA$scores[, i] <- -1 *PCATotal$PCA$scores[, i]
  }
  if(i == 4 & PCATotal$PCA$loadings["N", i] < 0){ # Force C4 to have positive nitrogen loading
    PCATotal$PCA$loadings[, i] <- -1 * PCATotal$PCA$loadings[, i]
    PCATotal$PCA$scores[, i] <- -1 *PCATotal$PCA$scores[, i]
  }
}

#### UNROTATED PCA:
PCATotal$PCANoVarimax <- psych::principal(scale(PCATotal$traits), nfactors=PCATotal$dimensions,
                                          rotate="none", covar = T) # Same PCA but without rotation for comparison
sqrtEigen2 <- sqrt(colSums(PCATotal$PCANoVarimax$loadings**2))
for(i in 1:PCATotal$dimensions){
  PCATotal$PCANoVarimax$scores[, i] <- PCATotal$PCANoVarimax$scores[, i] * sqrtEigen2[i] 
}
sqrtEigenMat <- matrix(rep(sqrtEigen2, nrow(PCATotal$PCANoVarimax$loadings)), byrow=T, 
                        nrow = nrow(PCATotal$PCANoVarimax$loadings))
PCATotal$PCANoVarimax$loadings <- PCATotal$PCANoVarimax$loadings / sqrtEigenMat
# NOTE: The results in PCATotal$PCANoVarimax are identical to those using princomp(scale(PCATotal$traits))
#### Make sure orientation of all vectors is as in the paper: (same orientation rules for unrotated PCA)
for(i in 1:PCATotal$dimensions){
  if(i == 1 & PCATotal$PCANoVarimax$loadings["ph", i] < 0){ #ph is positive
    PCATotal$PCANoVarimax$loadings[, i] <- -1 * PCATotal$PCANoVarimax$loadings[, i]
    PCATotal$PCANoVarimax$scores[, i] <- -1 *PCATotal$PCANoVarimax$scores[, i]
  }
  if(i == 2 & PCATotal$PCANoVarimax$loadings["sla", i] < 0){ #sla is positive
    PCATotal$PCANoVarimax$loadings[, i] <- -1 * PCATotal$PCANoVarimax$loadings[, i]
    PCATotal$PCANoVarimax$scores[, i] <- -1 *PCATotal$PCANoVarimax$scores[, i]
  }
  if(i == 3 & PCATotal$PCANoVarimax$loadings["SRL", i] > 0){ #SRL is negative
    PCATotal$PCANoVarimax$loadings[, i] <- -1 * PCATotal$PCANoVarimax$loadings[, i]
    PCATotal$PCANoVarimax$scores[, i] <- -1 *PCATotal$PCANoVarimax$scores[, i]
  }
  if(i == 4 & PCATotal$PCANoVarimax$loadings["N", i] < 0){ #N is ppositive
    PCATotal$PCANoVarimax$loadings[, i] <- -1 * PCATotal$PCANoVarimax$loadings[, i]
    PCATotal$PCANoVarimax$scores[, i] <- -1 * PCATotal$PCANoVarimax$scores[, i]
  }
}

PCATotal$traitsUse <- data.frame(PCATotal$PCA$scores[, 1:PCATotal$dimensions]) # Species scores in 4D functional space
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse))) # Bandwidth for 4D kernel density
PCATotal$TPDs <- TPDsMean_large(species = rownames(PCATotal$traitsUse), 
                          means = PCATotal$traitsUse, 
                          sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse)), 
                                       byrow=T, ncol=PCATotal$dimensions),
                          n_divisions = gridSize) # Calculate 4D TPDs (computationally intensive)

gridSize <- 200 # Higher resolution for 2D projections
PCATotal$traitsUse2D <- list() # Container for 2D trait combinations
PCATotal$TPDs2D <- list() # Container for 2D TPD calculations
# TPD for the different planes combining pairs of dimensions
PCATotal$traitsUse2D$Comp1_Comp2 <- data.frame(PCATotal$PCA$scores[, 1:2]) # C1-C2 plane (aboveground)
PCATotal$traitsUse2D$Comp1_Comp3 <- data.frame(PCATotal$PCA$scores[, c(1, 3)]) # C1-C3 plane 
PCATotal$traitsUse2D$Comp1_Comp4 <- data.frame(PCATotal$PCA$scores[, c(1, 4)]) # C1-C4 plane
PCATotal$traitsUse2D$Comp2_Comp3 <- data.frame(PCATotal$PCA$scores[, c(2, 3)]) # C2-C3 plane
PCATotal$traitsUse2D$Comp2_Comp4 <- data.frame(PCATotal$PCA$scores[, c(2, 4)]) # C2-C4 plane
PCATotal$traitsUse2D$Comp3_Comp4 <- data.frame(PCATotal$PCA$scores[, c(3, 4)]) # C3-C4 plane (belowground) 
#Comp1. Comp2 (Calculate TPD for aboveground plane)
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp1_Comp2))) # Bandwidth for C1-C2 plane
PCATotal$TPDs2D$Comp1_Comp2 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp1_Comp2), 
                          means = PCATotal$traitsUse2D$Comp1_Comp2, 
                          sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp1_Comp2)), 
                                       byrow=T, ncol=2),
                          alpha = alphaUse, # Use 95% confidence level
                          n_divisions = gridSize)
#Comp1. Comp3 (Calculate TPD for C1-C3 plane)
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp1_Comp3)))
PCATotal$TPDs2D$Comp1_Comp3 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp1_Comp3), 
                                        means = PCATotal$traitsUse2D$Comp1_Comp3, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp1_Comp3)),
                                                     byrow=T, ncol=2),
                                        alpha = alphaUse,
                                        n_divisions = gridSize)
#Comp1. Comp4 (Calculate TPD for C1-C4 plane)
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp1_Comp4)))
PCATotal$TPDs2D$Comp1_Comp4 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp1_Comp4), 
                                        means = PCATotal$traitsUse2D$Comp1_Comp4, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp1_Comp4)),
                                                     byrow=T, ncol=2),
                                        alpha = alphaUse,
                                        n_divisions = gridSize)

#Comp2. Comp3 (Calculate TPD for C2-C3 plane)
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp2_Comp3)))
PCATotal$TPDs2D$Comp2_Comp3 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp2_Comp3), 
                                        means = PCATotal$traitsUse2D$Comp2_Comp3, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp2_Comp3)),
                                                     byrow=T, ncol=2),
                                        alpha = alphaUse,
                                        n_divisions = gridSize)
#Comp2. Comp4 (Calculate TPD for C2-C4 plane)
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp2_Comp4)))
PCATotal$TPDs2D$Comp2_Comp4 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp2_Comp4), 
                                        means = PCATotal$traitsUse2D$Comp2_Comp4, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp2_Comp4)),
                                                     byrow=T, ncol=2),
                                        alpha = alphaUse,
                                        n_divisions = gridSize)

#Comp3. Comp4 (Calculate TPD for belowground plane)
sdTraits <- sqrt(diag(Hpi.diag(PCATotal$traitsUse2D$Comp3_Comp4)))
PCATotal$TPDs2D$Comp3_Comp4 <- TPDsMean(species = rownames(PCATotal$traitsUse2D$Comp3_Comp4), 
                                        means = PCATotal$traitsUse2D$Comp3_Comp4, 
                                        sds = matrix(rep(sdTraits, nrow(PCATotal$traitsUse2D$Comp3_Comp4)),
                                                     byrow=T, ncol=2),
                                        alpha = alphaUse,
                                        n_divisions = gridSize)
PCATotal$Readme <- "This object contains the functional space created by: 1. PCA (function psych::principal) using 301 species with complete information followed by varimax rotation. 2. PCA using the same set of species but without rotation. 3 Estimation of TPDs for 4 dimensions and between pairs of rotated components." # Documentation string

saveRDS(PCATotal, paste0("data/PCATotal_CompleteObs.rds")) # Save complete analysis results for 301 species with both above/below traits











