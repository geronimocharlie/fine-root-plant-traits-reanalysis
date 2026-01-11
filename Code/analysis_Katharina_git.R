# Import packages
library(dplyr)
library(psych)
library(ks)
library(TPD)
library(paran)


# Import data
rootTraits <- read.table("../data/Root_traits.txt")[, c("SRL", "D", "RTD", "N", "biomesKoeppenGroup")] # already log transformed
aboveTraits <- read.table("../data/Above_traits.txt")
aboveTraits <- log10(aboveTraits) # log transform 
aboveTaxonomy <- read.table("../data/Above_taxonomy.txt")
rootBiomes <- read.table("../data/Root_traits.txt")[, c("biomesKoeppenGroup")]

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
imputed_data <- readRDS("../data/imputedTraits.rds")
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
# 1. PCA imputed combined full
gridSize <- 200 # had to reduce due to RAM capacities, originally 200
PCA_combined_full <- list()
PCA_combined_full$traits <- imputed_data_clean
PCA_combined_full$dimensions <- paran(PCA_combined_full$traits)$Retained
PCA_combined_full$means <- apply(PCA_combined_full$traits, 2, mean)
PCA_combined_full$sds <- apply(PCA_combined_full$traits, 2, sd)
PCA_combined_full$PCA <- psych::principal(scale(PCA_combined_full$traits), nfactors=PCA_combined_full$dimensions, 
                                 rotate="varimax", covar = TRUE) # covar = TRUE but FALSE to save R RAM
PCA_combined_full$Variance <- PCA_combined_full$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_combined_full$PCA$loadings**2))
for(i in 1:PCA_combined_full$dimensions){
  PCA_combined_full$PCA$scores[, i] <- PCA_combined_full$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_combined_full$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_combined_full$PCA$loadings))
PCA_combined_full$PCA$loadings <- PCA_combined_full$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_combined_full$traitsUse <- data.frame(PCA_combined_full$PCA$scores[, 1:PCA_combined_full$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_combined_full$traitsUse)))
PCA_combined_full$TPDs <- TPDsMean(species = rownames(PCA_combined_full$traitsUse), 
                          means = PCA_combined_full$traitsUse, 
                          sds = matrix(rep(sdTraits, nrow(PCA_combined_full$traitsUse)), byrow=T, 
                                       ncol=PCA_combined_full$dimensions),
                          n_divisions = gridSize)
saveRDS(PCA_combined_full, paste0("../PCA_results/PCA_imputed_combined_full.rds"))


# 2. PCA imputed combined half
gridSize <- 200
PCA_combined_half <- list()
PCA_combined_half$traits <- imputed_data_half
PCA_combined_half$dimensions <- paran(PCA_combined_half$traits)$Retained
PCA_combined_half$means <- apply(PCA_combined_half$traits, 2, mean)
PCA_combined_half$sds <- apply(PCA_combined_half$traits, 2, sd)
PCA_combined_half$PCA <- psych::principal(scale(PCA_combined_half$traits), nfactors=PCA_combined_half$dimensions, 
                                          rotate="varimax", covar = T)
PCA_combined_half$Variance <- PCA_combined_half$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_combined_half$PCA$loadings**2))
for(i in 1:PCA_combined_half$dimensions){
  PCA_combined_half$PCA$scores[, i] <- PCA_combined_half$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_combined_half$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_combined_half$PCA$loadings))
PCA_combined_half$PCA$loadings <- PCA_combined_half$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_combined_half$traitsUse <- data.frame(PCA_combined_half$PCA$scores[, 1:PCA_combined_half$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_combined_half$traitsUse)))
PCA_combined_half$TPDs <- TPDsMean(species = rownames(PCA_combined_half$traitsUse), 
                                   means = PCA_combined_half$traitsUse, 
                                   sds = matrix(rep(sdTraits, nrow(PCA_combined_half$traitsUse)), byrow=T, 
                                                ncol=PCA_combined_half$dimensions),
                                   n_divisions = gridSize)
saveRDS(PCA_combined_half, paste0("../PCA_results/PCA_imputed_combined_half.rds"))


# 3. PCA imputed combined quarter
gridSize <- 200
PCA_combined_quarter <- list()
PCA_combined_quarter$traits <- imputed_data_quarter
PCA_combined_quarter$dimensions <- paran(PCA_combined_quarter$traits)$Retained
PCA_combined_quarter$means <- apply(PCA_combined_quarter$traits, 2, mean)
PCA_combined_quarter$sds <- apply(PCA_combined_quarter$traits, 2, sd)
PCA_combined_quarter$PCA <- psych::principal(scale(PCA_combined_quarter$traits), nfactors=PCA_combined_quarter$dimensions, 
                                          rotate="varimax", covar = T)
PCA_combined_quarter$Variance <- PCA_combined_quarter$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_combined_quarter$PCA$loadings**2))
for(i in 1:PCA_combined_quarter$dimensions){
  PCA_combined_quarter$PCA$scores[, i] <- PCA_combined_quarter$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_combined_quarter$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_combined_quarter$PCA$loadings))
PCA_combined_quarter$PCA$loadings <- PCA_combined_quarter$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_combined_quarter$traitsUse <- data.frame(PCA_combined_quarter$PCA$scores[, 1:PCA_combined_quarter$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_combined_quarter$traitsUse)))
PCA_combined_quarter$TPDs <- TPDsMean(species = rownames(PCA_combined_quarter$traitsUse), 
                                   means = PCA_combined_quarter$traitsUse, 
                                   sds = matrix(rep(sdTraits, nrow(PCA_combined_quarter$traitsUse)), byrow=T, 
                                                ncol=PCA_combined_quarter$dimensions),
                                   n_divisions = gridSize)
saveRDS(PCA_combined_quarter, paste0("../PCA_results/PCA_imputed_combined_quarter.rds"))


# 4. PCA imputed combined tenth
gridSize <- 200
PCA_combined_tenth <- list()
PCA_combined_tenth$traits <- imputed_data_tenth
PCA_combined_tenth$dimensions <- paran(PCA_combined_tenth$traits)$Retained
PCA_combined_tenth$means <- apply(PCA_combined_tenth$traits, 2, mean)
PCA_combined_tenth$sds <- apply(PCA_combined_tenth$traits, 2, sd)
PCA_combined_tenth$PCA <- psych::principal(scale(PCA_combined_tenth$traits), nfactors=PCA_combined_tenth$dimensions, 
                                             rotate="varimax", covar = T)
PCA_combined_tenth$Variance <- PCA_combined_tenth$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_combined_tenth$PCA$loadings**2))
for(i in 1:PCA_combined_tenth$dimensions){
  PCA_combined_tenth$PCA$scores[, i] <- PCA_combined_tenth$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_combined_tenth$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_combined_tenth$PCA$loadings))
PCA_combined_tenth$PCA$loadings <- PCA_combined_tenth$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_combined_tenth$traitsUse <- data.frame(PCA_combined_tenth$PCA$scores[, 1:PCA_combined_tenth$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_combined_tenth$traitsUse)))
PCA_combined_tenth$TPDs <- TPDsMean(species = rownames(PCA_combined_tenth$traitsUse), 
                                      means = PCA_combined_tenth$traitsUse, 
                                      sds = matrix(rep(sdTraits, nrow(PCA_combined_tenth$traitsUse)), byrow=T, 
                                                   ncol=PCA_combined_tenth$dimensions),
                                      n_divisions = gridSize)
saveRDS(PCA_combined_tenth, paste0("../PCA_results/PCA_imputed_combined_tenth.rds"))

#  now compare if PCA is stable, if not, no need to look at biome specific PCAs


# 2. PCA below
imputed_below <- subset(imputed_data_clean, select = -c(la, ln, ph, sla, ssd, sm))

gridSize <- 200
PCA_imputed_below <- list()
PCA_imputed_below$traits <- imputed_below
PCA_imputed_below$dimensions <- paran(PCA_imputed_below$traits)$Retained
PCA_imputed_below$means <- apply(PCA_imputed_below$traits, 2, mean)
PCA_imputed_below$sds <- apply(PCA_imputed_below$traits, 2, sd)
PCA_imputed_below$PCA <- psych::principal(scale(PCA_imputed_below$traits), nfactors=PCA_imputed_below$dimensions, 
                                           rotate="varimax", covar = T)
PCA_imputed_below$Variance <- PCA_imputed_below$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_imputed_below$PCA$loadings**2))
for(i in 1:PCA_imputed_below$dimensions){
  PCA_imputed_below$PCA$scores[, i] <- PCA_imputed_below$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_imputed_below$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_imputed_below$PCA$loadings))
PCA_imputed_below$PCA$loadings <- PCA_imputed_below$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_imputed_below$traitsUse <- data.frame(PCA_imputed_below$PCA$scores[, 1:PCA_imputed_below$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_imputed_below$traitsUse)))
PCA_imputed_below$TPDs <- TPDsMean(species = rownames(PCA_imputed_below$traitsUse), 
                                    means = PCA_imputed_below$traitsUse, 
                                    sds = matrix(rep(sdTraits, nrow(PCA_imputed_below$traitsUse)), byrow=T, 
                                                 ncol=PCA_imputed_below$dimensions),
                                    n_divisions = gridSize)
saveRDS(PCA_imputed_below, paste0("../PCA_results/PCA_imputed_below.rds"))


# 3. PCA above
imputed_above <- subset(imputed_data_clean, select = -c(SRL, D, RTD, N))# imputed data here less than not imputed data 

gridSize <- 200
PCA_imputed_above <- list()
PCA_imputed_above$traits <- imputed_above
PCA_imputed_above$dimensions <- paran(PCA_imputed_above$traits)$Retained
PCA_imputed_above$means <- apply(PCA_imputed_above$traits, 2, mean)
PCA_imputed_above$sds <- apply(PCA_imputed_above$traits, 2, sd)
PCA_imputed_above$PCA <- psych::principal(scale(PCA_imputed_above$traits), nfactors=PCA_imputed_above$dimensions, 
                                          rotate="varimax", covar = T)
PCA_imputed_above$Variance <- PCA_imputed_above$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_imputed_above$PCA$loadings**2))
for(i in 1:PCA_imputed_above$dimensions){
  PCA_imputed_above$PCA$scores[, i] <- PCA_imputed_above$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_imputed_above$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_imputed_above$PCA$loadings))
PCA_imputed_above$PCA$loadings <- PCA_imputed_above$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_imputed_above$traitsUse <- data.frame(PCA_imputed_above$PCA$scores[, 1:PCA_imputed_above$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_imputed_above$traitsUse)))
PCA_imputed_above$TPDs <- TPDsMean(species = rownames(PCA_imputed_above$traitsUse), 
                                   means = PCA_imputed_above$traitsUse, 
                                   sds = matrix(rep(sdTraits, nrow(PCA_imputed_above$traitsUse)), byrow=T, 
                                                ncol=PCA_imputed_above$dimensions),
                                   n_divisions = gridSize)
saveRDS(PCA_imputed_above, paste0("../PCA_results/PCA_imputed_above.rds"))


# 4. PCA für imputed biomes (seperated), falls PCA stabil
gridSize <- 200
PCA_imputed_continental <- list()
PCA_imputed_continental$traits <- imputed_subset_continental_nobiome
PCA_imputed_continental$dimensions <- paran(PCA_imputed_continental$traits)$Retained
PCA_imputed_continental$means <- apply(PCA_imputed_continental$traits, 2, mean)
PCA_imputed_continental$sds <- apply(PCA_imputed_continental$traits, 2, sd)
PCA_imputed_continental$PCA <- psych::principal(scale(PCA_imputed_continental$traits), nfactors=PCA_imputed_continental$dimensions, 
                                          rotate="varimax", covar = T)
PCA_imputed_continental$Variance <- PCA_imputed_continental$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_imputed_continental$PCA$loadings**2))
for(i in 1:PCA_imputed_continental$dimensions){
  PCA_imputed_continental$PCA$scores[, i] <- PCA_imputed_continental$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_imputed_continental$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_imputed_continental$PCA$loadings))
PCA_imputed_continental$PCA$loadings <- PCA_imputed_continental$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_imputed_continental$traitsUse <- data.frame(PCA_imputed_continental$PCA$scores[, 1:PCA_imputed_continental$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_imputed_continental$traitsUse)))
PCA_imputed_continental$TPDs <- TPDsMean(species = rownames(PCA_imputed_continental$traitsUse), 
                                   means = PCA_imputed_continental$traitsUse, 
                                   sds = matrix(rep(sdTraits, nrow(PCA_imputed_continental$traitsUse)), byrow=T, 
                                                ncol=PCA_imputed_continental$dimensions),
                                   n_divisions = gridSize)
saveRDS(PCA_imputed_continental, paste0("../PCA_results/PCA_imputed_continental.rds"))


gridSize <- 200
PCA_imputed_temperate <- list()
PCA_imputed_temperate$traits <- imputed_subset_temperate_nobiome
PCA_imputed_temperate$dimensions <- paran(PCA_imputed_temperate$traits)$Retained
PCA_imputed_temperate$means <- apply(PCA_imputed_temperate$traits, 2, mean)
PCA_imputed_temperate$sds <- apply(PCA_imputed_temperate$traits, 2, sd)
PCA_imputed_temperate$PCA <- psych::principal(scale(PCA_imputed_temperate$traits), nfactors=PCA_imputed_temperate$dimensions, 
                                                rotate="varimax", covar = T)
PCA_imputed_temperate$Variance <- PCA_imputed_temperate$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_imputed_temperate$PCA$loadings**2))
for(i in 1:PCA_imputed_temperate$dimensions){
  PCA_imputed_temperate$PCA$scores[, i] <- PCA_imputed_temperate$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_imputed_temperate$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_imputed_temperate$PCA$loadings))
PCA_imputed_temperate$PCA$loadings <- PCA_imputed_temperate$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_imputed_temperate$traitsUse <- data.frame(PCA_imputed_temperate$PCA$scores[, 1:PCA_imputed_temperate$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_imputed_temperate$traitsUse)))
PCA_imputed_temperate$TPDs <- TPDsMean(species = rownames(PCA_imputed_temperate$traitsUse), 
                                         means = PCA_imputed_temperate$traitsUse, 
                                         sds = matrix(rep(sdTraits, nrow(PCA_imputed_temperate$traitsUse)), byrow=T, 
                                                      ncol=PCA_imputed_temperate$dimensions),
                                         n_divisions = gridSize)
saveRDS(PCA_imputed_temperate, paste0("../PCA_results/PCA_imputed_temperate.rds"))


gridSize <- 200
PCA_imputed_tropical <- list()
PCA_imputed_tropical$traits <- imputed_subset_tropical_nobiome
PCA_imputed_tropical$dimensions <- paran(PCA_imputed_tropical$traits)$Retained
PCA_imputed_tropical$means <- apply(PCA_imputed_tropical$traits, 2, mean)
PCA_imputed_tropical$sds <- apply(PCA_imputed_tropical$traits, 2, sd)
PCA_imputed_tropical$PCA <- psych::principal(scale(PCA_imputed_tropical$traits), nfactors=PCA_imputed_tropical$dimensions, 
                                              rotate="varimax", covar = T)
PCA_imputed_tropical$Variance <- PCA_imputed_tropical$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_imputed_tropical$PCA$loadings**2))
for(i in 1:PCA_imputed_tropical$dimensions){
  PCA_imputed_tropical$PCA$scores[, i] <- PCA_imputed_tropical$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_imputed_tropical$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_imputed_tropical$PCA$loadings))
PCA_imputed_tropical$PCA$loadings <- PCA_imputed_tropical$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_imputed_tropical$traitsUse <- data.frame(PCA_imputed_tropical$PCA$scores[, 1:PCA_imputed_tropical$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_imputed_tropical$traitsUse)))
PCA_imputed_tropical$TPDs <- TPDsMean(species = rownames(PCA_imputed_tropical$traitsUse), 
                                       means = PCA_imputed_tropical$traitsUse, 
                                       sds = matrix(rep(sdTraits, nrow(PCA_imputed_tropical$traitsUse)), byrow=T, 
                                                    ncol=PCA_imputed_tropical$dimensions),
                                       n_divisions = gridSize)
saveRDS(PCA_imputed_tropical, paste0("../PCA_results/PCA_imputed_tropical.rds"))


# ==== NOT IMPUTED ====
# 1. PCA combined
gridSize <- 200
PCA_combined <- list()
PCA_combined$traits <- AllTraits_nobiomes
PCA_combined$dimensions <- paran(PCA_combined$traits)$Retained
PCA_combined$means <- apply(PCA_combined$traits, 2, mean)
PCA_combined$sds <- apply(PCA_combined$traits, 2, sd)
PCA_combined$PCA <- psych::principal(scale(PCA_combined$traits), nfactors=PCA_combined$dimensions, 
                                             rotate="varimax", covar = T)
PCA_combined$Variance <- PCA_combined$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_combined$PCA$loadings**2))
for(i in 1:PCA_combined$dimensions){
  PCA_combined$PCA$scores[, i] <- PCA_combined$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_combined$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_combined$PCA$loadings))
PCA_combined$PCA$loadings <- PCA_combined$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_combined$traitsUse <- data.frame(PCA_combined$PCA$scores[, 1:PCA_combined$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_combined$traitsUse)))
PCA_combined$TPDs <- TPDsMean(species = rownames(PCA_combined$traitsUse), 
                                      means = PCA_combined$traitsUse, 
                                      sds = matrix(rep(sdTraits, nrow(PCA_combined$traitsUse)), byrow=T, 
                                                   ncol=PCA_combined$dimensions),
                                      n_divisions = gridSize)
saveRDS(PCA_combined, paste0("../PCA_results/PCA_combined.rds"))


# 2. PCA below
gridSize <- 200
PCA_below <- list()
PCA_below$traits <- rootTraitsComplete_nobiomes
PCA_below$dimensions <- paran(PCA_below$traits)$Retained
PCA_below$means <- apply(PCA_below$traits, 2, mean)
PCA_below$sds <- apply(PCA_below$traits, 2, sd)
PCA_below$PCA <- psych::principal(scale(PCA_below$traits), nfactors=PCA_below$dimensions, 
                                     rotate="varimax", covar = T)
PCA_below$Variance <- PCA_below$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_below$PCA$loadings**2))
for(i in 1:PCA_below$dimensions){
  PCA_below$PCA$scores[, i] <- PCA_below$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_below$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_below$PCA$loadings))
PCA_below$PCA$loadings <- PCA_below$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_below$traitsUse <- data.frame(PCA_below$PCA$scores[, 1:PCA_below$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_below$traitsUse)))
PCA_below$TPDs <- TPDsMean(species = rownames(PCA_below$traitsUse), 
                              means = PCA_below$traitsUse, 
                              sds = matrix(rep(sdTraits, nrow(PCA_below$traitsUse)), byrow=T, 
                                           ncol=PCA_below$dimensions),
                              n_divisions = gridSize)
saveRDS(PCA_below, paste0("../PCA_results/PCA_below.rds"))


# 3. PCA above
gridSize <- 200
PCA_above <- list()
PCA_above$traits <- aboveTraitscomplete # more data here than in imputed data set
PCA_above$dimensions <- paran(PCA_above$traits)$Retained
PCA_above$means <- apply(PCA_above$traits, 2, mean)
PCA_above$sds <- apply(PCA_above$traits, 2, sd)
PCA_above$PCA <- psych::principal(scale(PCA_above$traits), nfactors=PCA_above$dimensions, 
                                  rotate="varimax", covar = T)
PCA_above$Variance <- PCA_above$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_above$PCA$loadings**2))
for(i in 1:PCA_above$dimensions){
  PCA_above$PCA$scores[, i] <- PCA_above$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_above$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_above$PCA$loadings))
PCA_above$PCA$loadings <- PCA_above$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_above$traitsUse <- data.frame(PCA_above$PCA$scores[, 1:PCA_above$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_above$traitsUse)))
PCA_above$TPDs <- TPDsMean(species = rownames(PCA_above$traitsUse), 
                           means = PCA_above$traitsUse, 
                           sds = matrix(rep(sdTraits, nrow(PCA_above$traitsUse)), byrow=T, 
                                        ncol=PCA_above$dimensions),
                           n_divisions = gridSize)
saveRDS(PCA_above, paste0("../PCA_results/PCA_above.rds"))


# 4. PCA per biome
gridSize <- 200
PCA_temperate <- list()
PCA_temperate$traits <- subset_temperate_nobiome
PCA_temperate$dimensions <- paran(PCA_temperate$traits)$Retained
PCA_temperate$means <- apply(PCA_temperate$traits, 2, mean)
PCA_temperate$sds <- apply(PCA_temperate$traits, 2, sd)
PCA_temperate$PCA <- psych::principal(scale(PCA_temperate$traits), nfactors=PCA_temperate$dimensions, 
                                  rotate="varimax", covar = T)
PCA_temperate$Variance <- PCA_temperate$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_temperate$PCA$loadings**2))
for(i in 1:PCA_temperate$dimensions){
  PCA_temperate$PCA$scores[, i] <- PCA_temperate$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_temperate$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_temperate$PCA$loadings))
PCA_temperate$PCA$loadings <- PCA_temperate$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_temperate$traitsUse <- data.frame(PCA_temperate$PCA$scores[, 1:PCA_temperate$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_temperate$traitsUse)))
PCA_temperate$TPDs <- TPDsMean(species = rownames(PCA_temperate$traitsUse), 
                           means = PCA_temperate$traitsUse, 
                           sds = matrix(rep(sdTraits, nrow(PCA_temperate$traitsUse)), byrow=T, 
                                        ncol=PCA_temperate$dimensions),
                           n_divisions = gridSize)
saveRDS(PCA_temperate, paste0("../PCA_results/PCA_temperate.rds"))


gridSize <- 200
PCA_tropical <- list()
PCA_tropical$traits <- subset_tropical_nobiome
PCA_tropical$dimensions <- paran(PCA_tropical$traits)$Retained
PCA_tropical$means <- apply(PCA_tropical$traits, 2, mean)
PCA_tropical$sds <- apply(PCA_tropical$traits, 2, sd)
PCA_tropical$PCA <- psych::principal(scale(PCA_tropical$traits), nfactors=PCA_tropical$dimensions, 
                                  rotate="varimax", covar = T)
PCA_tropical$Variance <- PCA_tropical$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_tropical$PCA$loadings**2))
for(i in 1:PCA_tropical$dimensions){
  PCA_tropical$PCA$scores[, i] <- PCA_tropical$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_tropical$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_tropical$PCA$loadings))
PCA_tropical$PCA$loadings <- PCA_tropical$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_tropical$traitsUse <- data.frame(PCA_tropical$PCA$scores[, 1:PCA_tropical$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_tropical$traitsUse)))
PCA_tropical$TPDs <- TPDsMean(species = rownames(PCA_tropical$traitsUse), 
                           means = PCA_tropical$traitsUse, 
                           sds = matrix(rep(sdTraits, nrow(PCA_tropical$traitsUse)), byrow=T, 
                                        ncol=PCA_tropical$dimensions),
                           n_divisions = gridSize)
saveRDS(PCA_tropical, paste0("../PCA_results/PCA_tropical.rds"))


gridSize <- 200
PCA_continental <- list()
PCA_continental$traits <- subset_continental_nobiome
PCA_continental$dimensions <- paran(PCA_continental$traits)$Retained
PCA_continental$means <- apply(PCA_continental$traits, 2, mean)
PCA_continental$sds <- apply(PCA_continental$traits, 2, sd)
PCA_continental$PCA <- psych::principal(scale(PCA_continental$traits), nfactors=PCA_continental$dimensions, 
                                  rotate="varimax", covar = T)
PCA_continental$Variance <- PCA_continental$PCA$Vaccounted[2,]
sqrtEigen <- sqrt(colSums(PCA_continental$PCA$loadings**2))
for(i in 1:PCA_continental$dimensions){
  PCA_continental$PCA$scores[, i] <- PCA_continental$PCA$scores[, i] * sqrtEigen[i] 
}
# Loadings from psych::principal are expressed as eigenvectors * sqrt(eigenvalues). 
# Let's express them as eigenvectors (without scale):
sqrtEigenMat <- matrix(rep(sqrtEigen, nrow(PCA_continental$PCA$loadings)), byrow=T, 
                       nrow = nrow(PCA_continental$PCA$loadings))
PCA_continental$PCA$loadings <- PCA_continental$PCA$loadings / sqrtEigenMat
# Check all are 1: colSums(PCAAbove$PCA$loadings**2)
PCA_continental$traitsUse <- data.frame(PCA_continental$PCA$scores[, 1:PCA_continental$dimensions])
sdTraits <- sqrt(diag(Hpi.diag(PCA_continental$traitsUse)))
PCA_continental$TPDs <- TPDsMean(species = rownames(PCA_continental$traitsUse), 
                           means = PCA_continental$traitsUse, 
                           sds = matrix(rep(sdTraits, nrow(PCA_continental$traitsUse)), byrow=T, 
                                        ncol=PCA_continental$dimensions),
                           n_divisions = gridSize)
saveRDS(PCA_continental, paste0("../PCA_results/PCA_continental.rds"))




# TRIAL PRETTIER PLOT (CHATGPT)
# Extract scores and loadings
#scores <- pca_combined$x
#loadings <- pca_combined$rotation

# Plot scores as points
#plot(scores[,1], scores[,2],
#     pch = 19, col = "grey40",
#     xlab = "PC1", ylab = "PC2",
#     main = "Combined PCA")

# Add loading arrows
#arrows(0, 0, 
#       x1 = loadings[,1] * 3, 
#       y1 = loadings[,2] * 3, 
#       length = 0.1, col = "red")

#text(loadings[,1] * 3, loadings[,2] * 3, 
#     labels = rownames(loadings),
#     col = "red", cex = 0.8, pos = 3)


# Combined by biome
#pca_temperate <- prcomp(subset_temperate_nobiome, scale. = TRUE, center = TRUE, retx = TRUE) 
#summary(pca_temperate)
#pca_temperate$rotation[1:10, 1:4]
#biplot(pca_temperate, main = "Combined", scale = 0)

#pca_tropical <- prcomp(subset_tropical_nobiome, scale. = TRUE, center = TRUE, retx = TRUE) 
#summary(pca_tropical)
#pca_tropical$rotation[1:10, 1:4]
#biplot(pca_tropical, main = "Combined", scale = 0)

#pca_continental <- prcomp(subset_continental_nobiome, scale. = TRUE, center = TRUE, retx = TRUE) 
#summary(pca_continental)
#pca_continental$rotation[1:10, 1:4]
#biplot(pca_continental, main = "Combined", scale = 0)