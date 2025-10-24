
cat(paste0("\n\n Starting script #1b \n\n"))
########################################################################
### 0. Loading trait and taxonomy data at species level
########################################################################
rootTraits <- read.table("data/Root_traits.txt")[, c("SRL", "D", "RTD", "N")]
aboveTraits <- read.table("data/Above_traits.txt")
aboveTraits <- log10(aboveTraits)
aboveTaxonomy <- read.table("data/Above_taxonomy.txt")

########################################################################
### 1. Selection of all species common to both datasets
########################################################################
spCommon <- intersect(rownames(rootTraits), rownames(aboveTraits))
rootTraits <- rootTraits[spCommon, ]
aboveTraits <- aboveTraits[spCommon, ]
allTraits <- cbind(aboveTraits, rootTraits) # 1719 species

# the pairwise correlation between traits
pwc <- cor(allTraits, use = "pairwise.complete")
# spectral decomposition of the matrix to find eigenvectors and eigenvalues (variance). This is equivalent to a PCA based on correlation.
eigPCA <- eigen(pwc)
dimnames(eigPCA$vectors) <- list(rownames(pwc), paste0("PC", 1:ncol(eigPCA$vectors)))
cumsum(eigPCA$values / sum(eigPCA$values)) #cumulative variance explained by dimensions.
round(eigPCA$vectors[, 1:4], 3)
##Apply varimax rotation:
variPCA <- varimax(eigPCA$vectors[, 1:4])$loadings
