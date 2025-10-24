### In this script we perform permanova analysis to estimate how much variation of the species in the aboveground plane, belowground plane and total functional space is explained by growth form (woodiness), family and biome: 
cat(paste0("\n\n Starting script #7 \n\n"))

data <- readRDS("data/PCATotal_ImputedObs.rds")

################################################################ 
################################################################ 
#### 1. Permanova Woodiness #### 
woodiness <- read.table("data/woodiness_1218sp.txt", header = T)

permWoodyTraits <- list()
for(i in 1:ncol(data$traits)){
  permWoodyTraits[[i]] <- 100 * adonis(dist(data$traits[, i]) ~ woodiness$woodiness,
                                 permutations = 0)$aov.tab$R2[1]
  names(permWoodyTraits)[i] <- colnames(data$traits)[i]
}
PermWoodyC1 <- 100 * adonis(dist(data$traitsUse[, 1]) ~ woodiness$woodiness, 
                            permutations = 0)$aov.tab$R2[1]
PermWoodyC2 <- 100 * adonis(dist(data$traitsUse[, 2]) ~ woodiness$woodiness, 
                            permutations = 0)$aov.tab$R2[1]
PermWoodyC3 <- 100 * adonis(dist(data$traitsUse[, 3]) ~ woodiness$woodiness, 
                            permutations = 0)$aov.tab$R2[1]
PermWoodyC4 <- 100 * adonis(dist(data$traitsUse[, 4]) ~ woodiness$woodiness, 
                            permutations = 0)$aov.tab$R2[1]

PermWoodyAbove <- 100 * adonis(dist(data$traitsUse[, 1:2]) ~ woodiness$woodiness, 
                         permutations = 0)$aov.tab$R2[1]
PermWoodyBelow <- 100 * adonis(dist(data$traitsUse[, 3:4]) ~ woodiness$woodiness, 
                         permutations = 0)$aov.tab$R2[1]
PermWoodyFull <- 100 * adonis(dist(data$traitsUse[, 1:4]) ~ woodiness$woodiness, 
                        permutations = 0)$aov.tab$R2[1]

#Permanova. Family:
permFamTraits <- list()
for(i in 1:ncol(data$traits)){
  permFamTraits[[i]] <- 100 * adonis(dist(data$traits[, i]) ~ data$AllInfo$family,
                               permutations = 0)$aov.tab$R2[1]
  names(permFamTraits)[i] <- colnames(data$traits)[i]
}

PermFamC1 <- 100 * adonis(dist(data$traitsUse[, 1]) ~ data$AllInfo$family, permutations = 0)$aov.tab$R2[1]
PermFamC2 <- 100 * adonis(dist(data$traitsUse[, 2]) ~ data$AllInfo$family, permutations = 0)$aov.tab$R2[1]
PermFamC3 <- 100 * adonis(dist(data$traitsUse[, 3]) ~ data$AllInfo$family, permutations = 0)$aov.tab$R2[1]
PermFamC4 <- 100 * adonis(dist(data$traitsUse[, 4]) ~ data$AllInfo$family, permutations = 0)$aov.tab$R2[1]

PermFamAbove <- 100 * adonis(dist(data$traitsUse[, 1:2]) ~ data$AllInfo$family, 
                       permutations = 0)$aov.tab$R2[1]
PermFamBelow <- 100 * adonis(dist(data$traitsUse[, 3:4]) ~ data$AllInfo$family, 
                       permutations = 0)$aov.tab$R2[1]
PermFamFull <- 100 * adonis(dist(data$traitsUse[, 1:4]) ~ data$AllInfo$family, 
                      permutations = 0)$aov.tab$R2[1]

#Permanova. Biome:
#There is more than 1 biome per species in some cases.
#We permanova by assigning biomes to species in proportion to the observations in each biome:
biomesMat <- read.table("data/biomes_1218sp.txt", header = T)
nreps <- 500
permBiomeAbove <- permBiomeBelow <- permBiomeFull <- 
  permBiomeC1 <- permBiomeC2 <- permBiomeC3 <- permBiomeC4 <- 
  permBiome_la <- permBiome_ln<-permBiome_ph<-permBiome_sla<-permBiome_ssd<-permBiome_sm<- 
  permBiome_SRL<-permBiome_D<-permBiome_RTD<-permBiome_N<- rep(NA, nreps)
biomesNames <- colnames(biomesMat)
for(rep in 1:nreps){
  cat(paste("\r Rep: ", rep, "\r"))
  biomesAux <- rep(NA, nrow(biomesMat))
  for(i in 1:nrow(biomesMat)){
    if(sum(biomesMat[i,]) > 0){
      biomesAux[i] <- sample(biomesNames, size = 1, prob = biomesMat[i,])  
    }  
  }
  spRemove <- which(is.na(biomesAux))
  traitsUseAux <- data$traitsUse
  if(length(spRemove) > 0){
    traitsUseAux <- data$traitsUse[- spRemove, ]
    biomesAux <- biomesAux[- spRemove]
    
  }  
  

  permBiome_ph[rep] <- 100 * adonis(dist(data$traits$ph) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  permBiome_ssd[rep] <- 100 * adonis(dist(data$traits$ssd) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  permBiome_sm[rep] <- 100 * adonis(dist(data$traits$sm) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  
  permBiome_la[rep] <-100 *  adonis(dist(data$traits$la) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  permBiome_ln[rep] <- 100 * adonis(dist(data$traits$ln) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  permBiome_sla[rep] <-100 *  adonis(dist(data$traits$sla) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  
  permBiome_SRL[rep] <-100 *  adonis(dist(data$traits$SRL) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  permBiome_D[rep] <- 100 * adonis(dist(data$traits$D) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  
  permBiome_RTD[rep] <-100 *  adonis(dist(data$traits$RTD) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  permBiome_N[rep] <- 100 * adonis(dist(data$traits$N) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  
  permBiomeC1[rep] <-100 *  adonis(dist(data$traitsUse[, 1]) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  permBiomeC2[rep] <-100 *  adonis(dist(data$traitsUse[, 2]) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  permBiomeC3[rep] <-100 *  adonis(dist(data$traitsUse[, 3]) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  permBiomeC4[rep] <- 100 * adonis(dist(data$traitsUse[, 4]) ~ biomesAux, permutations = 0)$aov.tab$R2[1]
  
  permBiomeAbove[rep] <- 100 * adonis(dist(traitsUseAux[, 1:2]) ~ biomesAux, 
                                permutations = 0)$aov.tab$R2[1]
  permBiomeBelow[rep] <- 100 * adonis(dist(traitsUseAux[, 3:4]) ~ biomesAux, 
                                permutations = 0)$aov.tab$R2[1]
  permBiomeFull[rep] <- 100 * adonis(dist(traitsUseAux[, 1:4]) ~ biomesAux, 
                               permutations = 0)$aov.tab$R2[1]
}

ResultsPermanova <- list()
ResultsPermanova$woodiness <- list(traits = permWoodyTraits,
                                   C1 = PermWoodyC1, C2 =PermWoodyC2, C3 =PermWoodyC3, 
                                   C4 =PermWoodyC4, 
                                   Above = PermWoodyAbove, Below = PermWoodyBelow, 
                                   Total = PermWoodyFull)
ResultsPermanova$families <- list(traits = permFamTraits,
                                  components = list(C1 = PermFamC1, C2 =PermFamC2, 
                                                    C3 = PermFamC3, C4 = PermFamC4), 
                                  Above = PermFamAbove, Below = PermFamBelow,
                                  Total = PermFamFull)
ResultsPermanova$biomes <- list(traits = list(ph= permBiome_ph, ssd = permBiome_ssd, 
                                              sm = permBiome_sm, la = permBiome_la, 
                                              ln = permBiome_ln, sla = permBiome_sla, 
                                              SRL = permBiome_SRL, D= permBiome_D, 
                                              RTD = permBiome_RTD, N = permBiome_N),
                                components = list(C1 = permBiomeC1, C2 = permBiomeC2, 
                                                  C3 = permBiomeC3, C4= permBiomeC4), 
                                Above = permBiomeAbove, Below = permBiomeBelow, 
                                Total = permBiomeFull)


saveRDS(ResultsPermanova, "data/ResultsPermanova.rds")