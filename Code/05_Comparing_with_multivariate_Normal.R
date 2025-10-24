### In this script we compare the distribution of species within the above- and fine-root plane and the full functional space with a null model consisting on the same number of species randomly drawn from a multivariate normal distribution.
cat(paste0("\n\n Starting script #5 \n\n"))

data <- readRDS("data/PCATotal_CompleteObs.rds")

####################################################################
####################################################################
#### NOTE: Running the following takes a very long time!!!!

RawResultsNormal <- RawResultsObserved <- list()
totalNull <- 499 ###TOTAL NUMBER OF REPETITIONS IN NULL MODELS
### For above and fine-root and all spectrum (4D), we  create 499 multivariate normal distributions with same means and variance-covariance and as many observations as species in the spectrum. Then, we  estimate functional richness, evenness and divergence and compare with the observed values.
densProfile <- seq(0.001, 0.999, by = 0.001)
RawResultsAbove <- RawResultsBelow <- RawResultsALL <- list()
RawResultsAbove[["densProfile"]] <- RawResultsBelow[["densProfile"]] <- RawResultsALL[["densProfile"]] <-
  matrix(NA, nrow = length(densProfile), ncol = totalNull + 1,
         dimnames = list(densProfile, c("Obs", paste0("Rand", 1:totalNull))))
RawResultsAbove[["FEve"]] <- RawResultsBelow[["FEve"]] <- RawResultsALL[["FEve"]] <- rep(NA, totalNull + 1)
names(RawResultsAbove[["FEve"]]) <- names(RawResultsBelow[["FEve"]]) <- names(RawResultsALL[["FEve"]]) <- c("Obs", paste0("Rand", 1:totalNull))
RawResultsAbove[["FDiv"]] <- RawResultsBelow[["FDiv"]] <- RawResultsALL[["FDiv"]] <- rep(NA, totalNull + 1)
names(RawResultsAbove[["FDiv"]]) <- names(RawResultsBelow[["FDiv"]]) <- names(RawResultsALL[["FDiv"]]) <- c("Obs", paste0("Rand", 1:totalNull))
RawResultsAbove[["FEve"]] <- RawResultsBelow[["FEve"]] <- RawResultsALL[["FEve"]] <- rep(NA, totalNull + 1)
names(RawResultsAbove[["FEve"]]) <- names(RawResultsBelow[["FEve"]]) <- names(RawResultsALL[["FEve"]]) <- c("Obs", paste0("Rand", 1:totalNull))
##### ABOVEGROUND: ####
commMatrixObs <- matrix(1, nrow = 1, ncol = nrow(data$traitsUse),
                        dimnames = list("ALLSP", names(data$TPDs2D$Comp1_Comp2$TPDs)))
TPDsObs <- data$TPDs2D$Comp1_Comp2
TPDcObs <- TPDc(TPDs = data$TPDs2D$Comp1_Comp2, sampUnit = commMatrixObs)
RENDObsAux <- REND(TPDc=TPDcObs)
RawResultsAbove[["densProfile"]][, 1]  <- densityProfileTPD(TPDcObs, probs = densProfile)
RawResultsAbove[["FEve"]][1] <- RENDObsAux$communities$FEvenness
RawResultsAbove[["FDiv"]][1] <- RENDObsAux$communities$FDivergence

traitsUSEAbove <- data$traitsUse[, 1:2]
meansAux <- colMeans(traitsUSEAbove)
covAux <- cov(traitsUSEAbove)
for(i in 1:totalNull){
  cat(paste("\r", "ABOVEGROUND Rep: ", i, "out of", totalNull, "\r"))
  normalTraits <- mvtnorm::rmvnorm(n = nrow(traitsUSEAbove), mean = meansAux, sigma=covAux)
  normalSD <- sqrt(diag(Hpi.diag(normalTraits)))
  TPDsNormal <- TPDsMean(species = names(TPDsObs$TPDs),
                           means =  normalTraits,
                           sds = matrix(rep(normalSD, nrow(normalTraits)),
                                        byrow=T, ncol=ncol(normalTraits)),
                           alpha = TPDsObs$data$alpha,
                           n_divisions = nrow(TPDsObs$data$evaluation_grid)^(1/ncol(normalTraits)))
  TPDcNormal <- TPDc(TPDs = TPDsNormal, sampUnit = commMatrixObs)
  RENDNormalAux <- REND(TPDc=TPDcNormal)
  RawResultsAbove[["densProfile"]][, 1 + i] <- densityProfileTPD(TPDcNormal, probs=densProfile)
  RawResultsAbove[["FEve"]][1 + i] <- RENDNormalAux$communities$FEvenness
  RawResultsAbove[["FDiv"]][1 + i] <- RENDNormalAux$communities$FDivergence
  TPDcNormal <- NULL
}
saveRDS(RawResultsAbove, "data/NullModelsMultivariateNormalAbove.rds")



##### fine-root: ####
commMatrixObs <- matrix(1, nrow = 1, ncol = nrow(data$traitsUse),
                        dimnames = list("ALLSP", names(data$TPDs2D$Comp3_Comp4$TPDs)))
TPDsObs <- data$TPDs2D$Comp3_Comp4
TPDcObs <- TPDc(TPDs = data$TPDs2D$Comp3_Comp4, sampUnit = commMatrixObs)
RENDObsAux <- REND(TPDc=TPDcObs)
RawResultsBelow[["densProfile"]][, 1]  <- densityProfileTPD(TPDcObs, probs = densProfile)
RawResultsBelow[["FEve"]][1] <- RENDObsAux$communities$FEvenness
RawResultsBelow[["FDiv"]][1] <- RENDObsAux$communities$FDivergence

traitsUSEBelow <- data$traitsUse[, 3:4]
meansAux <- colMeans(traitsUSEBelow)
covAux <- cov(traitsUSEBelow)
for(i in 1:totalNull){
  cat(paste("\r", "fine-root: Rep: ", i, "out of", totalNull, "\r"))
  normalTraits <- mvtnorm::rmvnorm(n = nrow(traitsUSEBelow), mean = meansAux, sigma=covAux)
  normalSD <- sqrt(diag(Hpi.diag(normalTraits)))
  TPDsNormal <- TPDsMean(species = names(TPDsObs$TPDs),
                         means =  normalTraits,
                         sds = matrix(rep(normalSD, nrow(normalTraits)),
                                      byrow=T, ncol=ncol(normalTraits)),
                         alpha = TPDsObs$data$alpha,
                         n_divisions = nrow(TPDsObs$data$evaluation_grid)^(1/ncol(normalTraits)))
  TPDcNormal <- TPDc(TPDs = TPDsNormal, sampUnit = commMatrixObs)
  RENDNormalAux <- REND(TPDc=TPDcNormal)
  RawResultsBelow[["densProfile"]][, 1 + i] <- densityProfileTPD(TPDcNormal, probs=densProfile)
  RawResultsBelow[["FEve"]][1 + i] <- RENDNormalAux$communities$FEvenness
  RawResultsBelow[["FDiv"]][1 + i] <- RENDNormalAux$communities$FDivergence
  TPDcNormal <- NULL
}
saveRDS(RawResultsBelow, "data/NullModelsMultivariateNormalBelow.rds")

##### ALL: ####
densProfile <- seq(0.001, 0.999, by = 0.001)

commMatrixObs <- matrix(1, nrow = 1, ncol = nrow(data$traitsUse),
                        dimnames = list("ALLSP", names(data$TPDs$TPDs)))
TPDsObs <- data$TPDs
TPDcObs <- TPDc_large(TPDs = data$TPDs, sampUnit = commMatrixObs)
RENDObsAux <- REND_large(TPDc=TPDcObs)
RawResultsALL[["densProfile"]][, 1]  <- densityProfileTPD_large(TPDcObs, probs = densProfile)
RawResultsALL[["FEve"]][1] <- RENDObsAux$communities$FEvenness
RawResultsALL[["FDiv"]][1] <- RENDObsAux$communities$FDivergence

traitsUSEALL <- data$traitsUse
meansAux <- colMeans(traitsUSEALL)
covAux <- cov(traitsUSEALL)
for(i in 1:totalNull){
  cat(paste("\n", "ALL: Rep: ", i, "out of", totalNull))
  normalTraits <- mvtnorm::rmvnorm(n = nrow(traitsUSEALL), mean = meansAux, sigma=covAux)
  normalSD <- sqrt(diag(Hpi.diag(normalTraits)))
  TPDsNormal <- TPDsMean_large(species = names(data$TPDs$TPDs),
                         means =  normalTraits,
                         sds = matrix(rep(normalSD, nrow(normalTraits)),
                                      byrow=T, ncol=ncol(normalTraits)),
                         alpha = TPDsObs$data$alpha,
                         n_divisions = nrow(data$TPDs$data$evaluation_grid)^(1/ncol(normalTraits)))
  TPDcNormal <- TPDc_large(TPDs = TPDsNormal, sampUnit = commMatrixObs)
  RENDNormalAux <- REND_large(TPDc=TPDcNormal)
  RawResultsALL[["densProfile"]][, 1 + i] <- densityProfileTPD_large(TPDcNormal, probs=densProfile)
  RawResultsALL[["FEve"]][1 + i] <- RENDNormalAux$communities$FEvenness
  RawResultsALL[["FDiv"]][1 + i] <- RENDNormalAux$communities$FDivergence
  TPDcNormal <- NULL
  if(i %% 10 ==0) saveRDS(RawResultsALL, "data/NullModelsMultivariateNormalTotal.rds")
}
saveRDS(RawResultsALL, "data/NullModelsMultivariateNormalTotal.rds")