### In this script we compare the positions of common species between the spaces created considering only above and only fine-root traits and the above- and fine-root planes of the spectrum considering the 10 traits. We do this by means of a procrustes test
cat(paste0("\n\n Starting script #4 \n\n"))

data <- readRDS("data/PCATotal_CompleteObs.rds")
dataAbove <- readRDS("data/PCAAboveONLY_CompleteObs.rds")
dataBelow <- readRDS("data/PCABelowONLY_CompleteObs.rds")


FullAboveS <- data$traitsUse[, 1:2]
OnlyAboveS <- dataAbove$traitsUse
FullInAbove <- which(rownames(FullAboveS) %in% rownames(OnlyAboveS))
FullAboveS <- FullAboveS[FullInAbove,]
AboveInFull <- which(rownames(OnlyAboveS) %in% rownames(FullAboveS))
OnlyAboveS <- OnlyAboveS[AboveInFull, ][rownames(FullAboveS),]
cat("\n\n --> Full aboveground (C1-C2) vs only aboveground\n\n")
print(ade4::procuste.rtest(as.data.frame(FullAboveS), as.data.frame(OnlyAboveS), nrepet = 9999))


FullBelowS <- data$traitsUse[, 3:4]
OnlyBelowS <- dataBelow$traitsUse
FullInBelow <- which(rownames(FullBelowS) %in% rownames(OnlyBelowS))
FullBelowS <- FullBelowS[FullInBelow,]
BelowInFull <- which(rownames(OnlyBelowS) %in% rownames(FullBelowS))
OnlyBelowS <- OnlyBelowS[BelowInFull, ][rownames(FullBelowS),]
cat("\n\n --> Full fine roots (C3-C4) vs only fine roots\n\n")
print(ade4::procuste.rtest(as.data.frame(FullBelowS), as.data.frame(OnlyBelowS), nrepet = 9999))


OnlyAboveS <- dataAbove$traitsUse
OnlyBelowS <- dataBelow$traitsUse
AboveInBelow <- which(rownames(OnlyAboveS) %in% rownames(OnlyBelowS))
OnlyAboveS <- OnlyAboveS[AboveInBelow,]
BelowInAbove <- which(rownames(OnlyBelowS) %in% rownames(OnlyAboveS))
OnlyBelowS <- OnlyBelowS[BelowInAbove, ][rownames(OnlyAboveS),]
cat("\n\n --> Only above vs only fine roots\n\n")
print(ade4::procuste.rtest(as.data.frame(OnlyAboveS), as.data.frame(OnlyBelowS), nrepet = 9999))
