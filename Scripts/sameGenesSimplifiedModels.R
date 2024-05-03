#######################################################
### Code for the comparative analysis of WT vs PlaB ###
#######################################################
source("/path/to/allInternalFunctions.R")

### Expression data
expressionLevels <- readRDS("/path/to/Results/Untreated/fullModel/expressionLevels.rds")

### Modeling
## Initial rates
INSPEcT_nascent <- read.table("/path/to/Files/regulated_genes_features.xls",sep="\t",header=TRUE)

synthesis=signif(median(INSPEcT_nascent[,"synthesis_0"]),2)
processing=signif(median(INSPEcT_nascent[,"processing_0"]),2)
degradation=signif(median(INSPEcT_nascent[,"degradation_0"]),2)

exampleRates <- c(k1=synthesis
                 ,k2=processing
                 ,k3=processing
                 ,k4=0.1*processing
                 ,k5=0.25*processing
                 ,k6=1.5*processing
                 ,k7=0.5*degradation
                 ,k8=degradation
                 ,k10=degradation)

initialRates <- list(10**ceiling(log10(exampleRates))*(1+seq_along(exampleRates)*1e-5)
                    ,10**round(log10(exampleRates))*(1+seq_along(exampleRates)*1e-5)
                    ,10**floor(log10(exampleRates))*(1+seq_along(exampleRates)*1e-5))

## Simplified expression levels
expressionLevels_yesChpNpP <- unname(lapply(expressionLevels$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP,"[[","mean"))
expressionLevels_yesChpNpP <- lapply(expressionLevels_yesChpNpP,function(i){i[i==1e-10] <- 0;i})
expressionLevels_yesChpNpNoP <- lapply(expressionLevels_yesChpNpP,function(i)i[,!(colnames(i)%in%c("pyp0.33","pyn0.33","pyn0","pyp0","pyt"))])

expressionLevels_yesChpPNoNp <- expressionLevels_yesChpNpP
expressionLevels_yesChpNoNpP <- expressionLevels_yesChpNpNoP

expressionLevels_yesChpPNoNp <- lapply(expressionLevels_yesChpPNoNp,function(i){i[,c("npp0.33","npn0.33","npn0","npp0","npt")] <- (i[,c("npp0.33","npn0.33","npn0","npp0","npt")] + i[,c("nmp0.33","nmn0.33","nmn0","nmp0","nmt")]);i})
expressionLevels_yesChpPNoNp <- lapply(expressionLevels_yesChpPNoNp,function(i)i[,!(colnames(i)%in%c("nmp0.33","nmn0.33","nmn0","nmp0","nmt"))])
expressionLevels_yesChpPNoNp <- lapply(expressionLevels_yesChpPNoNp,function(i){colnames(i) <- gsub("^np","n",colnames(i));i})

expressionLevels_yesChpNoNpP <- lapply(expressionLevels_yesChpNoNpP,function(i){i[,c("npp0.33","npn0.33","npn0","npp0","npt")] <- (i[,c("npp0.33","npn0.33","npn0","npp0","npt")] + i[,c("nmp0.33","nmn0.33","nmn0","nmp0","nmt")]);i})
expressionLevels_yesChpNoNpP <- lapply(expressionLevels_yesChpNoNpP,function(i)i[,!(colnames(i)%in%c("nmp0.33","nmn0.33","nmn0","nmp0","nmt"))])
expressionLevels_yesChpNoNpP <- lapply(expressionLevels_yesChpNoNpP,function(i){colnames(i) <- gsub("^np","n",colnames(i));i})

expressionLevels_yesPNoChp <- expressionLevels_yesChpPNoNp
expressionLevels_noChpP <- expressionLevels_yesChpNoNpP

expressionLevels_yesPNoChp <- lapply(expressionLevels_yesPNoChp,function(i){i[,c("chpp0.33","chpn0.33","chpn0","chpp0","chpt")] <- (i[,c("chpp0.33","chpn0.33","chpn0","chpp0","chpt")] + i[,c("chmp0.33","chmn0.33","chmn0","chmp0","chmt")]);i})
expressionLevels_yesPNoChp <- lapply(expressionLevels_yesPNoChp,function(i)i[,!(colnames(i)%in%c("chmp0.33","chmn0.33","chmn0","chmp0","chmt"))])
expressionLevels_yesPNoChp <- lapply(expressionLevels_yesPNoChp,function(i){colnames(i) <- gsub("^chp","ch",colnames(i));i})

expressionLevels_noChpP <- lapply(expressionLevels_noChpP,function(i){i[,c("chpp0.33","chpn0.33","chpn0","chpp0","chpt")] <- (i[,c("chpp0.33","chpn0.33","chpn0","chpp0","chpt")] + i[,c("chmp0.33","chmn0.33","chmn0","chmp0","chmt")]);i})
expressionLevels_noChpP <- lapply(expressionLevels_noChpP,function(i)i[,!(colnames(i)%in%c("chmp0.33","chmn0.33","chmn0","chmp0","chmt"))])
expressionLevels_noChpP <- lapply(expressionLevels_noChpP,function(i){colnames(i) <- gsub("^chp","ch",colnames(i));i})

expressionLevels_yesChpNpNoP <- lapply(expressionLevels_yesChpNpNoP,function(i){i[i==0] <- 1e-10;i})
expressionLevels_yesChpPNoNp <- lapply(expressionLevels_yesChpPNoNp,function(i){i[i==0] <- 1e-10;i})
expressionLevels_yesChpNoNpP <- lapply(expressionLevels_yesChpNoNpP,function(i){i[i==0] <- 1e-10;i})
expressionLevels_yesPNoChp <- lapply(expressionLevels_yesPNoChp,function(i){i[i==0] <- 1e-10;i})
expressionLevels_noChpP <- lapply(expressionLevels_noChpP,function(i){i[i==0] <- 1e-10;i})

inferedRatesUntreatedMerged_yesChpNpNoP_multi <- inferRates(expressionData=expressionLevels_yesChpNpNoP # Expression data of the genes to be modeled.
														   ,expressionDataDev=NULL # Standard deviations of the genes to be modeled.
														   ,simulatedDataset=NULL # Simulated dataset if this is the case (just to produce real rates correlations).
														   ,initialRates=initialRates # List of initial rates for optimization.
														   ,TauFractions=c(0,0.33) # Time points with cellular fractionation.
														   ,TauPoly=c(0,0.33) # Time points with polysomal profiling.
														   ,TauTotal=NULL # Time points with total RNA profiling.
														   ,cpus=24 # Number of cpus.
														   ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
														   ,lowB=1e-6 # Lower boundary for the rates.
														   ,upB=1e10 # Upper boundary for the rates.
														   ,FlagDev="FC" # Cost function.
														   ,lambda=0.05 # Regularization strength.
														   ,excludeSpecies=NULL # List of species to be excluded from the cost function.
														   ,parFixed=NULL) # List of parameters to be excluded from the optimization.
saveRDS(inferedRatesUntreatedMerged_yesChpNpNoP_multi,"inferedRatesUntreatedMerged_yesChpNpNoP_multi.rds")

inferedRatesUntreatedMerged_yesChpPNoNp_multi <- inferRates(expressionData=expressionLevels_yesChpPNoNp # Expression data of the genes to be modeled.
														   ,expressionDataDev=NULL # Standard deviations of the genes to be modeled.
														   ,simulatedDataset=NULL # Simulated dataset if this is the case (just to produce real rates correlations).
														   ,initialRates=initialRates # List of initial rates for optimization.
														   ,TauFractions=c(0,0.33) # Time points with cellular fractionation.
														   ,TauPoly=c(0,0.33) # Time points with polysomal profiling.
														   ,TauTotal=NULL # Time points with total RNA profiling.
														   ,cpus=24 # Number of cpus.
														   ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
														   ,lowB=1e-6 # Lower boundary for the rates.
														   ,upB=1e10 # Upper boundary for the rates.
														   ,FlagDev="FC" # Cost function.
														   ,lambda=0.05 # Regularization strength.
														   ,excludeSpecies=NULL # List of species to be excluded from the cost function.
														   ,parFixed=NULL) # List of parameters to be excluded from the optimization.
saveRDS(inferedRatesUntreatedMerged_yesChpPNoNp_multi,"inferedRatesUntreatedMerged_yesChpPNoNp_multi.rds")

inferedRatesUntreatedMerged_yesChpNoNpP_multi <- inferRates(expressionData=expressionLevels_yesChpNoNpP # Expression data of the genes to be modeled.
														   ,expressionDataDev=NULL # Standard deviations of the genes to be modeled.
														   ,simulatedDataset=NULL # Simulated dataset if this is the case (just to produce real rates correlations).
														   ,initialRates=initialRates # List of initial rates for optimization.
														   ,TauFractions=c(0,0.33) # Time points with cellular fractionation.
														   ,TauPoly=c(0,0.33) # Time points with polysomal profiling.
														   ,TauTotal=NULL # Time points with total RNA profiling.
														   ,cpus=24 # Number of cpus.
														   ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
														   ,lowB=1e-6 # Lower boundary for the rates.
														   ,upB=1e10 # Upper boundary for the rates.
														   ,FlagDev="FC" # Cost function.
														   ,lambda=0.05 # Regularization strength.
														   ,excludeSpecies=NULL # List of species to be excluded from the cost function.
														   ,parFixed=NULL) # List of parameters to be excluded from the optimization.
saveRDS(inferedRatesUntreatedMerged_yesChpNoNpP_multi,"inferedRatesUntreatedMerged_yesChpNoNpP_multi.rds")

inferedRatesUntreatedMerged_yesPNoChp_multi <- inferRates(expressionData=expressionLevels_yesPNoChp # Expression data of the genes to be modeled.
														 ,expressionDataDev=NULL # Standard deviations of the genes to be modeled.
														 ,simulatedDataset=NULL # Simulated dataset if this is the case (just to produce real rates correlations).
														 ,initialRates=initialRates # List of initial rates for optimization.
														 ,TauFractions=c(0,0.33) # Time points with cellular fractionation.
														 ,TauPoly=c(0,0.33) # Time points with polysomal profiling.
														 ,TauTotal=NULL # Time points with total RNA profiling.
														 ,cpus=24 # Number of cpus.
														 ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
														 ,lowB=1e-6 # Lower boundary for the rates.
														 ,upB=1e10 # Upper boundary for the rates.
														 ,FlagDev="FC" # Cost function.
														 ,lambda=0.05 # Regularization strength.
														 ,excludeSpecies=NULL # List of species to be excluded from the cost function.
														 ,parFixed=NULL) # List of parameters to be excluded from the optimization.
saveRDS(inferedRatesUntreatedMerged_yesPNoChp_multi,"inferedRatesUntreatedMerged_yesPNoChp_multi.rds")

inferedRatesUntreatedMerged_noChpP_multi <- inferRates(expressionData=expressionLevels_noChpP # Expression data of the genes to be modeled.
													  ,expressionDataDev=NULL # Standard deviations of the genes to be modeled.
													  ,simulatedDataset=NULL # Simulated dataset if this is the case (just to produce real rates correlations).
													  ,initialRates=initialRates # List of initial rates for optimization.
													  ,TauFractions=c(0,0.33) # Time points with cellular fractionation.
													  ,TauPoly=c(0,0.33) # Time points with polysomal profiling.
													  ,TauTotal=NULL # Time points with total RNA profiling.
													  ,cpus=24 # Number of cpus.
													  ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
													  ,lowB=1e-6 # Lower boundary for the rates.
													  ,upB=1e10 # Upper boundary for the rates.
													  ,FlagDev="FC" # Cost function.
													  ,lambda=0.05 # Regularization strength.
													  ,excludeSpecies=NULL # List of species to be excluded from the cost function.
													  ,parFixed=NULL) # List of parameters to be excluded from the optimization.
saveRDS(inferedRatesUntreatedMerged_noChpP_multi,"inferedRatesUntreatedMerged_noChpP_multi.rds")
