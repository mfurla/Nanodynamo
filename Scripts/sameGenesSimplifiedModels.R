### Libraries
library("GenomicAlignments")
require("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library("dplyr")
library("deSolve")
library("parallel")
library("ggplot2")
library("gridExtra")
library("pheatmap")
library("clusterProfiler")
library("org.Hs.eg.db")
library("biomaRt")

source("/path/to/allInternalFunctions.R")

### Expression levels
expressionLevelsMerged <- readRDS("/path/to/Results/Untreated/expressionLevelsMerged.rds")

### Modeling
## Initial rates
INSPEcT_nascent <- read.table("/path/to/regulated_genes_features.xls",sep="\t",header=TRUE)

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
expressionLevels_YesChpNpP <- expressionLevelsMerged$normalizedCountsYesChpNpP[["Sample_Merged"]]
expressionLevels_YesChpNpNoP <- expressionLevels_YesChpNpP[,!(colnames(expressionLevels_YesChpNpP)%in%c("pyp0.33","pyn0.33","pyn0","pyp0","pyt"))]

expressionLevels_YesChpPNoNp <- expressionLevels_YesChpNpP
expressionLevels_YesChpNoNpP <- expressionLevels_YesChpNpNoP

expressionLevels_YesChpPNoNp[,c("npp0.33","npn0.33","npn0","npp0","npt")] <- expressionLevels_YesChpNpP[,c("npp0.33","npn0.33","npn0","npp0","npt")] + expressionLevels_YesChpNpP[,c("nmp0.33","nmn0.33","nmn0","nmp0","nmt")]
expressionLevels_YesChpPNoNp <- expressionLevels_YesChpPNoNp[,!(colnames(expressionLevels_YesChpPNoNp)%in%c("nmp0.33","nmn0.33","nmn0","nmp0","nmt"))]
colnames(expressionLevels_YesChpPNoNp) <- gsub("^np","n",colnames(expressionLevels_YesChpPNoNp))

expressionLevels_YesChpNoNpP[,c("npp0.33","npn0.33","npn0","npp0","npt")] <- expressionLevels_YesChpNpNoP[,c("npp0.33","npn0.33","npn0","npp0","npt")] + expressionLevels_YesChpNpNoP[,c("nmp0.33","nmn0.33","nmn0","nmp0","nmt")]
expressionLevels_YesChpNoNpP <- expressionLevels_YesChpNoNpP[,!(colnames(expressionLevels_YesChpNoNpP)%in%c("nmp0.33","nmn0.33","nmn0","nmp0","nmt"))]
colnames(expressionLevels_YesChpNoNpP) <- gsub("^np","n",colnames(expressionLevels_YesChpNoNpP))

expressionLevels_YesPNoChp <- expressionLevels_YesChpPNoNp
expressionLevels_NoChpP <- expressionLevels_YesChpNoNpP

expressionLevels_YesPNoChp[,c("chpp0.33","chpn0.33","chpn0","chpp0","chpt")] <- expressionLevels_YesChpPNoNp[,c("chpp0.33","chpn0.33","chpn0","chpp0","chpt")] + expressionLevels_YesChpPNoNp[,c("chmp0.33","chmn0.33","chmn0","chmp0","chmt")]
expressionLevels_YesPNoChp <- expressionLevels_YesPNoChp[,!(colnames(expressionLevels_YesPNoChp)%in%c("chmp0.33","chmn0.33","chmn0","chmp0","chmt"))]
colnames(expressionLevels_YesPNoChp) <- gsub("^chp","ch",colnames(expressionLevels_YesPNoChp))

expressionLevels_NoChpP[,c("chpp0.33","chpn0.33","chpn0","chpp0","chpt")] <- expressionLevels_YesChpNoNpP[,c("chpp0.33","chpn0.33","chpn0","chpp0","chpt")] + expressionLevels_YesChpNoNpP[,c("chmp0.33","chmn0.33","chmn0","chmp0","chmt")]
expressionLevels_NoChpP <- expressionLevels_NoChpP[,!(colnames(expressionLevels_NoChpP)%in%c("chmp0.33","chmn0.33","chmn0","chmp0","chmt"))]
colnames(expressionLevels_NoChpP) <- gsub("^chp","ch",colnames(expressionLevels_NoChpP))

## Inference
inferedRatesUntreatedMerged_yesChpNpNoP <- inferRates(expressionData=expressionLevels_YesChpNpNoP # Expression data of the genes to be modeled.
													 ,expressionDataDev=NULL # Standard deviations of the genes to be modeled.
													 ,simulatedDataset=NULL # Simulated dataset if this is the case (just to produce real rates correlations).
													 ,initialRates=initialRates # List of initial rates for optimization.
													 ,TauFractions=c(0,0.33) # Time points with cellular fractionation.
													 ,TauPoly=c(0,0.33) # Time points with polysomal profiling.
													 ,TauTotal=NULL # Time points with total RNA profiling.
													 ,cpus=24 # Number of cpus.
													 ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
													 ,lowB=1e-6 # Lower boundary for the rates.
													 ,upB=1e6 # Upper boundary for the rates.
													 ,FlagDev="FC" # Cost function.
													 ,lambda=0.05 # Regularization strength.
													 ,excludeSpecies=NULL # List of species to be excluded from the cost function.
													 ,parFixed=NULL) # List of parameters to be excluded from the optimization.

inferedRatesUntreatedMerged_yesChpPNoNp <- inferRates(expressionData=expressionLevels_YesChpPNoNp # Expression data of the genes to be modeled.
													 ,expressionDataDev=NULL # Standard deviations of the genes to be modeled.
													 ,simulatedDataset=NULL # Simulated dataset if this is the case (just to produce real rates correlations).
													 ,initialRates=initialRates # List of initial rates for optimization.
													 ,TauFractions=c(0,0.33) # Time points with cellular fractionation.
													 ,TauPoly=c(0,0.33) # Time points with polysomal profiling.
													 ,TauTotal=NULL # Time points with total RNA profiling.
													 ,cpus=24 # Number of cpus.
													 ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
													 ,lowB=1e-6 # Lower boundary for the rates.
													 ,upB=1e6 # Upper boundary for the rates.
													 ,FlagDev="FC" # Cost function.
													 ,lambda=0.05 # Regularization strength.
													 ,excludeSpecies=NULL # List of species to be excluded from the cost function.
													 ,parFixed=NULL) # List of parameters to be excluded from the optimization.

inferedRatesUntreatedMerged_yesChpNoNpP <- inferRates(expressionData=expressionLevels_YesChpNoNpP # Expression data of the genes to be modeled.
													 ,expressionDataDev=NULL # Standard deviations of the genes to be modeled.
													 ,simulatedDataset=NULL # Simulated dataset if this is the case (just to produce real rates correlations).
													 ,initialRates=initialRates # List of initial rates for optimization.
													 ,TauFractions=c(0,0.33) # Time points with cellular fractionation.
													 ,TauPoly=c(0,0.33) # Time points with polysomal profiling.
													 ,TauTotal=NULL # Time points with total RNA profiling.
													 ,cpus=24 # Number of cpus.
													 ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
													 ,lowB=1e-6 # Lower boundary for the rates.
													 ,upB=1e6 # Upper boundary for the rates.
													 ,FlagDev="FC" # Cost function.
													 ,lambda=0.05 # Regularization strength.
													 ,excludeSpecies=NULL # List of species to be excluded from the cost function.
													 ,parFixed=NULL) # List of parameters to be excluded from the optimization.

inferedRatesUntreatedMerged_yesPNoChp <- inferRates(expressionData=expressionLevels_YesPNoChp # Expression data of the genes to be modeled.
												   ,expressionDataDev=NULL # Standard deviations of the genes to be modeled.
												   ,simulatedDataset=NULL # Simulated dataset if this is the case (just to produce real rates correlations).
												   ,initialRates=initialRates # List of initial rates for optimization.
												   ,TauFractions=c(0,0.33) # Time points with cellular fractionation.
												   ,TauPoly=c(0,0.33) # Time points with polysomal profiling.
												   ,TauTotal=NULL # Time points with total RNA profiling.
												   ,cpus=24 # Number of cpus.
												   ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
												   ,lowB=1e-6 # Lower boundary for the rates.
												   ,upB=1e6 # Upper boundary for the rates.
												   ,FlagDev="FC" # Cost function.
												   ,lambda=0.05 # Regularization strength.
												   ,excludeSpecies=NULL # List of species to be excluded from the cost function.
												   ,parFixed=NULL) # List of parameters to be excluded from the optimization.

inferedRatesUntreatedMerged_noChpP <- inferRates(expressionData=expressionLevels_NoChpP # Expression data of the genes to be modeled.
												,expressionDataDev=NULL # Standard deviations of the genes to be modeled.
												,simulatedDataset=NULL # Simulated dataset if this is the case (just to produce real rates correlations).
												,initialRates=initialRates # List of initial rates for optimization.
												,TauFractions=c(0,0.33) # Time points with cellular fractionation.
												,TauPoly=c(0,0.33) # Time points with polysomal profiling.
												,TauTotal=NULL # Time points with total RNA profiling.
												,cpus=24 # Number of cpus.
												,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
												,lowB=1e-6 # Lower boundary for the rates.
												,upB=1e6 # Upper boundary for the rates.
												,FlagDev="FC" # Cost function.
												,lambda=0.05 # Regularization strength.
												,excludeSpecies=NULL # List of species to be excluded from the cost function.
												,parFixed=NULL) # List of parameters to be excluded from the optimization.

### Analyses - Supplemental Methods Figures
inferedRatesUntreatedMerged_YesChpNpP <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreatedMerged_yesChpNpP_multi.rds")

mains=c(k1="Synthesis\n[fg/(M. cells*h)]",k3="Detachment of\nmature chromatin\n[h^-1]",k6="Export\n[h^-1]",k8="Degradation of\nuntranslated RNAs\n[h^-1]")
units=c(k1="fg/(M. cells*h)",k3="h^-1",k6="h^-1",k8="h^-1")

pheatmap(cor(log10(cbind(Chp_Np_P=inferedRatesUntreatedMerged_YesChpNpP$inferedRates[,"k1"]
						,Chp_Np=inferedRatesUntreatedMerged_YesChpNpNoP$inferedRates[,"k1"]
						,Chp_P=inferedRatesUntreatedMerged_YesChpPNoNp$inferedRates[,"k1"]
						,Chp=inferedRatesUntreatedMerged_YesChpNoNpP$inferedRates[,"k1"]
						,P=inferedRatesUntreatedMerged_YesPNoChp$inferedRates[,"k1"]
						,"-"=inferedRatesUntreatedMerged_NoChpP$inferedRates[,"k1"])),method="s")
	   ,breaks=seq(0,1,length.out=20)
	   ,color=colorRampPalette(c("white","red"))(20)
	   ,width=3.5
	   ,height=3.5
	   ,filename="synthesisHeatmap.pdf"
	   ,display_numbers=TRUE
	   ,cluster_rows=FALSE
	   ,cluster_cols=FALSE
	   ,main="synthesis"
	   ,legend_labels="Spearman correlation")

pheatmap(cor(log10(cbind(Chp_Np_P=inferedRatesUntreatedMerged_YesChpNpP$inferedRates[,"k8"]
						,Chp_Np=inferedRatesUntreatedMerged_YesChpNpNoP$inferedRates[,"k8"]
						,Chp_P=inferedRatesUntreatedMerged_YesChpPNoNp$inferedRates[,"k8"]
						,Chp=inferedRatesUntreatedMerged_YesChpNoNpP$inferedRates[,"k8"]
						,P=inferedRatesUntreatedMerged_YesPNoChp$inferedRates[,"k8"]
						,"-"=inferedRatesUntreatedMerged_NoChpP$inferedRates[,"k8"])),method="s")
	   ,breaks=seq(0,1,length.out=20)
	   ,color=colorRampPalette(c("white","red"))(20)
	   ,width=3.5
	   ,height=3.5
	   ,filename="degradationHeatmap.pdf"
	   ,display_numbers=TRUE
	   ,cluster_rows=FALSE
	   ,cluster_cols=FALSE
	   ,main="degradation"
	   ,legend_labels="Spearman correlation")

pheatmap(cor(log10(cbind(Chp_Np_P=inferedRatesUntreatedMerged_YesChpNpP$inferedRates[,"k6"]
						,Chp_Np=inferedRatesUntreatedMerged_YesChpNpNoP$inferedRates[,"k6"]
						,Chp_P=inferedRatesUntreatedMerged_YesChpPNoNp$inferedRates[,"k6"]
						,Chp=inferedRatesUntreatedMerged_YesChpNoNpP$inferedRates[,"k6"]
						,P=inferedRatesUntreatedMerged_YesPNoChp$inferedRates[,"k6"]
						,"-"=inferedRatesUntreatedMerged_NoChpP$inferedRates[,"k6"])),method="s")
	   ,breaks=seq(0,1,length.out=20)
	   ,color=colorRampPalette(c("white","red"))(20)
	   ,width=3.5
	   ,height=3.5
	   ,filename="exportHeatmap.pdf"
	   ,display_numbers=TRUE
	   ,cluster_rows=FALSE
	   ,cluster_cols=FALSE
	   ,main="export"
	   ,legend_labels="Spearman correlation")

pheatmap(cor(log10(cbind(Chp_Np_P=inferedRatesUntreatedMerged_YesChpNpP$inferedRates[,"k3"]
						,Chp_Np=inferedRatesUntreatedMerged_YesChpNpNoP$inferedRates[,"k3"]
						,Chp_P=inferedRatesUntreatedMerged_YesChpPNoNp$inferedRates[,"k3"]
						,Chp=inferedRatesUntreatedMerged_YesChpNoNpP$inferedRates[,"k3"]
						,P=inferedRatesUntreatedMerged_YesPNoChp$inferedRates[,"k3"]
						,"-"=inferedRatesUntreatedMerged_NoChpP$inferedRates[,"k3"])),method="s")
	   ,breaks=seq(0,1,length.out=20)
	   ,color=colorRampPalette(c("white","red"))(20)
	   ,width=3.5
	   ,height=3.5
	   ,filename="detachmentHeatmap.pdf"
	   ,display_numbers=TRUE
	   ,cluster_rows=FALSE
	   ,cluster_cols=FALSE
	   ,main="detachment"
	   ,legend_labels="Spearman correlation")

