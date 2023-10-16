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
# source("/path/to/allInternalFunctionsPatch_nucleoplasmicPrematureExport.R") # Execute this line to model Nucleoplasmic Premature RNA export; available for the Full Model only.
### BAMs_treat
bamPaths_treat <- c(Chr_treat1="/path/to/bam.RData"
				   ,Nuc_treat1="/path/to/bam.RData"
				   ,Cyt_treat1="/path/to/bam.RData"
				   ,Poly_treat1="/path/to/bam.RData"
				   ,Chr_treat2="/path/to/bam.RData"
				   ,Nuc_treat2="/path/to/bam.RData"
				   ,Cyt_treat2="/path/to/bam.RData"
				   ,Poly_treat2="/path/to/bam.RData"
				   )

### Nascent_treat
nascentPaths_treat <- c(Chr_treat1="/path/to/nanoID.modification.probabilities.RData"
					   ,Nuc_treat1="/path/to/nanoID.modification.probabilities.RData"
					   ,Cyt_treat1="/path/to/nanoID.modification.probabilities.RData"
					   ,Poly_treat1="/path/to/nanoID.modification.probabilities.RData"
					   ,Chr_treat2="/path/to/nanoID.modification.probabilities.RData"
					   ,Nuc_treat2="/path/to/nanoID.modification.probabilities.RData"
					   ,Cyt_treat2="/path/to/nanoID.modification.probabilities.RData"
					   ,Poly_treat2="/path/to/nanoID.modification.probabilities.RData"
					   )

### FASTQs_treat
fastqs_treat <-  c(Chr_treat1="/path/to/FASTQ_DNA.fastq"
				  ,Nuc_treat1="/path/to/FASTQ_DNA.fastq"
				  ,Cyt_treat1="/path/to/FASTQ_DNA.fastq"
				  ,Poly_treat1="/path/to/FASTQ_DNA.fastq"
				  ,Chr_treat2="/path/to/FASTQ_DNA.fastq"
				  ,Nuc_treat2="/path/to/FASTQ_DNA.fastq"
				  ,Cyt_treat2="/path/to/FASTQ_DNA.fastq"
				  ,Poly_treat2="/path/to/FASTQ_DNA.fastq"
				  )

### Normalization data
## Number of reads
nr_treat = c()
for(i in 1:length(fastqs_treat)){nr_treat[i] <- system(command = paste0("cat ", fastqs_treat[i], " | /path/to/miniconda3/bin/seqtk seq -A - | grep \"^>\" | wc -l"), intern = TRUE)}
nr_treat <- as.numeric(nr_treat)

## Cells and RNA
tbl_treat <- matrix(0,nrow=length(bamPaths_treat),ncol=3)
rownames(tbl_treat) <- names(fastqs_treat)
colnames(tbl_treat) <- c("cells","PolyA","reads")
tbl_treat[,3] <- nr_treat
tbl_treat[,2] <- # Amount of polyA RNA. Example for Untreated condition: c(247.5,268.6,302.5,1298.5,684,1600,1134.2,516)
tbl_treat[,1] <- c(16,16,16,29,21,21,21,13)

### Experimental design
expDesign <- c(Chr_treat1="chr"
			  ,Nuc_treat1="nuc"
			  ,Cyt_treat1="cyt"
			  ,Poly_treat1="poly"
			  ,Chr_treat2="chr"
			  ,Nuc_treat2="nuc"
			  ,Cyt_treat2="cyt"
			  ,Poly_treat2="poly"
)

### Expression data
## Mean of replicates
expressionLevels <- expressionDataEstimation(bamPaths=bamPaths_treat
											,nascentPaths=nascentPaths_treat
											,expDesign=expDesign
											,tbl=tbl_treat
											,labelingTime=0.33
											,labelingTimePoly=0.33
											,txdb=txdb
											,countsTh=TRUE
											,cytPTh=NULL
											,minoverlap_I=10
											,minoverlap_E=10
											,spikeInsConcentrations=NULL
											,yeastCounts=c(Chr_treat1=127234,Nuc_treat1=54340,Cyt_treat1=73272) # This is required only for the Sample 1 of the Untreated condition because it was aligned to the human genome without ENO2, for all the other samples these data are inferred from the bam; set to NULL.
											,mergeSamples=FALSE)

## Merge of replicates
expressionLevelsMerged <- expressionDataEstimation(bamPaths=bamPaths_treat
												  ,nascentPaths=nascentPaths_treat
												  ,expDesign=expDesign
												  ,tbl=tbl_treat
												  ,labelingTime=0.33
												  ,labelingTimePoly=0.33
												  ,txdb=txdb
												  ,countsTh=TRUE
												  ,cytPTh=NULL
												  ,minoverlap_I=10
												  ,minoverlap_E=10
												  ,spikeInsConcentrations=NULL
												  ,yeastCounts=c(Chr_treat1=127234,Nuc_treat1=54340,Cyt_treat1=73272) # This is required only for the Sample 1 of the Untreated condition because it was aligned to the human genome without ENO2, for all the other samples these data are inferred from the bam; set to NULL.
												  ,mergeSamples=TRUE)

### Single initial condition modeling
initialRates <- list(c(k1=1,k2=1,k3=1,k4=1,k5=1,k6=1,k7=1,k8=1,k10=1))

inferedRatestreat1_yesChpNpP_single <- inferRates(expressionData=expressionLevels$normalizedCountsYesChpNpP[["Sample_treat1"]] # Expression data of the genes to be modeled.
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

inferedRatestreat2_yesChpNpP_single <- inferRates(expressionData=expressionLevels$normalizedCountsYesChpNpP[["Sample_treat2"]] # Expression data of the genes to be modeled.
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

inferedRatestreatMerged_yesChpNpP_single <- inferRates(expressionData=expressionLevelsMerged$normalizedCountsYesChpNpP[["Sample_Merged"]] # Expression data of the genes to be modeled.
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

### Multiple initial condition modeling
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

inferedRatestreat1_yesChpNpP_multi <- inferRates(expressionData=expressionLevels$normalizedCountsYesChpNpP[["Sample_treat1"]] # Expression data of the genes to be modeled.
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

inferedRatestreat2_yesChpNpP_multi <- inferRates(expressionData=expressionLevels$normalizedCountsYesChpNpP[["Sample_treat2"]] # Expression data of the genes to be modeled.
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

inferedRatestreatMerged_yesChpNpP_multi <- inferRates(expressionData=expressionLevelsMerged$normalizedCountsYesChpNpP[["Sample_Merged"]] # Expression data of the genes to be modeled.
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