### Libraries
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
tbl_treat <- cbind("cells"=# Number of cells, example for the Untreated condition: c(3,3,3,2.8,2.8,2.8,29,13)
			   ,"PolyA"=# Amount of polyA RNA, example for the Untreated condition: c(128.4,122,335,109.2,88.48,308,1298.5,516)
			   )
rownames(tbl_treat) <- # Experimental design for yield estimation, example for the Untreated condition: c("Chr_WTA","Nuc_WTA","Cyt_WTA","Chr_WTB","Nuc_WTB","Cyt_WTB","Poly_WTC","Poly_WTD")

### Expression data
expressionLevels <- expressionDataEstimation(bamPaths=bamPaths_treat
											,nascentPaths=nascentPaths_treat
											,tbl=tbl_treat
											,labelingTime=0.33
											,labelingTimePoly=0.33
											,txdb=txdb
											,minoverlap_I=10
											,minoverlap_E=10
											,saveDataDistributions=FALSE
											,cpus=1)

### Single initial condition modeling
initialRates <- list(c(k1=1,k2=1,k3=1,k4=1,k5=1,k6=1,k7=1,k8=1,k10=1))

inferedRatesTreatMerged_YesChpNpP_single <- inferRates(expressionData=list(expressionLevels$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT1$mean
																	  ,expressionLevels$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT2$mean) # Expression data of the genes to be modeled.
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
saveRDS(inferedRatesTreatMerged_YesChpNpP_single,"inferedRatesTreatMerged_YesChpNpP_single.rds")

# Genes sub-sampling based on replicate 1
expressionDataListUntreated1 <- expressionLevels$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT1$mean
expressionDataListUntreated1 <- expressionDataListUntreated1[apply(expressionDataListUntreated1[,!grepl("n0$",colnames(expressionDataListUntreated1))]>1e-10,1,all),]

inferedRatesTreat1_YesChpNpP_single <- inferRates(expressionData=list(expressionDataListUntreated1) # Expression data of the genes to be modeled.
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
saveRDS(inferedRatesTreat1_YesChpNpP_single,"inferedRatesTreat1_YesChpNpP_single.rds")

# Genes sub-sampling based on replicate 2
expressionDataListUntreated2 <- expressionLevels$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT2$mean
expressionDataListUntreated2 <- expressionDataListUntreated2[apply(expressionDataListUntreated2[,!grepl("n0$",colnames(expressionDataListUntreated2))]>1e-10,1,all),]

inferedRatesTreat2_YesChpNpP_single <- inferRates(expressionData=list(expressionDataListUntreated2) # Expression data of the genes to be modeled.
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
saveRDS(inferedRatesTreat2_YesChpNpP_single,"inferedRatesTreat2_YesChpNpP_single.rds")

### Multiple initial condition modeling
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

inferedRatesTreatMerged_YesChpNpP_multi <- inferRates(expressionData=list(expressionLevels$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT1$mean
																		 ,expressionLevels$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT2$mean) # Expression data of the genes to be modeled.
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
saveRDS(inferedRatesTreatMerged_YesChpNpP_multi,"inferedRatesTreatMerged_YesChpNpP_multi.rds")

# Genes sub-sampling based on replicate 1
expressionDataListUntreated1 <- expressionLevels$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT1$mean
expressionDataListUntreated1 <- expressionDataListUntreated1[apply(expressionDataListUntreated1[,!grepl("n0$",colnames(expressionDataListUntreated1))]>1e-10,1,all),]

inferedRatesTreat1_YesChpNpP_multi <- inferRates(expressionData=list(expressionDataListUntreated1) # Expression data of the genes to be modeled.
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
saveRDS(inferedRatesTreat1_YesChpNpP_multi,"inferedRatesTreat1_YesChpNpP_multi.rds")

# Genes sub-sampling based on replicate 2
expressionDataListUntreated2 <- expressionLevels$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT2$mean
expressionDataListUntreated2 <- expressionDataListUntreated2[apply(expressionDataListUntreated2[,!grepl("n0$",colnames(expressionDataListUntreated2))]>1e-10,1,all),]

inferedRatesTreat2_YesChpNpP_multi <- inferRates(expressionData=list(expressionDataListUntreated2) # Expression data of the genes to be modeled.
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
saveRDS(inferedRatesTreat2_YesChpNpP_multi,"inferedRatesTreat2_YesChpNpP_multi.rds")