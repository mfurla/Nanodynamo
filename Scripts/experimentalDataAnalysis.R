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

sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-conda-linux-gnu (64-bit)
# Running under: Rocky Linux 8.7 (Green Obsidian)
# 
# Matrix products: default
# BLAS/LAPACK:
# 
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
# [1] TxDb.Hsapiens.UCSC.hg38.knownGene_3.16.0
# [2] GenomicFeatures_1.50.4                  
# [3] AnnotationDbi_1.60.2                    
# [4] Biobase_2.58.0                          
# [5] GenomicRanges_1.50.2                    
# [6] GenomeInfoDb_1.34.9                     
# [7] IRanges_2.32.0                          
# [8] S4Vectors_0.36.2                        
# [9] BiocGenerics_0.44.0                     
# 
# loaded via a namespace (and not attached):
#  [1] lattice_0.21-8              prettyunits_1.1.1          
#  [3] png_0.1-8                   Rsamtools_2.14.0           
#  [5] Biostrings_2.66.0           digest_0.6.33              
#  [7] utf8_1.2.3                  BiocFileCache_2.6.1        
#  [9] R6_2.5.1                    RSQLite_2.3.1              
# [11] httr_1.4.6                  pillar_1.9.0               
# [13] zlibbioc_1.44.0             rlang_1.1.1                
# [15] progress_1.2.2              curl_4.3.2                 
# [17] blob_1.2.4                  Matrix_1.6-0               
# [19] BiocParallel_1.32.6         stringr_1.5.0              
# [21] RCurl_1.98-1.12             bit_4.0.5                  
# [23] biomaRt_2.54.1              DelayedArray_0.24.0        
# [25] compiler_4.2.1              rtracklayer_1.58.0         
# [27] pkgconfig_2.0.3             tidyselect_1.2.0           
# [29] KEGGREST_1.38.0             SummarizedExperiment_1.28.0
# [31] tibble_3.2.1                GenomeInfoDbData_1.2.9     
# [33] codetools_0.2-19            matrixStats_1.0.0          
# [35] XML_3.99-0.14               fansi_1.0.4                
# [37] crayon_1.5.2                dplyr_1.1.2                
# [39] dbplyr_2.3.3                GenomicAlignments_1.34.1   
# [41] bitops_1.0-7                rappdirs_0.3.3             
# [43] grid_4.2.1                  lifecycle_1.0.3            
# [45] DBI_1.1.3                   magrittr_2.0.3             
# [47] cli_3.6.1                   stringi_1.7.12             
# [49] cachem_1.0.8                XVector_0.38.0             
# [51] xml2_1.3.3                  filelock_1.0.2             
# [53] generics_0.1.3              vctrs_0.6.3                
# [55] rjson_0.2.21                restfulr_0.0.15            
# [57] tools_4.2.1                 bit64_4.0.5                
# [59] glue_1.6.2                  hms_1.1.3                  
# [61] MatrixGenerics_1.10.0       parallel_4.2.1             
# [63] fastmap_1.1.1               yaml_2.3.7                 
# [65] memoise_2.0.1               BiocIO_1.8.0               
