### External arguments extracted from simulateGene.sh
args = commandArgs(trailingOnly=TRUE)
 
for(v in args)
{
       vTmp <- strsplit(v,"=")[[1]]
       assign(vTmp[[1]],vTmp[[2]])
}

reps <- (unlist(strsplit(reps,",")))
TauFractions <- (unlist(strsplit(TauFractions,",")))
TauPoly <- (unlist(strsplit(TauPoly,",")))
Flags <- unlist(strsplit(Flags,","))
rates <- (unlist(strsplit(rates,",")))

reps <- as.numeric(reps)
CVs <- as.numeric(CVs)
translationCoeff <- as.numeric(translationCoeff)
TauFractions <- (as.numeric(TauFractions))
TauPoly <- (as.numeric(TauPoly))
nGenes <- as.numeric(nGenes)
cpus <- as.numeric(cpus)
ZeroThresh <- as.numeric(ZeroThresh)
MultFact <- as.numeric(MultFact)

### Libraries
library(parallel)
library(deSolve)
library(pheatmap)

source("/path/to/allInternalFunctions.R")
# source("path/to/allInternalFunctionsPatch_nuclearPrematureDecay.R") # Execute this line to simulate data with Nucleoplasmic Premature RNA degradation.
### Average gene
INSPEcT_nascent <- read.table("/path/to/regulated_genes_features.xls",sep="\t",header=TRUE)

### Median kinetic rates
synthesis=signif(median(INSPEcT_nascent[,"synthesis_0"]),2)
processing=signif(median(INSPEcT_nascent[,"processing_0"]),2)
degradation=signif(median(INSPEcT_nascent[,"degradation_0"]),2)

### Mean rates
exampleRates <- c(k1=synthesis
                 ,k2=processing
                 ,k3=processing
                 ,k4=0.1*processing
                 ,k5=0.25*processing
                 ,k6=1.5*processing
                 ,k7=translationCoeff*degradation
                 ,k8=degradation
                 ,k9=0.5*degradation
                 ,k10=degradation)

initialRates <- list(10**ceiling(log10(exampleRates))*(1+seq_along(exampleRates)*1e-5)
                    ,10**round(log10(exampleRates))*(1+seq_along(exampleRates)*1e-5)
                    ,10**floor(log10(exampleRates))*(1+seq_along(exampleRates)*1e-5))

exampleRates <- exampleRates[names(exampleRates)%in%rates]
initialRates <- lapply(initialRates,function(i)i[names(i)%in%rates])

### Simulation (loop)
## Number of replicates
for(rep in reps)
{
    print(rep)
    
    ## Variation coefficient
    for(CV in CVs)
    {
        print(CV)

        simulatedDataset <- simulateData(exampleRates=exampleRates # Mean rates to sample.
                                        ,TauFractions=TauFractions # Time points with cellular fractionation.
                                        ,TauPoly=TauPoly # Time points with polysomal profiling.
                                        ,TauTotal=NULL # Time points with total RNA profiling (FIXED to null for our purposes).
                                        ,noise=TRUE # TRUE to add noise to data.
                                        ,CV=CV # Variation Coefficient for noise.
                                        ,Reps=rep # Number of replicates.
                                        ,nGenes=nGenes # Number of genes to be simulated.
                                        ,seed=1 # Seed for reproducibility.
                                        ,ZeroThresh=ZeroThresh # Minimum expression value.
                                        ,MultFact=MultFact) # Number of nGenes to be modeled.
        ## Save simulated data
        # saveRDS(simulatedDataset,file=paste0("simulatedDataset_",simulatedDataset$name,".rds"))

        ## Cost functions
        for(flag in Flags){
            
            ## Inference
            inferedNumericalModel <- inferRates(expressionData=simulatedDataset$exampleData # Expression data of the genes to be modeled.
                                               ,expressionDataDev=simulatedDataset$DevDataTmp # Standard deviations of the genes to be modeled.
                                               ,simulatedDataset=simulatedDataset # Simulated dataset if this is the case (just to produce real rates correlations).
                                               ,initialRates=initialRates # List of initial rates for optimization.
                                               ,TauFractions=TauFractions # Time points with cellular fractionation.
                                               ,TauPoly=TauPoly # Time points with polysomal profiling.
                                               ,TauTotal=NULL # Time points with total RNA profiling.
                                               ,cpus=cpus # Number of cpus.
                                               ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
                                               ,lowB=1e-6 # Lower boundary for the rates.
                                               ,upB=1e4 # Upper boundary for the rates.
                                               ,FlagDev=flag # Cost function.
                                               ,lambda=0.05 # Regularization strength.
                                               ,excludeSpecies=NULL # List of species to be excluded from the cost function.
                                               ,parFixed=NULL) # List of parameters to be excluded from the optimization.
            
            ## Save inferred data
            saveRDS(inferedNumericalModel,file=paste0("inferedNumericalModel_",simulatedDataset$name,".rds"))
        }
    }
}

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
