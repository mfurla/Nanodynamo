### Libraries
library("INSPEcT")
library("pheatmap")

### Expression data
expressionLevelsUntreated <- readRDS("/path/to/Results/Untreated/expressionLevelsMerged.rds")

### NanoID results
nanoID <- readRDS("/path/to/Results/nanoIDRates.rds")

### Data aggregation
mergedExpressionLevelsUntreated <- sapply(c("normalizedCountsYesChpNpP"
										   ,"normalizedCountsYesChpPNoNp"
										   ,"normalizedCountsYesPNoChp"
										   ,"normalizedCountsYesChpNpNoP"
										   ,"normalizedCountsYesChpNoNpP"
										   ,"normalizedCountsNoChpP")
,function(i){
	print(i)

	j <- expressionLevelsUntreated[[i]][["Sample_Merged"]]
	j <- j[,!grepl("^py",colnames(j))]
	j <- j[,!grepl("0$",colnames(j))]
	j <- j[,!grepl("t$",colnames(j))]

	pTmp <- j[,grep("p0.33",colnames(j))]
	nTmp <- j[,grep("n0.33",colnames(j))]

	ppTmp <- tryCatch(apply(as.matrix(j[,grep("pp0.33",colnames(j))],nrow=nrow(j)),1,sum),error=function(e){rep(0,nrow(j))})
	pnTmp <- tryCatch(apply(as.matrix(j[,grep("pn0.33",colnames(j))],nrow=nrow(j)),1,sum),error=function(e){rep(0,nrow(j))})
	
	cbind("pt"=ppTmp+pnTmp,"pn"=pnTmp,"tt"=apply(pTmp,1,sum)+apply(nTmp,1,sum),"tn"=apply(nTmp,1,sum))
})
mergedExpressionLevelsUntreated <- do.call("rbind",mergedExpressionLevelsUntreated)

### INSPEcT
## Expression levels quantification
nascent <- quantifyExpressionsFromTrAbundance(trAbundaces = list(exonsAbundances=cbind("Rep1"=mergedExpressionLevelsUntreated[,"tn"]
																					  ,"Rep2"=1e-10*mergedExpressionLevelsUntreated[,"tn"] # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																					  ,"Rep3"=mergedExpressionLevelsUntreated[,"tn"]
																					  ,"Rep4"=1e-10*mergedExpressionLevelsUntreated[,"tn"]) # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																,intronsAbundances=cbind("Rep1"=mergedExpressionLevelsUntreated[,"pn"]
																						,"Rep2"=1e-10*mergedExpressionLevelsUntreated[,"pn"] # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																						,"Rep3"=mergedExpressionLevelsUntreated[,"pn"]
																						,"Rep4"=1e-10*mergedExpressionLevelsUntreated[,"pn"])) # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
											, experimentalDesign = c(Rep1="Untreated",Rep2="Untreated",Rep3="UntreatedBis",Rep4="UntreatedBis"))

total <- quantifyExpressionsFromTrAbundance(trAbundaces = list(exonsAbundances=cbind("Rep1"=mergedExpressionLevelsUntreated[,"tt"]
																					,"Rep2"=1e-10*mergedExpressionLevelsUntreated[,"tt"] # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																					,"Rep3"=mergedExpressionLevelsUntreated[,"tt"]
																					,"Rep4"=1e-10*mergedExpressionLevelsUntreated[,"tt"]) # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																,intronsAbundances=cbind("Rep1"=mergedExpressionLevelsUntreated[,"pt"]
																						,"Rep2"=1e-10*mergedExpressionLevelsUntreated[,"pt"] # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																						,"Rep3"=mergedExpressionLevelsUntreated[,"pt"]
																						,"Rep4"=1e-10*mergedExpressionLevelsUntreated[,"pt"])) # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
											, experimentalDesign = c(Rep1="Untreated",Rep2="Untreated",Rep3="UntreatedBis",Rep4="UntreatedBis"))

## Main INSPEcT function
nascentInspObj <- newINSPEcT(tpts=c("Untreated","UntreatedBis")
							,labeling_time=20/60
							,nascentExpressions=nascent
							,matureExpressions=total)

### Comparative analysis
## INSPEcT rates
INSPEcT_synthesis <- log10(ratesFirstGuess(nascentInspObj,"synthesis")[,"synthesis_Untreated"])
INSPEcT_processing <- log10(ratesFirstGuess(nascentInspObj,"processing")[,"processing_Untreated"])
INSPEcT_degradation <- log10(ratesFirstGuess(nascentInspObj,"degradation")[,"degradation_Untreated"])

## Nanodynamo rates
inferedRatesUntreated_NoChpP <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_noChpP.rds")
inferedRatesUntreated_YesChpNoNpP <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_yesChpNoNpP.rds")
inferedRatesUntreated_YesChpNpNoP <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_yesChpNpNoP.rds")
inferedRatesUntreated_YesChpNpP <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreatedMerged_yesChpNpP_multi.rds")
inferedRatesUntreated_YesChpPNoNp <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_yesChpPNoNp.rds")
inferedRatesUntreated_YesPNoChp <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_yesPNoChp.rds")

Nanodynamo_synthesis <- list("WithP"=log10(c(inferedRatesWT_YesChpNpP$inferedRates[,"k1"]
											,inferedRatesWT_YesChpPNoNp$inferedRates[,"k1"]
											,inferedRatesWT_YesPNoChp$inferedRates[,"k1"]))
							,"WithOutP"=log10(c(inferedRatesWT_YesChpNpNoP$inferedRates[,"k1"]
											,inferedRatesWT_YesChpNoNpP$inferedRates[,"k1"]
											,inferedRatesWT_NoChpP$inferedRates[,"k1"]))
							,"NanoID"=log10(nanoID[,"Synthesis"]))

Nanodynamo_degradation <- list("WithP"=log10(c(inferedRatesUntreated_YesChpNpP$inferedRates[,"k8"]
						  					  ,inferedRatesUntreated_YesChpPNoNp$inferedRates[,"k8"]
						  					  ,inferedRatesUntreated_YesPNoChp$inferedRates[,"k8"]))
							  ,"WithOutP"=log10(c(inferedRatesUntreated_YesChpNpNoP$inferedRates[,"k8"]
											     ,inferedRatesUntreated_YesChpNoNpP$inferedRates[,"k8"]
											     ,inferedRatesUntreated_NoChpP$inferedRates[,"k8"]))
							  ,"NanoID"=log10(nanoID[,"Degradation"]))

sapply(names(Nanodynamo_synthesis),function(j)
{
	i <- Nanodynamo_synthesis[[j]]
	commonGenesTmp <- intersect(names(i),names(INSPEcT_synthesis))
	print(length(commonGenesTmp))
	round(cor(INSPEcT_synthesis[commonGenesTmp],i[commonGenesTmp],method="s",use="c"),2)
})

# WithP WithOutP NanoID
#  0.95     0.81   0.99

sapply(names(Nanodynamo_degradation),function(j)
{
	i <- Nanodynamo_degradation[[j]]
	commonGenesTmp <- intersect(names(i),names(INSPEcT_degradation))
	print(length(commonGenesTmp))
	round(cor(INSPEcT_degradation[commonGenesTmp],i[commonGenesTmp],method="s",use="c"),2)
})

# WithP WithOutP NanoID 
#  0.11     0.66   0.97 

sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-conda-linux-gnu (64-bit)
# Running under: Rocky Linux 8.7 (Green Obsidian)
# 
# Matrix products: default
# BLAS/LAPACK: /home/mfurlan/miniconda3/envs/nuclearDynamics/lib/libopenblasp-r0.3.23.so
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
