### Libraries
library(ggplot2)

### Data loading - from the supplemental material of de Pretis et al. Genome Research 2017
load("path/to/3T9MycER_rpkms.RData")

### CV profiling
## Total RNA
totalRNAMeanTmp <- apply(totalRpkmsExons,1,mean)
totalRNASdTmp <- apply(totalRpkmsExons,1,sd)
totalRNACvTmp <- totalRNASdTmp/totalRNAMeanTmp

## Premature RNA
prematureRNAMeanTmp <- apply(totalRpkmsIntrons,1,mean)
prematureRNASdTmp <- apply(totalRpkmsIntrons,1,sd)
prematureRNACvTmp <- prematureRNASdTmp/prematureRNAMeanTmp

## Total Nascent RNA
totalNascentRNAMeanTmp <- apply(labeledRpkmsExons,1,mean)
totalNascentRNASdTmp <- apply(labeledRpkmsExons,1,sd)
totalNascentRNACvTmp <- totalNascentRNASdTmp/totalNascentRNAMeanTmp

## Premature Nascent RNA
prematureNascentRNAMeanTmp <- apply(labeledRpkmsIntrons,1,mean)
prematureNascentRNASdTmp <- apply(labeledRpkmsIntrons,1,sd)
prematureNascentRNACvTmp <- prematureNascentRNASdTmp/prematureNascentRNAMeanTmp

### Plot - Figure S2
pdf("CV_distribution.pdf",width=4,height=4)
ggplot(data.frame(CV=c(totalRNACvTmp,prematureRNACvTmp,totalNascentRNACvTmp,prematureNascentRNACvTmp)
				 ,RNA_species=rep(c("Total_RNA","Premature_RNA","Total_Nascent_RNA","Premature_Nascent_RNA")
				 ,each=length(totalRNACvTmp))), aes(x = CV, fill = RNA_species)) + geom_density(alpha = 0.5) + theme(legend.position = c(0.65, 0.65))
dev.off()

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
