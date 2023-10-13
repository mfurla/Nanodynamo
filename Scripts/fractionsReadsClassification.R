##############################################################
### Code for the comparison of nascent RNA among fractions ###
##############################################################

#### Untreated
expressionLevels <- readRDS("/path/to/Results/Untreated/expressionLevelsMerged.rds")

### Nascent RNA
## First replicate
chpn1 <- expressionLevels$chromatin[[1]][,"chpn0.33"]/expressionLevels$chromatin[[1]][,"chpt"]
names(chpn1) <- rownames(expressionLevels$chromatin[[1]])
chmn1 <- expressionLevels$chromatin[[1]][,"chmn0.33"]/expressionLevels$chromatin[[1]][,"chmt"]
names(chmn1) <- rownames(expressionLevels$chromatin[[1]])

npn1 <- expressionLevels$nucleoplasmic[[1]][,"npn0.33"]/expressionLevels$nucleoplasmic[[1]][,"npt"]
names(npn1) <- rownames(expressionLevels$nucleoplasmic[[1]])
nmn1 <- expressionLevels$nucleoplasmic[[1]][,"nmn0.33"]/expressionLevels$nucleoplasmic[[1]][,"nmt"]
names(nmn1) <- rownames(expressionLevels$nucleoplasmic[[1]])

cyn1 <- expressionLevels$cytoplasmic[[1]][,"cyn0.33"]/expressionLevels$cytoplasmic[[1]][,"cyt"]
names(cyn1) <- rownames(expressionLevels$cytoplasmic[[1]])
pyn1 <- expressionLevels$polysomal[[1]][,"pyn0.33"]/expressionLevels$polysomal[[1]][,"pyt"]
names(pyn1) <- rownames(expressionLevels$polysomal[[1]])

### Premature RNA
ch1 <- expressionLevels$chromatin[[1]][,"chpt"]/(expressionLevels$chromatin[[1]][,"chpt"]+expressionLevels$chromatin[[1]][,"chmt"])
names(ch1) <- rownames(expressionLevels$chromatin[[1]])
n1 <- expressionLevels$nucleoplasmic[[1]][,"npt"]/(expressionLevels$nucleoplasmic[[1]][,"npt"]+expressionLevels$nucleoplasmic[[1]][,"nmt"])
names(n1) <- rownames(expressionLevels$nucleoplasmic[[1]])
cy1 <- expressionLevels$cytoplasmicP[[1]][,"cypt"]/(expressionLevels$cytoplasmicP[[1]][,"cypt"]+expressionLevels$cytoplasmicP[[1]][,"cymt"])
names(cy1) <- rownames(expressionLevels$cytoplasmicP[[1]])
py1 <- expressionLevels$polysomalP[[1]][,"pypt"]/(expressionLevels$polysomalP[[1]][,"pypt"]+expressionLevels$polysomalP[[1]][,"pymt"])
names(py1) <- rownames(expressionLevels$polysomalP[[1]])

### Significance
wilcox.test(cy1,py1)$p.value
 # [1] 2.403187e-38

### Plots - Figure S10
par(mfrow=c(1,2))
boxplot(list(Chp=chpn1[expressionLevels$expressedGenesYesChpNpP[[1]]],Np=npn1[expressionLevels$expressedGenesYesChpNpP[[1]]],Chm=chmn1[expressionLevels$expressedGenesYesChpNpP[[1]]],Nm=nmn1[expressionLevels$expressedGenesYesChpNpP[[1]]],C=cyn1[expressionLevels$expressedGenesYesChpNpP[[1]]],P=pyn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),outline=FALSE,ylab="Nascent RNA [%]",ylim=c(0,1),varwidth=TRUE,las=2)
boxplot(list(Ch=ch1[expressionLevels$expressedGenesYesChpNpP[[1]]],N=n1[expressionLevels$expressedGenesYesChpNpP[[1]]],C=cy1[expressionLevels$expressedGenesYesChpNpP[[1]]],P=py1[expressionLevels$expressedGenesYesChpNpP[[1]]]),outline=FALSE,ylab="Premature RNA [%]",ylim=c(0,0.6),varwidth=TRUE,las=2)

### Numerical values
## Premature
round(median(chpn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.39
round(median(chmn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.32
round(median(npn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.33
round(median(nmn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.24
round(median(cyn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.25
round(median(pyn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.29

round(median(ch1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.15
round(median(n1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.07
round(median(cy1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.06
round(median(py1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.07

#### Pladienolide B
expressionLevels <- readRDS("/path/to/Results/PladienolideB/expressionLevelsMerged.rds")

### Nascent RNA
## First replicate
chpn1 <- expressionLevels$chromatin[[1]][,"chpn0.33"]/expressionLevels$chromatin[[1]][,"chpt"]
names(chpn1) <- rownames(expressionLevels$chromatin[[1]])
chmn1 <- expressionLevels$chromatin[[1]][,"chmn0.33"]/expressionLevels$chromatin[[1]][,"chmt"]
names(chmn1) <- rownames(expressionLevels$chromatin[[1]])

npn1 <- expressionLevels$nucleoplasmic[[1]][,"npn0.33"]/expressionLevels$nucleoplasmic[[1]][,"npt"]
names(npn1) <- rownames(expressionLevels$nucleoplasmic[[1]])
nmn1 <- expressionLevels$nucleoplasmic[[1]][,"nmn0.33"]/expressionLevels$nucleoplasmic[[1]][,"nmt"]
names(nmn1) <- rownames(expressionLevels$nucleoplasmic[[1]])

cyn1 <- expressionLevels$cytoplasmic[[1]][,"cyn0.33"]/expressionLevels$cytoplasmic[[1]][,"cyt"]
names(cyn1) <- rownames(expressionLevels$cytoplasmic[[1]])
pyn1 <- expressionLevels$polysomal[[1]][,"pyn0.33"]/expressionLevels$polysomal[[1]][,"pyt"]
names(pyn1) <- rownames(expressionLevels$polysomal[[1]])

### Premature RNA
ch1 <- expressionLevels$chromatin[[1]][,"chpt"]/(expressionLevels$chromatin[[1]][,"chpt"]+expressionLevels$chromatin[[1]][,"chmt"])
names(ch1) <- rownames(expressionLevels$chromatin[[1]])
n1 <- expressionLevels$nucleoplasmic[[1]][,"npt"]/(expressionLevels$nucleoplasmic[[1]][,"npt"]+expressionLevels$nucleoplasmic[[1]][,"nmt"])
names(n1) <- rownames(expressionLevels$nucleoplasmic[[1]])
cy1 <- expressionLevels$cytoplasmicP[[1]][,"cypt"]/(expressionLevels$cytoplasmicP[[1]][,"cypt"]+expressionLevels$cytoplasmicP[[1]][,"cymt"])
names(cy1) <- rownames(expressionLevels$cytoplasmicP[[1]])
py1 <- expressionLevels$polysomalP[[1]][,"pypt"]/(expressionLevels$polysomalP[[1]][,"pypt"]+expressionLevels$polysomalP[[1]][,"pymt"])
names(py1) <- rownames(expressionLevels$polysomalP[[1]])

### Plots - Figure S16
pdf("fractionsClassificationPlaBMerged_fullModelGenes.pdf",width=5,height=4)
par(mfrow=c(1,2))
boxplot(list(Chp=chpn1[expressionLevels$expressedGenesYesChpNpP[[1]]],Np=npn1[expressionLevels$expressedGenesYesChpNpP[[1]]],Chm=chmn1[expressionLevels$expressedGenesYesChpNpP[[1]]],Nm=nmn1[expressionLevels$expressedGenesYesChpNpP[[1]]],C=cyn1[expressionLevels$expressedGenesYesChpNpP[[1]]],P=pyn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),outline=FALSE,ylab="Nascent RNA [%]",ylim=c(0,1),varwidth=TRUE,las=2)
boxplot(list(Ch=ch1[expressionLevels$expressedGenesYesChpNpP[[1]]],N=n1[expressionLevels$expressedGenesYesChpNpP[[1]]],C=cy1[expressionLevels$expressedGenesYesChpNpP[[1]]],P=py1[expressionLevels$expressedGenesYesChpNpP[[1]]]),outline=FALSE,ylab="Premature RNA [%]",ylim=c(0,0.6),varwidth=TRUE,las=2)
dev.off()

### Numerical values
## Premature
round(median(chpn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.22
round(median(chmn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.16
round(median(npn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.25
round(median(nmn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.16
round(median(cyn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.21
round(median(pyn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.33

round(median(ch1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.25
round(median(n1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.10
round(median(cy1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.09
round(median(py1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.10

#### PlaB vs Untreated - Figure 3B
## Chp
expressionLevelsUntreated <- readRDS("/path/to/Results/Untreated/expressionLevelsMerged.rds")
expressionLevelsPlaB <- readRDS("/path/to/Results/PladienolideB/expressionLevelsMerged.rds")

chUntreated <- expressionLevelsUntreated$chromatin[[1]][,"chpt"]/(expressionLevelsUntreated$chromatin[[1]][,"chpt"]+expressionLevelsUntreated$chromatin[[1]][,"chmt"])
names(chUntreated) <- rownames(expressionLevelsUntreated$chromatin[[1]])

chPlaB <- expressionLevelsPlaB$chromatin[[1]][,"chpt"]/(expressionLevelsPlaB$chromatin[[1]][,"chpt"]+expressionLevelsPlaB$chromatin[[1]][,"chmt"])
names(chPlaB) <- rownames(expressionLevelsPlaB$chromatin[[1]])

ks.test(chUntreated[expressionLevelsUntreated$expressedGenesYesChpNpP[[1]]],chPlaB[expressionLevelsPlaB$expressedGenesYesChpNpP[[1]]])$p.value # 0
wilcox.test(chUntreated[expressionLevelsUntreated$expressedGenesYesChpNpP[[1]]],chPlaB[expressionLevelsPlaB$expressedGenesYesChpNpP[[1]]],alternative="less")$p.value # 5.467804e-68

allGenes <- intersect(names(chUntreated),names(chPlaB))
allGenesFull <- intersect(expressionLevelsUntreated$expressedGenesYesChpNpP[[1]],expressionLevelsPlaB$expressedGenesYesChpNpP[[1]])

smoothScatter(log2(chUntreated[allGenes]),log2(chPlaB[allGenes]),xlab="Untreated",ylab="Pladienolide B",main="Log2 Proportion of\nPremature RNA reads",pch=".")
abline(0,1,lwd=2,col=1);abline(1,1,lwd=2,col=2);abline(2,1,lwd=2,col=3)

## Nuclear nascent
chpnUntreated<- expressionLevelsUntreated$chromatin[[1]][,"chpn0.33"]/expressionLevelsUntreated$chromatin[[1]][,"chpt"]
names(chpnUntreated) <- rownames(expressionLevelsUntreated$chromatin[[1]])
chmnUntreated<- expressionLevelsUntreated$chromatin[[1]][,"chmn0.33"]/expressionLevelsUntreated$chromatin[[1]][,"chmt"]
names(chmnUntreated) <- rownames(expressionLevelsUntreated$chromatin[[1]])
npnUntreated<- expressionLevelsUntreated$nucleoplasmic[[1]][,"npn0.33"]/expressionLevelsUntreated$nucleoplasmic[[1]][,"npt"]
names(npnUntreated) <- rownames(expressionLevelsUntreated$nucleoplasmic[[1]])
nmnUntreated<- expressionLevelsUntreated$nucleoplasmic[[1]][,"nmn0.33"]/expressionLevelsUntreated$nucleoplasmic[[1]][,"nmt"]
names(nmnUntreated) <- rownames(expressionLevelsUntreated$nucleoplasmic[[1]])

chpnPlaB<- expressionLevelsPlaB$chromatin[[1]][,"chpn0.33"]/expressionLevelsPlaB$chromatin[[1]][,"chpt"]
names(chpnPlaB) <- rownames(expressionLevelsPlaB$chromatin[[1]])
chmnPlaB<- expressionLevelsPlaB$chromatin[[1]][,"chmn0.33"]/expressionLevelsPlaB$chromatin[[1]][,"chmt"]
names(chmnPlaB) <- rownames(expressionLevelsPlaB$chromatin[[1]])
npnPlaB<- expressionLevelsPlaB$nucleoplasmic[[1]][,"npn0.33"]/expressionLevelsPlaB$nucleoplasmic[[1]][,"npt"]
names(npnPlaB) <- rownames(expressionLevelsPlaB$nucleoplasmic[[1]])
nmnPlaB<- expressionLevelsPlaB$nucleoplasmic[[1]][,"nmn0.33"]/expressionLevelsPlaB$nucleoplasmic[[1]][,"nmt"]
names(nmnPlaB) <- rownames(expressionLevelsPlaB$nucleoplasmic[[1]])

ks.test(chpnUntreated[expressionLevelsUntreated$expressedGenesYesChpNpP[[1]]],chpnPlaB[expressionLevelsPlaB$expressedGenesYesChpNpP[[1]]])$p.value # 0
ks.test(chmnUntreated[expressionLevelsUntreated$expressedGenesYesChpNpP[[1]]],chmnPlaB[expressionLevelsPlaB$expressedGenesYesChpNpP[[1]]])$p.value # 0
ks.test(npnUntreated[expressionLevelsUntreated$expressedGenesYesChpNpP[[1]]],npnPlaB[expressionLevelsPlaB$expressedGenesYesChpNpP[[1]]])$p.value # 0
ks.test(nmnUntreated[expressionLevelsUntreated$expressedGenesYesChpNpP[[1]]],nmnPlaB[expressionLevelsPlaB$expressedGenesYesChpNpP[[1]]])$p.value # 0

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