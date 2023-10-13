### Libraries
require("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

### Custom functions
perGeneAccuracy <- function(trainBamU # Unlabeled training sample reads bam
						   ,trainBamL # Labeled training sample reads bam
						   ,sampleP # Unlabeled reads modification probabilities
						   ,sampleBam # Sample of interest bam
						   ,bedFile # Genes BED file
						   ,trueValue # Expected probability value
						   ,trainReads=NULL # Reads used for training
						   ,title=NULL) # Plot title
{
	# Read to gene annotation
	find_gene_counts <- function(alignments, genes) {
		hits <- findOverlaps(query = alignments, subject = genes)
		hits_unique <- hits[which(isUnique(queryHits(hits)))]
		alignments <- alignments[queryHits(hits_unique)]
		mcols(alignments)$gene_name <- names(genes)[subjectHits(hits_unique)]
		counts <- sort(table(mcols(alignments)$gene_name), decreasing = TRUE)
	}

	sampleBam <- sampleBam[!names(sampleBam)%in%trainReads]
	sampleP <- sampleP[!names(sampleP)%in%trainReads]

	print(paste0("Test reads: ",length(sampleP)))

	if(!identical(names(sampleBam),names(sampleP))){print("Names issue!")}

	## Miss-classified reads
	if(trueValue==0)
	{
		sampleBamErr <- sampleBam[names(which(sampleP > 0.5))]
	}else{
		sampleBamErr <- sampleBam[names(which(sampleP < 0.5))]
	}
	
	## Reads annotation
	trainBamUCounts <- find_gene_counts(trainBamU, bedFile)
	trainBamLCounts <- find_gene_counts(trainBamL, bedFile)
	
	sampleBamCounts <- find_gene_counts(sampleBam, bedFile)
	sampleBamErrCounts <- find_gene_counts(sampleBamErr, bedFile)

	## Error frequency
	sampleErrFreq <- sampleBamErrCounts/sampleBamCounts[names(sampleBamErrCounts)]
	
	sampleErrFreqAll <- rep(0,length(sampleBamCounts))
	names(sampleErrFreqAll) <- names(sampleBamCounts)
	sampleErrFreqAll[names(sampleErrFreq)] <- sampleErrFreq

	## Final dataframe
	allGenes <- union(names(trainBamUCounts),names(trainBamLCounts))
	trainBamUCountsAll <- trainBamLCountsAll <- rep(0,length(allGenes))
	names(trainBamUCountsAll) <- names(trainBamLCountsAll) <- allGenes

	trainBamUCountsAll[names(trainBamUCounts)] <- trainBamUCounts
	trainBamLCountsAll[names(trainBamLCounts)] <- trainBamLCounts

	print(paste0("Genes with classification "))
	commonGenes <- intersect(names(trainBamLCountsAll),names(sampleErrFreqAll))

	dataFrame <- data.frame("gene"=commonGenes
						   ,"acc"=(1-as.numeric(sampleErrFreqAll[commonGenes]))
						   ,"counts0"=log10(as.integer(trainBamUCountsAll[commonGenes]))
						   ,"counts24"=log10(as.integer(trainBamLCountsAll[commonGenes])))

    if(!is.null(title))
    {
    		plotOut <- ggplot(dataFrame,aes(x=counts0,y=counts24))+geom_point(aes(color=acc))+scale_color_gradient(low="blue",high="red")+ geom_abline(intercept = 0)

			ggsave(paste0(title,".pdf"),plotOut)
    }
    outTmp <- list(dataFrame=dataFrame, sampleGenes=length(sampleBamCounts), allGenes=length(allGenes), commonGenes=length(commonGenes))
}

### Annotation
bedFile <- genes(txdb)
seqlevelsStyle(bedFile) <- "ENSEMBL"

### Reads alignments
## Unlabeled
load("/path/to/unlabeled/reads/bam.RData")
bamUnlabeled <- bam
rm(bam)

## Labeled
load("/path/to/labeled/reads/bam.RData")
bamLabeled <- bam
rm(bam)

### Training reads
trainReads <- readRDS("/path/to/train.reads.rds")
trainReads <- c(trainReads$train_reads_unlabeled,trainReads$train_reads_labeled)

### Unlabeled
## Classifications
load("/path/to/nanoID.modification.probabilities.unlabeled.RData")
modificationProbabilitiesUnlabeled <- probabilities4sU0h
rm(probabilities4sU0h)

## Accuracy
unlabeled_reads_results <- perGeneAccuracy(trainBamU=bamUnlabeled # Unlabeled training sample reads bam
										  ,trainBamL=bamLabeled # Labeled training sample reads bam
										  ,sampleP=modificationProbabilitiesUnlabeled # Unlabeled reads modification probabilities
										  ,sampleBam=bamUnlabeled # Sample of interest bam
										  ,bedFile=bedFile # Genes BED file
										  ,trueValue=0 # Expected probability value
										  ,trainReads=trainReads # Reads used for training
										  ,title=NULL) # Plot title

### Fully-labeled
## Classification
load("/path/to/nanoID.modification.probabilities.labeled.RData")
modificationProbabilitiesLabeled <- probabilities4sU24h
rm(probabilities4sU24h)

## Accuracy
labeled_reads_results <- perGeneAccuracy(trainBamU=bamUnlabeled # Unlabeled training sample reads bam
										,trainBamL=bamLabeled # Labeled training sample reads bam
										,sampleP=modificationProbabilitiesLabeled # Unlabeled reads modification probabilities
										,sampleBam=bamLabeled # Sample of interest bam
										,bedFile=bedFile # Genes BED file
										,trueValue=1 # Expected probability value
										,trainReads=trainReads # Reads used for training
										,title=NULL) # Plot title

### Plot - Figure S1
pdf("geneLevelAccuracy.pdf",width=4,height=4)
boxplot(list("Unlabeled"=unlabeled_reads_results[[1]][,"acc"]
			,"Labeled"=labeled_reads_results[[1]][,"acc"])
	   ,ylab="Accuracy",main="Gene level",ylim=c(0,1),las=1)
dev.off()

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
