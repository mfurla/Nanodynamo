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
load("/path/to/Files/nanoID.modification.probabilities.unlabeled.RData")
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
load("/path/to/Files/nanoID.modification.probabilities.labeled.RData")
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

### Plot - Figure S2
boxplot(list("Unlabeled"=unlabeled_reads_results[[1]][,"acc"]
			,"Labeled"=labeled_reads_results[[1]][,"acc"])
	   ,ylab="Accuracy",main="Gene level",ylim=c(0,1),las=1)