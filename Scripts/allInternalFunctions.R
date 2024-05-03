### Loading libraries
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
library("DESeq2")
library("vioplot")
library("rcartocolor")
library("e1071")
library("DescTools")
library("igraph")
library("RColorBrewer")

### Loop functions
mcsapply <- function( X, FUN, ... ) do.call('cbind', mclapply( X, FUN, ... ))

### Function for smoothScatters color scale
fudgeit <- function()
{
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  fields::image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =F,legend.width=3,legend.lab="density",legend.cex=0.85)
}

### Genes counts profiling
##  - From BAM to counts,
##  - From BAMs to expression data.

genesCountsProfiling <- function(bamPath
                                 ,nascentPath
                                 ,labelingTime
                                 ,txdb
                                 ,label
                                 ,prematureFlag
                                 ,countsTh=0
                                 ,minoverlap_I=minoverlap_I
                                 ,minoverlap_E=minoverlap_E)
{
  # BAM loading
  bamTmp <- get(load(bamPath))
  
  # Definition of regions classified as exonic in at-least one isoform
  exonsDB <- reduce(exonsBy(txdb ,'gene'))
  exonsDB <- exonsDB[elementNROWS(range(exonsDB))==1]
  
  # Definition of gaps as intronic regions
  intronsDB <- psetdiff(unlist(range(exonsDB)),exonsDB)
  intronsDB <- intronsDB[elementNROWS(intronsDB)>0]
  
  # Correct chromosomes names
  seqlevelsStyle(exonsDB) <- "ENSEMBL"
  seqlevelsStyle(intronsDB) <- "ENSEMBL"
  
  # Overlap between reads and exonic regions
  geneOverlaps <- findOverlaps(grglist(bamTmp,drop.D.ranges=TRUE),unlist(range(exonsDB)),minoverlap=minoverlap_E)
  
  # Unique overlaps
  geneOverlaps <- geneOverlaps[isUnique(queryHits(geneOverlaps)),]
  bamTmp <- bamTmp[names(bamTmp)[queryHits(geneOverlaps)]]
  
  # Overlap between reads and intronic regions
  intronicOverlaps <- findOverlaps(grglist(bamTmp,drop.D.ranges=TRUE),intronsDB,minoverlap=minoverlap_I)
  intronicOverlaps <- intronicOverlaps[isUnique(queryHits(intronicOverlaps)),]
  intronicReads <- names(bamTmp)[queryHits(intronicOverlaps)]
  names(intronicReads) <- names(intronsDB)[subjectHits(intronicOverlaps)]
  
  bamTmp <- bamTmp[!(names(bamTmp)%in%intronicReads)]
  
  # Annotation of no-intronic reads
  exonicOverlaps <- findOverlaps(grglist(bamTmp,drop.D.ranges=TRUE),exonsDB,minoverlap=minoverlap_E)
  exonicOverlaps <- exonicOverlaps[isUnique(queryHits(exonicOverlaps)),]
  exonicReads <- names(bamTmp)[queryHits(exonicOverlaps)]
  names(exonicReads) <- names(exonsDB)[subjectHits(exonicOverlaps)]
  
  # PreExisting/Nascent classification probabilities
  if(grepl(".RData$",nascentPath))
  {
    nascentTmp <- get(load(nascentPath))
  }else if(grepl(".rds$",nascentPath))
  {
    nascentTmp <- readRDS(nascentPath)
  }else{
    print("Strange nascent data type!");return(NULL)
  }
  
  # Selection of reads with a complete classification
  commonReadsTmp <- intersect(names(nascentTmp),c(intronicReads,exonicReads))
  if(length(commonReadsTmp)<1){print("No reads with full classification; check input data!");return(NULL)}
  
  # Unique dataframe
  reportTmp <- data.frame("X.Read_ID"=c(intronicReads,exonicReads)
                          ,"Gene_id"=c(names(intronicReads),names(exonicReads))
                          ,"Mature_status"=c(rep("Premature",length(intronicReads)),rep("Mature",length(exonicReads)))
                          ,"Nascent_status"=nascentTmp[c(intronicReads,exonicReads)])
  rownames(reportTmp) <- reportTmp$X.Read_ID
  reportTmp <- reportTmp[commonReadsTmp,]
  reportTmp$Nascent_status <- c("N","P")[as.integer(nascentTmp[commonReadsTmp]<0.5)+1]
  
  if(prematureFlag)
  {
    # If premature RNA is required, classification of reads as Premature/Mature and Nascent/Preexisting
    reportTmp <- as.data.frame(reportTmp %>% group_by(Gene_id) %>% summarize(pp = length(intersect(grep("Premature",Mature_status),grep("P",Nascent_status)))
                                                                             ,pn = length(intersect(grep("Premature",Mature_status),grep("N",Nascent_status)))
                                                                             ,mp = length(intersect(grep("Mature",Mature_status),grep("P",Nascent_status)))
                                                                             ,mn = length(intersect(grep("Mature",Mature_status),grep("N",Nascent_status)))))
    
    rownames(reportTmp) <- reportTmp[,"Gene_id"]
    reportTmp <- reportTmp[,-1]
    
    # If required, expression threshold
    if(is.numeric(countsTh))
    {
      selectionTmp <- apply(reportTmp,1,function(i){all(i>=countsTh)})
      reportTmp <- reportTmp[selectionTmp,]
    }
    
    # Steady state values
    reportTmp <- reportTmp %>% mutate(mn0=0,pn0=0,pp0=pp+pn,mp0=mp+mn)
    reportTmp <- reportTmp %>% mutate(pt=pn0+pp0,mt=mn0+mp0)
  }else{
    # If premature RNA is not required, classification of reads as Nascent/Preexisting
    reportTmp <- as.data.frame(reportTmp %>% group_by(Gene_id) %>%summarize(p = length(grep("P",Nascent_status))
                                                                            ,n = length(grep("N",Nascent_status))))
    
    rownames(reportTmp) <- reportTmp[,"Gene_id"]
    reportTmp <- reportTmp[,-1]
    
    # If required, expression threshold
    if(is.numeric(countsTh))
    {
      selectionTmp <- apply(reportTmp,1,function(i){all(i>=countsTh)})
      reportTmp <- reportTmp[selectionTmp,]
    }
    
    # Steady state values
    reportTmp <- reportTmp %>% mutate(n0=0,p0=p+n)
    reportTmp <- reportTmp %>% mutate(t=n0+p0)
  }
  
  # Colnames definition
  colnamesTmp <- colnames(reportTmp)
  colnamesTmp[!(grepl("Gene_id",colnamesTmp)|grepl("0$",colnamesTmp)|grepl("t$",colnamesTmp))] <- paste0(label,colnamesTmp[!(grepl("Gene_id",colnamesTmp)|grepl("0$",colnamesTmp)|grepl("t$",colnamesTmp))],labelingTime)
  
  colnamesTmp[(grepl("0$",colnamesTmp)|grepl("t$",colnamesTmp))] <- paste0(label,colnamesTmp[(grepl("0$",colnamesTmp)|grepl("t$",colnamesTmp))])
  colnames(reportTmp) <- colnamesTmp
  
  return(reportTmp)
}

expressionDataEstimation <- function(bamPaths # Paths to BAM files.
                                     ,nascentPaths # Paths to nascent RNA probabilities.
                                     ,tbl # Matrix of sequencing statistics.
                                     ,labelingTime # Labeling time of the fractionation samples.
                                     ,labelingTimePoly=NULL # Labeling time of the polysomal samples.
                                     ,txdb # Annotation object.
                                     ,minoverlap_I=10 # Minimum overlap for intronic reads.
                                     ,minoverlap_E=10 # Minimum overlap for exonic reads.
                                     ,cpus=1 # Number of cpus.
                                     ,saveDataDistributions=FALSE # If TRUE, save gene expression levels data distributions.
)
{
  # Is polysomal RNA part of the dataset?
  polysomalFlag <- any(grepl("Poly",names(bamPaths))&!is.null(labelingTimePoly))
  if(polysomalFlag){print("Polysomal RNA configuration.")}else{print("No Polysomal RNA configuration.")}
  
  # Fractions and replicates
  fractionsNames <- unique(sapply(strsplit(names(bamPaths),"_"),"[[",1))
  replicatesNames <- unique(sapply(strsplit(names(bamPaths),"_"),"[[",2))
  
  # For all the cellular fractions, we annotate the reads as premature, mature, nascent, and pre-existing.
  # This function should be extended to process samples with variable labeling time. For now, we have only
  # fractions vs polysomal labeling.
  
  # Identification of premature/mature and nascent/pre-existing reads for all the samples
  readsAnnotation <- mclapply(names(bamPaths),function(i)
  {
    print(i)
    
    bamPathTmp <- bamPaths[[i]]
    nascentPathTmp <- nascentPaths[[i]]
    
    suppressWarnings(genesCountsProfiling(bamPath=bamPathTmp
                                          ,nascentPath=nascentPathTmp
                                          ,label=switch(strsplit(i,"_")[[1]][[1]],Chr = "ch",Nuc = "n",Cyt = "cy",Poly = "py")
                                          ,labelingTime=switch(strsplit(i,"_")[[1]][[1]],Poly = labelingTimePoly,labelingTime)
                                          ,prematureFlag=TRUE
                                          ,countsTh=0
                                          ,txdb=txdb
                                          ,minoverlap_I=minoverlap_I
                                          ,minoverlap_E=minoverlap_E))
    
  },mc.cores=cpus)
  # saveRDS(readsAnnotation,"readsAnnotation.rds")
  
  names(readsAnnotation) <- names(bamPaths)
  readsAnnotationBackup <- readsAnnotation
  
  # Fractions of premature/mature and nascent/pre-existing reads for all the samples
  yieldScalingFactors <- sapply(readsAnnotation,colSums)
  yieldScalingFactors <- apply(yieldScalingFactors[grepl(labelingTime,rownames(yieldScalingFactors)),],2,function(i)i/sum(i))
  rownames(yieldScalingFactors) <- gsub("ch","",rownames(yieldScalingFactors))
  
  yieldScalingFactors <- rbind(yieldScalingFactors
                               ,"pp0"=colSums(yieldScalingFactors[grepl("^p",rownames(yieldScalingFactors)),])
                               ,"pn0"=rep(0,ncol(yieldScalingFactors))
                               ,"mp0"=colSums(yieldScalingFactors[grepl("^m",rownames(yieldScalingFactors)),])
                               ,"mn0"=rep(0,ncol(yieldScalingFactors))
                               ,"pt"=colSums(yieldScalingFactors[grepl("^p",rownames(yieldScalingFactors)),])
                               ,"mt"=colSums(yieldScalingFactors[grepl("^m",rownames(yieldScalingFactors)),]))
  
  # From reads fractions to yield scaling factors
  normFactor <- (tbl[,"PolyA"]/tbl[,"cells"])
  normFactor <- sapply(fractionsNames,function(i)median(normFactor[grepl(i,names(normFactor))]))
  normFactor <- sapply(colnames(yieldScalingFactors),function(i)
  {
    yieldScalingFactors[,i]*normFactor[[strsplit(i,"_")[[1]][[1]]]]
  })
  meanNormFactor <- sapply(fractionsNames,function(i)apply(normFactor[,grepl(i,colnames(normFactor))],1,mean))
  colnames(meanNormFactor) <- paste0(colnames(meanNormFactor),"_Mean")
  normFactor <- cbind(normFactor,meanNormFactor)
  
  # Counts re-framing - missing species are set to zero
  readsAnnotation <- split(readsAnnotation,sapply(strsplit(names(readsAnnotation),"_"),"[[",1))
  readsAnnotation <- lapply(readsAnnotation,function(i)
  {
    conditionsTmp <- names(i)
    allGenes <- sort(unique(unlist(sapply(i,rownames))))
    i <- lapply(i,function(j)
    {
      tmp <- matrix(0,nrow=length(allGenes),ncol=ncol(j))
      colnames(tmp) <- colnames(j)
      rownames(tmp) <- allGenes
      tmp[rownames(j),] <- as.matrix(j)
      tmp
    })
    i <- do.call("cbind",i)
    colnames(i) <- paste0(sapply(conditionsTmp,function(l)rep(l,ncol(i)/2)),"-",colnames(i))
    i
  })
  
  # Independently for each fraction, we run DESeq2 to estimate gene counts distributions.
  DESeq2Results <- lapply(readsAnnotation,function(i)
  {
    fractionTmp <- strsplit(colnames(i)[[1]],"_")[[1]][[1]]
    
    outTmp <- lapply(unique(sapply(strsplit(colnames(i),"[-]"),"[[",2)),function(j)
    {
      l <- i[,grep(paste0(j,"$"),colnames(i),value=T)]
      if(!(all(l==0)))
      {
        coldataTmp <- data.frame(condition=sapply(strsplit(colnames(l),"[-]"),"[[",2))
        rownames(coldataTmp) <- colnames(l)
        
        ddsTmp <- DESeqDataSetFromMatrix(countData = l,
                                         colData = coldataTmp,
                                         design = ~ 1)
        ddsTmp <- estimateSizeFactors(ddsTmp)
        
        # sizeFactors(ddsTmp) <- colSums(l)/1e6
        
        # j2 <- gsub(switch(fractionTmp,"Chr"="^ch","Nuc"="^n","Cyt"="^cy","Poly"="^py"),"",j)
        # print(j2)
        # print(sapply(strsplit(colnames(l),"-"),"[[",1))
        # sizeFactors(ddsTmp) <- colSums(l)/(1e6*normFactor[j2,sapply(strsplit(colnames(l),"-"),"[[",1)])
        
        ddsTmp <- estimateDispersions(ddsTmp)
        
        ddsTmp
      }else{NULL}
    })
    names(outTmp) <- sapply(strsplit(colnames(i),"-"),"[[",2)[1:length(outTmp)]
    outTmp
  })
  
  # Independently for each fraction, we retrieve RNA species' normalized counts.
  DESeq2NormalizedCounts <- lapply(DESeq2Results,function(i)
  {
    normalizedCountsTmp <- sapply(i,function(j)
    {
      if(!is.null(j))counts(j,norm=T)
    })
    
    for(k in names(normalizedCountsTmp))
    {
      if(is.null(normalizedCountsTmp[[k]]))
      {
        normalizedCountsTmp[[k]] <- matrix(0,nrow=nrow(normalizedCountsTmp[[1]])
                                           ,ncol=ncol(normalizedCountsTmp[[1]]))
        colnames(normalizedCountsTmp[[k]]) <- paste0(sapply(strsplit(colnames(normalizedCountsTmp[[1]]),"[-]"),"[[",1),"-",k)
      }
    }
    
    do.call("cbind",normalizedCountsTmp)
  })
  
  # Independently for each fraction, we retrieve RNA species' dispersions.
  DESeq2Dispersions <- lapply(DESeq2Results,function(i)
  {
    dispersionsTmp <- sapply(i,function(j)
    {
      if(!is.null(j))dispersions(j)
    })
    
    dispersionsTmp[sapply(dispersionsTmp,is.null)] <- NA
    
    do.call("cbind",dispersionsTmp)
  })
  
  for(i in names(DESeq2NormalizedCounts))rownames(DESeq2Dispersions[[i]]) <- rownames(DESeq2NormalizedCounts[[i]])
  
  # For each fraction, we retrieve RNA species' mean expression levels.
  DESeq2MeanNormalizedCounts <- lapply(DESeq2Results,function(i){sapply(i,function(j){if(!is.null(j)){outTmp <- mcols(j)[,"baseMean"];names(outTmp)<-rownames(mcols(j));outTmp}})})
  DESeq2MeanNormalizedCounts <- lapply(DESeq2MeanNormalizedCounts,function(i){sapply(i,function(j)if(is.null(j)){rep(0,length(i[[1]]))}else{j})})
  
  # DESeq2MeanDispersions <- lapply(DESeq2Results,function(i){sapply(i,function(j){if(!is.null(j)){outTmp <- mcols(j)[,"dispMAP"];names(outTmp)<-rownames(mcols(j));outTmp}})})
  # DESeq2MeanDispersions <- lapply(DESeq2MeanDispersions,function(i){sapply(i,function(j)if(is.null(j)){rep(0,length(i[[1]]))}else{j})})
  
  # Normalized counts for replicates + mean value
  DESeq2NormalizedCounts <- lapply(fractionsNames,function(i)
  {
    repsNormalizedCounts <- DESeq2NormalizedCounts[[i]]
    meanNormalizedCounts <- DESeq2MeanNormalizedCounts[[i]]
    
    replicatesTmp <- unique(sapply(strsplit(colnames(repsNormalizedCounts),"[-]"),"[[",1))
    
    repsNormalizedCounts <- lapply(replicatesTmp,function(j)
    {
      outTmp <- repsNormalizedCounts[,grepl(j,colnames(repsNormalizedCounts))]
      colnames(outTmp) <- gsub(paste0(j,"[-]"),"",colnames(outTmp))
      outTmp
    })
    repsNormalizedCounts <- append(list(meanNormalizedCounts),repsNormalizedCounts)
    names(repsNormalizedCounts) <- c(paste0(i,"_Mean"),replicatesTmp)
    repsNormalizedCounts
  })
  names(DESeq2NormalizedCounts) <- fractionsNames
  # DESeq2Dispersions=list("Mean"=DESeq2MeanDispersions,"Replicates"=DESeq2Dispersions)
  
  # One expression matrix for each sample + mean
  DESeq2NormalizedCounts <- lapply(names(DESeq2NormalizedCounts[[1]]),function(i)
  {
    i <- sapply(strsplit(i,"_"),"[[",2)
    j <- sapply(DESeq2NormalizedCounts,function(m)unname(m[grepl(i,names(m))]))
    allGenes <- sort(unique(unlist(sapply(j,rownames))))
    k <- lapply(j,function(l)
    {
      tmp <- matrix(0,nrow=length(allGenes),ncol=ncol(l))
      colnames(tmp) <- colnames(l)
      rownames(tmp) <- allGenes
      tmp[rownames(l),] <- as.matrix(l)
      tmp
    })
    k <- do.call("cbind",k)
    k
  })
  names(DESeq2NormalizedCounts) <- c("Mean",replicatesNames)
  
  # DESeq2Dispersions <- lapply(DESeq2Dispersions,function(i)
  # {
  #   allGenes <- sort(unique(unlist(sapply(i,rownames))))
  #   k <- lapply(i,function(l)
  #   {
  #     tmp <- matrix(NA,nrow=length(allGenes),ncol=ncol(l))
  #     colnames(tmp) <- colnames(l)
  #     rownames(tmp) <- allGenes
  #     tmp[rownames(l),] <- as.matrix(l)
  #     tmp
  #   })
  #   k <- do.call("cbind",k)
  #   k
  # })
  
  # One dispersion matrix for each sample + mean
  allGenes <- sort(unique(unlist(sapply(DESeq2Dispersions,rownames))))
  DESeq2Dispersions <- lapply(DESeq2Dispersions,function(l)
  {
    tmp <- matrix(NA,nrow=length(allGenes),ncol=ncol(l))
    colnames(tmp) <- colnames(l)
    rownames(tmp) <- allGenes
    tmp[rownames(l),] <- as.matrix(l)
    tmp
  })
  DESeq2Dispersions <- do.call("cbind",DESeq2Dispersions)
  
  # Definition of normalized expression levels + Confidence Intervals
  normalizedCountsStatistics <- mclapply(names(DESeq2NormalizedCounts),function(i)
  {
    print(i)
    muTmp <- DESeq2NormalizedCounts[[i]]
    # if(i=="Mean"){alphaTmp <- DESeq2Dispersions[["Mean"]]}else{alphaTmp <- DESeq2Dispersions[["Replicates"]]}
    alphaTmp <- DESeq2Dispersions
    
    bottomCITmp <- upCITmp <- medianTmp <- meanTmp <- matrix(0,nrow=nrow(muTmp),ncol=ncol(muTmp))
    colnames(bottomCITmp) <- colnames(upCITmp) <- colnames(medianTmp) <- colnames(meanTmp) <- colnames(muTmp)
    rownames(bottomCITmp) <- rownames(upCITmp) <- rownames(medianTmp) <- rownames(meanTmp) <- rownames(muTmp)
    
    if(saveDataDistributions){dataTmp <- list()}else{dataTmp=NULL}
    
    for(j in colnames(muTmp))
    {
      # print(j)
      if(any(muTmp[,j]!=0))
      {
        fractionTmp <- paste0(c("Chr","Nuc","Cyt","Poly")[sapply(c("^ch","^n","^cy","^py"),function(l)grepl(l,j))],"_",i)
        speciesTmp <- sapply(rownames(normFactor),function(l)grepl(paste0(l,"$"),j))
        set.seed(1)
        simDataTmp <- t(mcsapply(names(muTmp[,j]),function(l)rnbinom(1000,mu=muTmp[l,j],size=1/alphaTmp[l,j]),mc.cores=max(1,floor((cpus-length(DESeq2NormalizedCounts))/length(DESeq2NormalizedCounts)))))
        simDataTmp[!is.finite(simDataTmp)] <- 0
        rownames(simDataTmp) <- names(muTmp[,j])
        simDataTmp <- apply(simDataTmp,2,function(l)(1e6*normFactor[speciesTmp,fractionTmp]*(l/sum(l))))
        simDataCI <- t(apply(simDataTmp,1,function(j)quantile(j,c(0.025,0.5,0.975))))
        bottomCITmp[names(muTmp[,j]),j] <- simDataCI[,1]
        medianTmp[names(muTmp[,j]),j] <- simDataCI[,2]
        upCITmp[names(muTmp[,j]),j] <- simDataCI[,3]
        meanTmp[names(muTmp[,j]),j] <- 1e6*normFactor[speciesTmp,fractionTmp]*(muTmp[,j]/sum(muTmp[,j]))
        
        if(saveDataDistributions){dataTmp[[j]] <- simDataTmp}
      }
    }
    
    list(bottom=bottomCITmp,up=upCITmp,median=medianTmp,mean=meanTmp,distributions=dataTmp)
  },mc.cores=min(cpus,length(DESeq2NormalizedCounts)))
  names(normalizedCountsStatistics) <- names(DESeq2NormalizedCounts)
  
  # We remove the RNA species discarded during modeling (i.e. Cytoplasmic and Polysomal premature RNA)
  normalizedCountsStatisticsSmall <- lapply(normalizedCountsStatistics,function(i)
  {
    lapply(i[1:4],function(j)
    {
      j <- j[,!(grepl("^cyp",colnames(j))|grepl("^pyp",colnames(j)))]
      colnames(j) <- gsub("^cym","cy",colnames(j))
      colnames(j) <- gsub("^pym","py",colnames(j))
      j[j<1e-10] <- 1e-10
      j
    })
  })
  
  splitGenes <- function(normalizedCountsTmp)
  {
    
    j <- lapply(normalizedCountsTmp,"[[",4)
    
    expressedGenesYesChpn <- apply(sapply(j,function(i)i[,paste0("chpn",labelingTime)]>1e-10),1,any)
    expressedGenesYesChpp <- apply(sapply(j,function(i)i[,paste0("chpp",labelingTime)]>1e-10),1,any)
    
    expressedGenesYesChmn <- apply(sapply(j,function(i)i[,paste0("chmn",labelingTime)]>1e-10),1,any)
    expressedGenesYesChmp <- apply(sapply(j,function(i)i[,paste0("chmp",labelingTime)]>1e-10),1,any)
    
    expressedGenesYesNpn <- apply(sapply(j,function(i)i[,paste0("npn",labelingTime)]>1e-10),1,any)
    expressedGenesYesNpp <- apply(sapply(j,function(i)i[,paste0("npp",labelingTime)]>1e-10),1,any)
    
    expressedGenesYesNmn <- apply(sapply(j,function(i)i[,paste0("nmn",labelingTime)]>1e-10),1,any)
    expressedGenesYesNmp <- apply(sapply(j,function(i)i[,paste0("nmp",labelingTime)]>1e-10),1,any)
    
    expressedGenesYesCn <- apply(sapply(j,function(i)i[,paste0("cyn",labelingTime)]>1e-10),1,any)
    expressedGenesYesCp <- apply(sapply(j,function(i)i[,paste0("cyp",labelingTime)]>1e-10),1,any)
    
    if(polysomalFlag)
    {
      expressedGenesYesPn <- apply(sapply(j,function(i)i[,paste0("pyn",labelingTime)]>1e-10),1,any)
      expressedGenesYesPp <- apply(sapply(j,function(i)i[,paste0("pyp",labelingTime)]>1e-10),1,any)
    }else{
      expressedGenesYesPp <- expressedGenesYesPn <- rep(FALSE,nrow(j[[1]]))
      names(expressedGenesYesPp) <- names(expressedGenesYesPn) <- rownames(j[[1]])
    }
    
    # Expressed genes with chromatin premature RNA, nucleoplasmic premature RNA, and polysomal RNA
    expressedGenesYesChpNpP <- names(which(expressedGenesYesChpp&expressedGenesYesChpn
                                           &expressedGenesYesChmp&expressedGenesYesChmn
                                           &expressedGenesYesNpp&expressedGenesYesNpn
                                           &expressedGenesYesNmp&expressedGenesYesNmn
                                           &expressedGenesYesCp&expressedGenesYesCn
                                           &expressedGenesYesPp&expressedGenesYesPn))
    
    # Expressed genes with chromatin premature RNA, nucleoplasmic premature RNA, and without polysomal RNA
    expressedGenesYesChpNpNoP <- names(which(expressedGenesYesChpp&expressedGenesYesChpn
                                             &expressedGenesYesChmp&expressedGenesYesChmn
                                             &expressedGenesYesNpp&expressedGenesYesNpn
                                             &expressedGenesYesNmp&expressedGenesYesNmn
                                             &expressedGenesYesCp&expressedGenesYesCn
                                             &(!expressedGenesYesPp|!expressedGenesYesPn)))
    
    # Expressed genes with chromatin premature RNA, polysomal RNA, and without nucleoplasmic premature RNA
    expressedGenesYesChpPNoNp <- names(which(expressedGenesYesChpp&expressedGenesYesChpn
                                             &expressedGenesYesChmp&expressedGenesYesChmn
                                             &(!expressedGenesYesNpp|!expressedGenesYesNpn)
                                             &expressedGenesYesNmp&expressedGenesYesNmn
                                             &expressedGenesYesCp&expressedGenesYesCn
                                             &expressedGenesYesPp&expressedGenesYesPn))
    
    # Expressed genes with chromatin premature RNA, and without nucleoplasmic premature RNA and polysomal RNA
    expressedGenesYesChpNoNpP <- names(which(expressedGenesYesChpp&expressedGenesYesChpn
                                             &expressedGenesYesChmp&expressedGenesYesChmn
                                             &(!expressedGenesYesNpp|!expressedGenesYesNpn)
                                             &expressedGenesYesNmp&expressedGenesYesNmn
                                             &expressedGenesYesCp&expressedGenesYesCn
                                             &(!expressedGenesYesPp|!expressedGenesYesPn)))
    
    # Expressed genes with polysomal RNA, and without chromatin premature RNA
    expressedGenesYesPNoChp <- names(which((!expressedGenesYesChpp|!expressedGenesYesChpn)
                                           &expressedGenesYesChmp&expressedGenesYesChmn
                                           &(!expressedGenesYesNpp|!expressedGenesYesNpn)
                                           &expressedGenesYesNmp&expressedGenesYesNmn
                                           &expressedGenesYesCp&expressedGenesYesCn
                                           &expressedGenesYesPp&expressedGenesYesPn))
    
    # Expressed genes without polysomal RNA, and chromatin premature RNA
    expressedGenesNoChpP <- names(which((!expressedGenesYesChpp|!expressedGenesYesChpn)
                                        &expressedGenesYesChmp&expressedGenesYesChmn
                                        &(!expressedGenesYesNpp|!expressedGenesYesNpn)
                                        &expressedGenesYesNmp&expressedGenesYesNmn
                                        &expressedGenesYesCp&expressedGenesYesCn
                                        &(!expressedGenesYesPp|!expressedGenesYesPn)))
    
    normalizedCountsTmpYesChpNpP <- lapply(normalizedCountsTmp,function(i)lapply(i,function(j)j[expressedGenesYesChpNpP,]))
    normalizedCountsTmpYesChpNpNoP <- lapply(normalizedCountsTmp,function(i)lapply(i,function(j)
    {
      j <- j[expressedGenesYesChpNpNoP,]
      j <- j[,!(grepl("^py",colnames(j)))]
      j
    }))
    
    normalizedCountsTmpYesChpPNoNp <- lapply(normalizedCountsTmp,function(i)lapply(i,function(j)
    {
      j <- j[expressedGenesYesChpPNoNp,]
      j <- j[,!(grepl("^np",colnames(j)))]
      colnames(j) <- gsub("^nm","n",colnames(j))
      j
    }))
    
    normalizedCountsTmpYesChpNoNpP <- lapply(normalizedCountsTmp,function(i)lapply(i,function(j)
    {
      j <- j[expressedGenesYesChpNoNpP,]
      j <- j[,!(grepl("^py",colnames(j)))]
      j <- j[,!(grepl("^np",colnames(j)))]
      colnames(j) <- gsub("^nm","n",colnames(j))
      j
    }))
    
    normalizedCountsTmpYesPNoChp <- lapply(normalizedCountsTmp,function(i)lapply(i,function(j)
    {
      j <- j[expressedGenesYesPNoChp,]
      j <- j[,!(grepl("^np",colnames(j)))]
      colnames(j) <- gsub("^nm","n",colnames(j))
      j <- j[,!(grepl("^chp",colnames(j)))]
      colnames(j) <- gsub("^chm","ch",colnames(j))
      j
    }))
    
    normalizedCountsTmpNoChpP <- lapply(normalizedCountsTmp,function(i)lapply(i,function(j)
    {
      j <- j[expressedGenesNoChpP,]
      j <- j[,!(grepl("^py",colnames(j)))]
      j <- j[,!(grepl("^np",colnames(j)))]
      colnames(j) <- gsub("^nm","n",colnames(j))
      j <- j[,!(grepl("^chp",colnames(j)))]
      colnames(j) <- gsub("^chm","ch",colnames(j))
      j
    }))
    
    list(normalizedCountsTmpYesChpNpP=normalizedCountsTmpYesChpNpP,
         normalizedCountsTmpYesChpNpNoP=normalizedCountsTmpYesChpNpNoP,
         normalizedCountsTmpYesChpPNoNp=normalizedCountsTmpYesChpPNoNp,
         normalizedCountsTmpYesChpNoNpP=normalizedCountsTmpYesChpNoNpP,
         normalizedCountsTmpYesPNoChp=normalizedCountsTmpYesPNoChp,
         normalizedCountsTmpNoChpP=normalizedCountsTmpNoChpP)
  }
  
  normalizedCountsStatisticsSplitted <- splitGenes(normalizedCountsStatisticsSmall[2:3])
  
  list(rawCounts=readsAnnotationBackup
       ,DESeq2Objects=DESeq2Results
       ,DESeq2NormalizedCounts=DESeq2NormalizedCounts
       ,DESeq2Dispersions=DESeq2Dispersions
       ,normalizedCountsStatistics=normalizedCountsStatistics
       ,normalizedCountsStatisticsSplitted=normalizedCountsStatisticsSplitted)
}

spikesQuantification <- function(bamPath
                                 ,spikes="ERCC")
{
  bamTmp = get(load(bamPath))
  boolTmp = sapply(spikes,function(i)grepl(i,as.character(seqnames(bamTmp))))
  
  print(apply(boolTmp,2,sum))
  
  boolTmp = apply(boolTmp,1,any)
  bamTmp = bamTmp[boolTmp]
  return(unlist(as.list(table(as.character(seqnames(bamTmp))))))
}

### Models fit
##  - Cost function to be minimized,
##  - Infer rates function,
##  - Models simulation,
##  - Data generation.


## Cost function
costFunction <- function(par # Set of rates.
                         ,data # Experimental data.
                         ,dev # Experimental data standard deviation.
                         ,TauFractions = TauFractions # Time points with cellular fractionation.
                         ,TauPoly = TauPoly # Time points with polysomal profiling.
                         ,TauTotal = TauTotal # Time points with total RNA profiling.
                         ,dataGeneration # Data generation function.
                         ,parFixed=NULL # List of rates to be excluded from the optimization.
                         ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
                         ,excludeSpecies=NULL # List of species to be excluded from the cost function.
                         ,lowB=1e-6 # Lower bound for rates.
                         ,upB=1e4 # Upper bound for rates.
                         ,FlagDev="FC" # Cost function selection.
                         ,lambda=0.05) # Regularization strength.
{
  # If any standard deviation is null assume CV = 1 (true only for initial conditions where nascent RNA is absent by definition).
  dev[which(dev==0)] <- data[which(dev==0)]
  
  # If any parameter is fixed add it to the par list.
  if(!is.null(parFixed)){
    par <- c(par,parFixed)
    par <- par[sort(names(par))]
  }
  
  # From logarithmic to linear rates if needed.
  if(logOptim)
  {
    parL <- exp(par)
  }else{
    parL <- par
  }
  
  # Soft boundaries on rates values.
  if(any(parL<lowB)|any(parL>upB)|any(!is.finite(parL)))
  {
    return(1e8)
  }else{
    
    # Simulation of data according to the set of rates.
    simulatedData <- dataGeneration(exampleRates=parL # Parameters.
                                    ,TauFractions = TauFractions # Time points with cellular fractionation.
                                    ,TauPoly = TauPoly # Time points with polysomal profiling.
                                    ,TauTotal = TauTotal) # Time points with total RNA profiling.
    
    # Removal of species (e.g. Nascent Polysomal RNA)
    if(!is.null(excludeSpecies))
    {
      simulatedData <- simulatedData[!(names(simulatedData)%in%excludeSpecies)]
      data <- lapply(data,function(i)i[names(simulatedData)])
      dev <- lapply(dev,function(i)i[names(simulatedData)])
    }
    
    simulatedData <- simulatedData[!grepl("n0$",names(simulatedData))]
    
    # Exclusion of configurations with significantly negative expression levels.
    simulatedData[abs(simulatedData)<1e-6] <- 1e-6
    if(any(simulatedData<0))
    {
      return(1e8)
    }else{
      
      # Selection of the cost function.
      dev <- unlist(dev)
      data <- unlist(data)
      
      dev <- dev[data!=1e-10]
      data <- data[data!=1e-10]
      
      if(!all(names(simulatedData)%in%names(data))){print(names(simulatedData)[!(names(simulatedData)%in%names(data))])}
      
      simulatedData <- simulatedData[names(data)]
      
      x <- abs(log(simulatedData/data))
      return(switch(FlagDev
                    , M=sum(((simulatedData-data)**2)/data) + lambda*norm(par, type="2")
                    , M2=sum(((simulatedData-data)**2)/data**2) + lambda*norm(par, type="2")
                    , D=sum(((simulatedData-data)**2)/dev**2) + lambda*norm(par, type="2")
                    , FC=sum(x)+lambda*norm(par, type="2")
                    , LL=-sum(2*pnorm(-abs(simulatedData-data),0,dev,log.p=TRUE)))) # LL is not regularized because its is designed to be computed a posteriori on the optimal rates.
    }
  }
}

inferRates <- function(expressionData # Expression data of the genes to be modeled.
                       ,expressionDataDev # Standard deviations of the genes to be modeled.
                       ,simulatedDataset=NULL # Simulated dataset if this is the case (just to produce real rates correlations).
                       ,initialRates # List of initial rates for optimization.
                       ,TauFractions # Time points with cellular fractionation.
                       ,TauPoly # Time points with polysomal profiling.
                       ,TauTotal # Time points with total nascent and pre-existing RNA profiling.
                       ,cpus # Number of cpus.
                       ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
                       ,lowB=1e-6 # Lower boundary for the rates.
                       ,upB=1e4 # Upper boundary for the rates.
                       ,FlagDev="FC" # Cost function.
                       ,lambda=0.05 # Regularization strength.
                       ,excludeSpecies=NULL # List of species to be excluded from the cost function.
                       ,parFixed=NULL) # List of parameters to be excluded from the optimization.
{
  # Definition of the dataset as simulated or not.
  if(!is.null(simulatedDataset))
  {
    simulatedDatasetFlag <- TRUE
    initialRates <- lapply(initialRates,function(i)i[colnames(simulatedDataset$exampleRates)])
  }else{
    simulatedDatasetFlag <- FALSE
  }
  
  #  expressionDataTmp <- matrix(1e-10,nrow=length(allGenes),ncol=(ncol(expressionData[[1]])*length(expressionData)))
  #  rownames(expressionDataTmp) <- allGenes
  #  colnames(expressionDataTmp) <- unname(unlist(lapply(names(expressionData),function(i)paste0(i,"-",colnames(expressionData[[i]])))))
  #  
  #  for(i in names(expressionData))
  #  {
  #    print(i)
  #    expressionDataTmp[rownames(expressionData[[i]]),paste0(i,"-",colnames(expressionData[[i]]))] <- expressionData[[i]]
  #  }
  #
  #  expressionData <- expressionDataTmp
  
  # Sorting of expression levels.
  colOrderTmp <- sort(colnames(expressionData[[1]]))
  expressionData <- lapply(expressionData,function(i)i[,colOrderTmp])
  expressionDataDev <- lapply(expressionDataDev,function(i)i[,colOrderTmp])
  
  # If standard deviation are missing assume a CV of 1.
  if(is.null(expressionDataDev)){expressionDataDev<-expressionData}
  
  # If the optimization in the logarithmic space is required transform the initial rates.
  if(logOptim){initialRates <- lapply(initialRates,log)}
  
  # Selection of the configuration to be simulated according to the provided example rates or the RNA species.
  # This means definition of the modelSimulation function (ode implementation), of dataGeneration function (ode solution),
  # and in case of polysomal RNA profiling exclusion of nascent RNA from the species for the cost function optimization.
  
  if("k9"%in%names(unlist(initialRates))
     &any(grepl("^chpp",colnames(expressionData[[1]])))
     &any(grepl("^npp",colnames(expressionData[[1]])))
     &any(grepl("^py",colnames(expressionData[[1]]))))
  {
    dataGeneration <- dataGenerationFullNDM
    print("Inference mode: Mature nuclear RNA decay.")
  }else if("k11"%in%names(unlist(initialRates))
           &any(grepl("^chpp",colnames(expressionData[[1]])))
           &any(grepl("^npp",colnames(expressionData[[1]])))
           &any(grepl("^py",colnames(expressionData[[1]]))))
  {
    dataGeneration <- dataGenerationFullNDP
    print("Inference mode: Premature nuclear RNA decay.")
  }else if(any(grepl("^chpp",colnames(expressionData[[1]])))
           &any(grepl("^npp",colnames(expressionData[[1]])))
           &!any(grepl("^py",colnames(expressionData[[1]]))))
  {
    print("Inference mode: No polysomal RNA, Chromatin premature RNA, Nucleoplasmic premature RNA.")
    dataGeneration <- dataGenerationNoPoly
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k7","k9","k11","k10")]))
  }else if(any(grepl("^chpp",colnames(expressionData[[1]])))
           &any(grepl("^npp",colnames(expressionData[[1]])))
           &any(grepl("^py",colnames(expressionData[[1]]))))
  {
    print("Inference mode: Full model")
    dataGeneration <- dataGenerationFull
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k9","k11")]))
  }else if(any(grepl("^chpp",colnames(expressionData[[1]])))
           &!any(grepl("^npp",colnames(expressionData[[1]])))
           &!any(grepl("^py",colnames(expressionData[[1]]))))
  {
    print("Inference mode: No polysomal RNA, Chromatin premature RNA.")
    dataGeneration <- dataGenerationNoPolyNoNucP
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k4","k5","k7","k9","k11","k10")]))
  }else if(any(grepl("^chpp",colnames(expressionData[[1]])))
           &!any(grepl("^npp",colnames(expressionData[[1]])))
           &any(grepl("^py",colnames(expressionData[[1]]))))
  {
    print("Inference mode: Yes polysomal RNA, Chromatin premature RNA.")
    dataGeneration <- dataGenerationNoNucP
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k4","k5","k9","k11")]))
  }else if(!any(grepl("^chpp",colnames(expressionData[[1]])))
           &!any(grepl("^npp",colnames(expressionData[[1]])))
           &!any(grepl("^py",colnames(expressionData[[1]]))))
  {
    print("Inference mode: No polysomal RNA, No premature RNA.")
    dataGeneration <- dataGenerationNoPolyNoP
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k2","k4","k5","k7","k9","k11","k10")]))
  }else if(!any(grepl("^chpp",colnames(expressionData[[1]])))
           &!any(grepl("^npp",colnames(expressionData[[1]])))
           &any(grepl("^py",colnames(expressionData[[1]]))))
  {
    print("Inference mode: Yes polysomal RNA, No premature RNA.")
    dataGeneration <- dataGenerationNoP
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k2","k4","k5","k9","k11")]))
  }else{
    print("Issues in the input data!")
    return(NULL)
  }
  
  # Inference loop for any gene (row) in expressionData.
  inferedRates <- t(mcsapply(rownames(expressionData[[1]]),function(i)
  {
    print(i)
    # Expression levels isolation.
    exampleDataTmp <- lapply(expressionData,function(j)j[i,])
    exampleDataDevTmp <- lapply(expressionDataDev,function(j)j[i,])
    
    # Loop over three different algorithms.
    inferedRatesTmp <- sapply(c("L-BFGS-B","BFGS","Nelder-Mead"),function(m)
    {
      # Loop over all the initial rates configurations.
      foe <- sapply(initialRates,function(initialRatesTmp)
      {
        tryCatch(optim(par=initialRatesTmp # Initial rates.
                       ,FlagDev=FlagDev # Cost function selection.
                       ,fn=costFunction # Error function to be minimized.
                       ,dataGeneration=dataGeneration # Data generation function.
                       ,parFixed=parFixed # Fixed parameters.
                       ,data=exampleDataTmp # Experimental data.
                       ,dev=exampleDataDevTmp # Experimental data standard deviation.
                       ,lowB=lowB # Lower bound for rates.
                       ,upB=upB # Upper bound for rates.
                       ,TauFractions=TauFractions # Time points with cellular fractionation.
                       ,TauPoly=TauPoly # Time points with polysomal profiling.
                       ,TauTotal=TauTotal # Time points with total RNA profiling.
                       ,logOptim=logOptim # TRUE to optimize the model parameters in the Log space.
                       ,method=m # Optimization algorithm.
                       ,excludeSpecies=excludeSpecies # List of species to be excluded from the cost function.
                       ,lambda=lambda # Regularization strength.
                       ,control=list(maxit=1e4,reltol=1e-50)) # Optimization algorithm FIXED parameter.
                 ,error=function(e){list("par"=rep(NaN,length(initialRatesTmp))
                                         ,"value"=Inf
                                         ,"counts"=NaN
                                         ,"convergence"=NaN
                                         ,"message"="Error!")})
      })
      # Selection of the best model among the different initial conditions.
      foe[,which.min(foe["value",])]
    })
    
    # Selection of the best model among the different algorithms.
    winnerTmp <- names(which.min(inferedRatesTmp["value",]))
    inferedRatesTmp <- inferedRatesTmp[,winnerTmp]
    
    return(inferedRatesTmp$par)
  },mc.cores=cpus))
  
  # From logarithmic to linear rates.
  if(logOptim){inferedRates <- exp(inferedRates)}
  rownames(inferedRates) <- rownames(expressionData[[1]])
  
  # Generation of the expression levels for the optimal set of rates.
  inferedData <- t(sapply(1:nrow(inferedRates),function(i)dataGeneration(exampleRates=inferedRates[i,] # Rates.
                                                                         ,TauFractions=TauFractions # Time points with cellular fractionation.
                                                                         ,TauPoly=TauPoly # Time points with polysomal profiling.
                                                                         ,TauTotal=TauTotal))) # Time points with total nascent and pre-existing RNA profiling.
  rownames(inferedData) <- rownames(expressionData[[1]])
  
  # Loop to estimate for each gene the model Log Likelihood as a specific cost function.
  #  LogL <- t(sapply(1:nrow(expressionData),function(x){costFunction(par=inferedRates[x,] # Set of rates.
  #                                                                  ,data=expressionData[x,] # Experimental data.
  #                                                                  ,dev=expressionDataDev[x,] # Experimental data standard deviation.
  #                                                                  ,TauFractions = TauFractions # Time points with cellular fractionation.
  #                                                                  ,TauPoly = TauPoly # Time points with polysomal profiling.
  #                                                                  ,TauTotal = TauTotal # Time points with total RNA profiling.
  #                                                                  ,dataGeneration=dataGeneration # Data generation function.
  #                                                                  ,parFixed=NULL # List of rates to be excluded from the optimization.
  #                                                                  ,logOptim=FALSE # Always FALSE because we previously converted (potential) logarithmic rates.
  #                                                                  ,excludeSpecies=NULL # List of species to be excluded from the cost function.
  #                                                                  ,lowB=1e-6 # Lower bound for rates.
  #                                                                  ,upB=1e4 # Upper bound for rates.
  #                                                                  ,FlagDev="LL" # We are interested in LogLikelihood.
  #                                                                  ,lambda=0.05) # Regularization strength.
  #                                                      }))
  
  # AIC and BIC goodness of fit metrics.
  # AIC <- 2*ncol(inferedRates) - 2*LogL
  # BIC <- ncol(inferedRates)*log(ncol(expressionData))- 2*LogL
  
  # All goodness of fit metrics.
  # metrics <- rbind(LL=LogL,AIC=AIC,BIC=BIC)
  # rownames(metrics) <- c("LL","AIC","BIC")
  # colnames(metrics) <- rownames(expressionData)
  metrics <- matrix(1)
  
  # If data are simulated real rates correlations estimation.
  if(simulatedDatasetFlag)
  {
    realRates <- simulatedDataset$exampleRates
    
    print("Rates correlations:")
    print(round(cor(inferedRates,realRates,method="s"),2))
    
    cc <- round(cor(inferedRates,realRates,method="s"),2)
    
    relativeErrors <- cbind(inferedRates,realRates,(1-inferedRates/realRates))
    
    colnames(relativeErrors)[1:ncol(realRates)] <- paste0(colnames(realRates),".inf")
    colnames(relativeErrors)[(ncol(realRates)+1):(2*ncol(realRates))] <- paste0(colnames(realRates),".real")
    colnames(relativeErrors)[((2*ncol(realRates))+1):(3*ncol(realRates))] <- paste0(colnames(realRates),".err")
    
    idxs <- list(grep("inf",colnames(relativeErrors))
                 ,grep("real",colnames(relativeErrors))
                 ,grep("err",colnames(relativeErrors)))
    ratesCorrelations=round(sapply(idxs[[1]],function(i)cor(relativeErrors[,i],relativeErrors[,(max(idxs[[1]])+i)],method="s",use="c")),2)
  }else{relativeErrors=ratesCorrelations=NULL}
  
  # Output.
  list(inferedRates=inferedRates, inferedData=inferedData, expressionData=expressionData, metrics=t(metrics), relativeErrors=relativeErrors, ratesCorrelations=ratesCorrelations)
}

compareRates <- function(inferedRatesControl
                         ,inferedRatesTreatment
                         ,confidenceIntervalsControl
                         ,confidenceIntervalsTreatment
                         ,TauFractions
                         ,TauPoly
                         ,TauTotal
                         ,cpus
                         ,idx=NULL
                         ,width=8
                         ,height=6
                         ,fileNameLabel="")
{
  # Inferred data
  expressionDataControl <- inferedRatesControl$expressionData
  inferedDataControl <- inferedRatesControl$inferedData
  inferedDataControl[inferedDataControl==0] <- 1e-10
  inferedRatesControl <- inferedRatesControl$inferedRates
  
  if(!is.null(inferedRatesTreatment))
  {
    expressionDataTreatment <- inferedRatesTreatment$expressionData
    inferedDataTreatment <- inferedRatesTreatment$inferedData
    inferedDataTreatment[inferedDataTreatment==0] <- 1e-10
    inferedRatesTreatment <- inferedRatesTreatment$inferedRates
  }else{
    expressionDataTreatment <- expressionDataControl 
    inferedDataTreatment <- inferedRatesControl
    confidenceIntervalsTreatment <- confidenceIntervalsControl
    
    expressionDataTreatment[expressionDataTreatment!=+Inf] <- NaN
    inferedDataTreatment[inferedDataTreatment!=+Inf] <- NaN
    confidenceIntervalsTreatment <- lapply(confidenceIntervalsTreatment,function(i){i[i!=+Inf] <- NaN;i})
  }
  
  if(!is.null(idx))
  {
    expressionDataControlTmp <- expressionDataControl[idx,]
    expressionDataControlTmp <- expressionDataControlTmp[!grepl("0$",names(expressionDataControlTmp))]
    
    expressionDataTreatmentTmp <- expressionDataTreatment[idx,]
    expressionDataTreatmentTmp <- expressionDataTreatmentTmp[!grepl("0$",names(expressionDataTreatmentTmp))]
    
    inferedDataControlTmp <- inferedDataControl[idx,]
    inferedDataControlTmp <- inferedDataControlTmp[!grepl("0$",names(inferedDataControlTmp))]
    
    inferedDataTreatmentTmp <- inferedDataTreatment[idx,]
    inferedDataTreatmentTmp <- inferedDataTreatmentTmp[!grepl("0$",names(inferedDataTreatmentTmp))]
    
    confidenceIntervalsControl <- sapply(confidenceIntervalsControl,function(i)i[idx,!grepl("0$",colnames(i))])
    confidenceIntervalsTreatment <- sapply(confidenceIntervalsTreatment,function(i)i[idx,!grepl("0$",colnames(i))])
    
    x <- matrix(NaN,nrow=6,ncol=nrow(confidenceIntervalsControl))
    x[1,] <- confidenceIntervalsControl[,"mean"]
    x[3,] <- confidenceIntervalsTreatment[,"mean"]
    colTmp <- rep(1,nrow(x));colTmp[3] <- 2
    x <- log10(c(x))
    
    x2 <- matrix(NaN,nrow=6,ncol=length(inferedDataControlTmp))
    x2[1,] <- inferedDataControlTmp[rownames(confidenceIntervalsControl)]
    x2[3,] <- inferedDataTreatmentTmp[rownames(confidenceIntervalsTreatment)]
    colTmp <- rep(1,nrow(x2));colTmp[3] <- 2
    x2 <- log10(c(x2))
    
    labelsTmp <- c(sapply(rownames(confidenceIntervalsControl),function(i){c("",i,rep("",4))}))
    
    y <- matrix(NaN,nrow=6,ncol=nrow(confidenceIntervalsControl))
    y[1,] <- confidenceIntervalsControl[,"bottom"]
    y[3,] <- confidenceIntervalsTreatment[,"bottom"]
    y <- log10(c(y))
    
    y2 <- (y==-10)
    y[y2] <- 0.9*min(y[!y2],na.rm=TRUE)
    
    z <- matrix(NaN,nrow=6,ncol=nrow(confidenceIntervalsControl))
    z[1,] <- confidenceIntervalsControl[,"up"]
    z[3,] <- confidenceIntervalsTreatment[,"up"]
    z <- log10(c(z))
    
    lineTypeTmp <- rep(1,length(y))
    lineTypeTmp[y2] <- 2
    
    pdf(paste0("CIs_",idx,fileNameLabel,".pdf"),width=width,height=height)
    par(mar=c(7.1, 4.1, 4.1, 2.1))
    plot(x,pch=16,col=colTmp,ylim=c(min(c(y[is.finite(y)],x2[is.finite(x2)])),max(c(z[is.finite(z)],x2[is.finite(x2)]))),xaxt = 'n',xlab="",main=idx,ylab="Log10 [fg/M cells]")
    points(x2,pch=8,col=colTmp)
    axis(1, at=seq_along(x), labels = FALSE, lwd.ticks = 1*(labelsTmp!=""))
    text(seq_along(x), par("usr")[3] - 0.2, labels = labelsTmp, srt = 90, pos = 1, xpd = TRUE)
    arrows(x0=seq_along(x),y0=x,x1=seq_along(x),y1=y,length=0,col=colTmp,lty=lineTypeTmp)
    arrows(x0=seq_along(x),y0=x,x1=seq_along(x),y1=z,length=0,col=colTmp)
    dev.off()
    return(NULL)
  }else{
    if("k9"%in%colnames(inferedRatesControl))
    { #Full Model with nuclear mature RNA decay
      modelSimulation <- fullModelSimulationNDM
      dataGeneration <- dataGenerationFullNDM
      print("Full Model with nuclear mature RNA decay")
    }else if("k11"%in%colnames(inferedRatesControl))
    { #Full Model with nuclear premature RNA decay
      modelSimulation <- fullModelSimulationNDP
      dataGeneration <- dataGenerationFullNDP
      print("Full Model with nuclear premature RNA decay")
    }else if(all(c("k1","k2","k3","k4","k5","k6","k7","k8","k10")%in%colnames(inferedRatesControl)))
    { #Full Model
      modelSimulation <- fullModelSimulation
      dataGeneration <- dataGenerationFull
      print("Full Model")
    }else if(all(c("k1","k2","k3","k4","k5","k6","k8")%in%colnames(inferedRatesControl)))
    { #Full Model without translation
      modelSimulation <- modelSimulationNoPoly
      dataGeneration <- dataGenerationNoPoly
      print("Full Model without translation")
    }else if(all(c("k1","k2","k3","k6","k7","k8","k10")%in%colnames(inferedRatesControl)))
    { #No Nucleoplasmic Premature - yes tranlation
      modelSimulation <- modelSimulationNoNucP
      dataGeneration <- dataGenerationNoNucP
      print("No Nucleoplasmic Premature - yes tranlation")
    }else if(all(c("k1","k2","k3","k6","k8")%in%colnames(inferedRatesControl)))
    { #No Nucleoplasmic Premature - no tranlation
      modelSimulation <- modelSimulationNoPolyNoNucP
      dataGeneration <- dataGenerationNoPolyNoNucP
      print("No Nucleoplasmic Premature - no tranlation")
    }else if(all(c("k1","k3","k6","k7","k8","k10")%in%colnames(inferedRatesControl)))
    { #No premature - yes translation
      modelSimulation <- modelSimulationNoP
      dataGeneration <- dataGenerationNoP
      print("No premature - yes translation")
    }else if(all(c("k1","k3","k6","k8")%in%colnames(inferedRatesControl)))
    { #No premature - no translation
      modelSimulation <- modelSimulationNoPolyNoP
      dataGeneration <- dataGenerationNoPolyNoP
      print("No premature - no translation")
    }else{
      print("Error in configuration!")
    }
    
    # Confidence intervals re-shaping
    # confidenceIntervalsControl <- lapply(confidenceIntervalsControl,function(i)
    # {
    #     i <- i[,!(grepl("^cyp",colnames(i))|grepl("^pyp",colnames(i)))]
    #     colnames(i) <- gsub("^cym","cy",colnames(i))
    #     colnames(i) <- gsub("^pym","py",colnames(i))
    #   i[rownames(inferedDataControl),colnames(inferedDataControl)]
    # })
    
    # confidenceIntervalsTreatment <- lapply(confidenceIntervalsTreatment,function(i)
    # {
    #     i <- i[,!(grepl("^cyp",colnames(i))|grepl("^pyp",colnames(i)))]
    #     colnames(i) <- gsub("^cym","cy",colnames(i))
    #     colnames(i) <- gsub("^pym","py",colnames(i))
    #     i[rownames(inferedDataTreatment),colnames(inferedDataTreatment)]
    # })
    
    # Genes selection
    selectedGenesControl <- names(which(apply((inferedDataControl>=confidenceIntervalsControl$bottom[,colnames(inferedDataControl)]
                                               &inferedDataControl<=confidenceIntervalsControl$up[,colnames(inferedDataControl)]),1,all)))
    
    selectedGenesTreatment <- names(which(apply((inferedDataTreatment>=confidenceIntervalsTreatment$bottom[,colnames(inferedDataTreatment)]
                                                 &inferedDataTreatment<=confidenceIntervalsTreatment$up[,colnames(inferedDataTreatment)]),1,all)))
    
    # foe <- (confidenceIntervalsControl$up - confidenceIntervalsControl$bottom)/confidenceIntervalsControl$mean
    # foe <- foe[,!grepl("0$",colnames(foe))]
    
    commonGenes <- intersect(selectedGenesControl,selectedGenesTreatment)
    
    print(paste0(length(commonGenes)," genes profiled in both treatments."))
    
    # Differential rates
    differentialRatesSimple <- lapply(commonGenes,function(i)
    {
      sapply(colnames(inferedRatesTreatment),function(j)
      {
        inferedRatesTmp <- inferedRatesControl[i,]
        inferedRatesTmp[j] <- inferedRatesTreatment[i,j]
        
        foe <- dataGeneration(inferedRatesTmp,TauFractions,TauPoly,TauTotal)
        foe <- (foe<confidenceIntervalsControl$bottom[i,names(foe)]|foe>confidenceIntervalsControl$up[i,names(foe)])
        foe[!grepl("0$",names(foe))]
      })
    })
    names(differentialRatesSimple) <- commonGenes
    
    differentialRates <- sapply(colnames(inferedRatesTreatment),function(j)
    {
      print(j)
      inferedRatesTmp <- inferedRatesControl[commonGenes,]
      inferedRatesTmp[commonGenes,j] <- inferedRatesTreatment[commonGenes,j]
      
      foe <- inferRatesCI(expressionData=expressionDataControl[commonGenes,]
                          ,expressionDataCIDown=confidenceIntervalsControl$bottom[commonGenes,colnames(expressionDataControl)]
                          ,expressionDataCIUp=confidenceIntervalsControl$up[commonGenes,colnames(expressionDataControl)]
                          ,initialRates=inferedRatesTmp[commonGenes,]
                          ,TauFractions=TauFractions
                          ,TauPoly=TauPoly
                          ,TauTotal=TauTotal
                          ,cpus=cpus
                          ,logOptim=TRUE
                          ,lowB=1e-6
                          ,upB=1e10
                          ,FlagDev="FC"
                          ,lambda=0.05
                          ,excludeSpecies=NULL
                          ,parFixed=j)
      
    })
    names(differentialRates) <- unlist(lapply(colnames(inferedRatesTreatment),function(i)paste0(i,c("_rates","_expressions"))))
    differentialRates <- lapply(differentialRates,function(i){i[i==0]<-1e-10;i})
    
    differentialRatesTable <- lapply(differentialRates[grepl("expressions",names(differentialRates))],function(l)
    {
      (l<confidenceIntervalsControl$bottom[rownames(l),colnames(l)]
       |l>confidenceIntervalsControl$up[rownames(l),colnames(l)])
    })
    names(differentialRatesTable) <- sapply(strsplit(names(differentialRatesTable),"_"),"[[",1)
    
    table(sapply(differentialRatesSimple,any))
    # FALSE  TRUE 
    #    27    60 
    
    table(apply(Reduce('+',differentialRatesTable)!=0,1,any))
    # FALSE  TRUE 
    #    48    39 
    
    # saveRDS(differentialRatesSimple,"differentialRatesSimple.rds")
    # saveRDS(differentialRates,"differentialRates.rds")
    
    differentialRatesTable <- lapply(differentialRatesTable,function(l)l[,!grepl("0$",colnames(l))])
    
    return(list(differentialRatesSimple=differentialRatesSimple
                ,differentialRates=differentialRates
                ,differentialRatesTable=differentialRatesTable
                ,selectedGenesControl=selectedGenesControl
                ,selectedGenesTreatment=selectedGenesTreatment
                ,commonGenes=commonGenes))
  }
}

inferRatesCI <- function(expressionData # Expression data of the genes to be modeled.
                         ,expressionDataDev # Standard deviations of the genes to be modeled.
                         ,initialRates # List of initial rates for optimization.
                         ,TauFractions # Time points with cellular fractionation.
                         ,TauPoly # Time points with polysomal profiling.
                         ,TauTotal # Time points with total nascent and pre-existing RNA profiling.
                         ,cpus # Number of cpus.
                         ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
                         ,lowB=1e-6 # Lower boundary for the rates.
                         ,upB=1e4 # Upper boundary for the rates.
                         ,FlagDev="FC" # Cost function.
                         ,lambda=0.05 # Regularization strength.
                         ,excludeSpecies=NULL # List of species to be excluded from the cost function.
                         ,parFixed=NULL) # List of parameters to be excluded from the optimization.
{
  # Sorting of expression levels.
  colOrderTmp <- sort(colnames(expressionData))
  expressionData <- expressionData[,colOrderTmp]
  expressionDataDev <- expressionDataDev[,colOrderTmp]
  
  # If standard deviation are missing assume a CV of 1.
  if(is.null(expressionDataDev)){expressionDataDev<-expressionData}
  
  # If the optimization in the logarithmic space is required transform the initial rates.
  if(logOptim){initialRates <- log(initialRates)}
  
  # Selection of the configuration to be simulated according to the provided example rates or the RNA species.
  # This means definition of the modelSimulation function (ode implementation), of dataGeneration function (ode solution),
  # and in case of polysomal RNA profiling exclusion of nascent RNA from the species for the cost function optimization.
  
  if("k9"%in%names(unlist(initialRates))
     &any(grepl("^chpp",colnames(expressionData)))
     &any(grepl("^npp",colnames(expressionData)))
     &any(grepl("^py",colnames(expressionData))))
  {
    dataGeneration <- dataGenerationFullNDM
    # excludeSpecies=colnames(expressionData[,c(grep("^pyn",colnames(expressionData))
    #                                          ,grep("^pyp",colnames(expressionData)))])
    print("Inference mode: Nuclear mature RNA decay.")
  }else if("k11"%in%names(unlist(initialRates))
           &any(grepl("^chpp",colnames(expressionData)))
           &any(grepl("^npp",colnames(expressionData)))
           &any(grepl("^py",colnames(expressionData))))
  {
    dataGeneration <- dataGenerationFullNDP
    print("Inference mode: Nuclear premature RNA decay.")
  }else if(any(grepl("^chpp",colnames(expressionData)))
           &any(grepl("^npp",colnames(expressionData)))
           &!any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: No polysomal RNA, Chromatin premature RNA, Nucleoplasmic premature RNA.")
    dataGeneration <- dataGenerationNoPoly
    initialRates <- initialRates[,!colnames(initialRates)%in%c("k7","k9","k10")]
  }else if(any(grepl("^chpp",colnames(expressionData)))
           &any(grepl("^npp",colnames(expressionData)))
           &any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: Full model")
    dataGeneration <- dataGenerationFull
    initialRates <- initialRates[,!colnames(initialRates)%in%c("k9")]
    # excludeSpecies=colnames(expressionData[,c(grep("^pyn",colnames(expressionData))
    #                                          ,grep("^pyp",colnames(expressionData)))])
  }else if(any(grepl("^chpp",colnames(expressionData)))
           &!any(grepl("^npp",colnames(expressionData)))
           &!any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: No polysomal RNA, Chromatin premature RNA.")
    dataGeneration <- dataGenerationNoPolyNoNucP
    initialRates <- initialRates[,!colnames(initialRates)%in%c("k4","k5","k7","k9","k10")]
  }else if(any(grepl("^chpp",colnames(expressionData)))
           &!any(grepl("^npp",colnames(expressionData)))
           &any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: Yes polysomal RNA, Chromatin premature RNA.")
    dataGeneration <- dataGenerationNoNucP
    # excludeSpecies=colnames(expressionData[,c(grep("^pyn",colnames(expressionData))
    #                                          ,grep("^pyp",colnames(expressionData)))])
    initialRates <- initialRates[,!colnames(initialRates)%in%c("k4","k5","k9")]
  }else if(!any(grepl("^chpp",colnames(expressionData)))
           &!any(grepl("^npp",colnames(expressionData)))
           &!any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: No polysomal RNA, No premature RNA.")
    dataGeneration <- dataGenerationNoPolyNoP
    initialRates <- initialRates[,!colnames(initialRates)%in%c("k2","k4","k5","k7","k9","k10")]
  }else if(!any(grepl("^chpp",colnames(expressionData)))
           &!any(grepl("^npp",colnames(expressionData)))
           &any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: Yes polysomal RNA, No premature RNA.")
    dataGeneration <- dataGenerationNoP
    # excludeSpecies=colnames(expressionData[,c(grep("^pyn",colnames(expressionData))
    #                                          ,grep("^pyp",colnames(expressionData)))])
    initialRates <- initialRates[,!colnames(initialRates)%in%c("k2","k4","k5","k9")]
  }else{
    print("Issues in the input data!")
    return(NULL)
  }
  
  # Inference loop for any gene (row) in expressionData.
  inferedRates <- t(mcsapply(1:nrow(expressionData),function(i)
  {
    print(i)
    
    # Expression levels isolation.
    exampleDataTmp <- unlist(expressionData[i,])
    exampleDataDevTmp <- unlist(expressionDataDev[i,])
    
    initialRatesTmp <- initialRates[i,]
    
    # Loop over all the initial rates configurations.
    foe <- sapply(c("L-BFGS-B","BFGS","Nelder-Mead"),function(m)
    {
      tryCatch(optim(par=initialRatesTmp[names(initialRatesTmp)!=parFixed] # Initial rates.
                     ,FlagDev=FlagDev # Cost function selection.
                     ,fn=costFunction # Error function to be minimized.
                     ,dataGeneration=dataGeneration # Data generation function.
                     ,parFixed=initialRatesTmp[names(initialRatesTmp)%in%parFixed] # Fixed parameters.
                     ,data=exampleDataTmp # Experimental data.
                     ,dev=exampleDataDevTmp # Experimental data standard deviation.
                     ,lowB=lowB # Lower bound for rates.
                     ,upB=upB # Upper bound for rates.
                     ,TauFractions=TauFractions # Time points with cellular fractionation.
                     ,TauPoly=TauPoly # Time points with polysomal profiling.
                     ,TauTotal=TauTotal # Time points with total RNA profiling.
                     ,logOptim=logOptim # TRUE to optimize the model parameters in the Log space.
                     ,method=m # Optimization algorithm.
                     ,excludeSpecies=excludeSpecies # List of species to be excluded from the cost function.
                     ,lambda=lambda # Regularization strength.
                     ,control=list(maxit=1e9,reltol=1e-50)) # Optimization algorithm FIXED parameter.
                     ,error=function(e){list("par"=rep(NaN,length(initialRatesTmp))
                                       ,"value"=Inf
                                       ,"counts"=NaN
                                       ,"convergence"=NaN
                                       ,"message"="Error!")})
    })
    
    # Selection of the best model among the different algorithms.
    winnerTmp <- names(which.min(foe["value",]))
    foe <- foe[,winnerTmp]
    out <- c(foe$par,initialRatesTmp[names(initialRatesTmp)%in%parFixed])
    out <- out[sort(names(out),index.return=TRUE)[[2]]]
    return(out)
  },mc.cores=cpus))
  
  # From logarithmic to linear rates.
  if(logOptim){inferedRates <- exp(inferedRates)}
  rownames(inferedRates) <- rownames(expressionData)
  
  # Generation of the expression levels for the optimal set of rates.
  inferedData <- t(sapply(1:nrow(inferedRates),function(i)dataGeneration(exampleRates=inferedRates[i,] # Rates.
                                                                         ,TauFractions=TauFractions # Time points with cellular fractionation.
                                                                         ,TauPoly=TauPoly # Time points with polysomal profiling.
                                                                         ,TauTotal=TauTotal))) # Time points with total nascent and pre-existing RNA profiling.
  rownames(inferedData) <- rownames(inferedRates)
  
  # Output.
  list(inferedRates=inferedRates, inferedData=inferedData)
}

## Full model with nuclear mature decay
fullModelSimulationNDM <- function(t,y,parms)
{
  k1=parms[["k1"]] # Synthesis.
  k2=parms[["k2"]] # Co-transcriptional splicing.
  k3=parms[["k3"]] # Detachment of spliced chromatin associated RNA.
  k4=parms[["k4"]] # Detachment of unspliced chromatin associated RNA.
  k5=parms[["k5"]] # Post-transcriptional splicing.
  k6=parms[["k6"]] # Export.
  k7=parms[["k7"]] # Translation.
  k8=parms[["k8"]] # Cytoplasmic decay
  k9=parms[["k9"]] # Mature nuclear RNA decay.
  k10=parms[["k10"]] # Polysomal decay
  
  # Implementation of the ODE model.
  # System of equations for NASCENT RNA
  dy1=k1 - (k2+k4)*y[1]                     #y1=Chp
  dy2=k2*y[1] - k3*y[2]                     #y2=Chm
  dy3=k4*y[1] - k5*y[3]                     #y3=Np
  dy4=k3*y[2] + k5*y[3] - (k6+k9)*y[4]      #y4=Nm
  dy5=k6*y[4] - (k7+k8)*y[5]                #y5=C
  dy6=k7*y[5] - k10*y[6]                     #y6=P             
  
  nascentTmp <- c(dy1,dy2,dy3,dy4,dy5,dy6)
  
  # System of equations for PRE-EXISTING RNA
  dy7=0 - (k2+k4)*y[7]                       #y7=Chp
  dy8=k2*y[7] - k3*y[8]                      #y8=Chm
  dy9=k4*y[7] - k5*y[9]                      #y9=Np
  dy10=k3*y[8] + k5*y[9] - (k6+k9)*y[10]     #y10=Nm
  dy11=k6*y[10] - (k7+k8)*y[11]              #y11=C
  dy12=k7*y[11] - k10*y[12]                   
  
  preExistingTmp <- c(dy7,dy8,dy9,dy10,dy11,dy12)
  
  # Output.
  list(c(nascentTmp,preExistingTmp))
}

dataGenerationFullNDM <- function(exampleRates # Set of rate.
                                  ,TauFractions # Time points with cellular fractionation.
                                  ,TauPoly # Time points with polysomal profiling.
                                  ,TauTotal) # Time points with total RNA profiling.
{
  k1=exampleRates[["k1"]] # Synthesis.
  k2=exampleRates[["k2"]] # Co-transcriptional splicing.
  k3=exampleRates[["k3"]] # Detachment of spliced chromatin associated RNA.
  k4=exampleRates[["k4"]] # Detachment of unspliced chromatin associated RNA.
  k5=exampleRates[["k5"]] # Post-transcriptional splicing.
  k6=exampleRates[["k6"]] # Export.
  k7=exampleRates[["k7"]] # Translation.
  k8=exampleRates[["k8"]] # Cytoplasmic decay
  k9=exampleRates[["k9"]] # Mature nuclear RNA decay.
  k10=exampleRates[["k10"]] # Polysomal decay
  
  # Initial conditions.
  yini <- c(y1=0,y2=0,y3=0,y4=0,y5=0,y6=0
            ,y7= k1/(k2+k4)
            ,y8= k1*k2/(k3*(k2+k4))
            ,y9= k4*k1/(k5*(k2+k4))
            ,y10= k1/(k6+k9)
            ,y11= k1*k6/((k7+k8)*(k6+k9))
            ,y12= k7*k1*k6/(k10*(k6+k9)*(k7+k8)))
  
  # All time-points to be simulated
  TauGlobals <- unique(c(TauFractions,TauPoly,TauTotal))
  
  # Numerical solution of the ODE model.
  exampleData <- ode(times=TauGlobals # Time points.
                     ,y=yini # Initial conditions.
                     ,func=fullModelSimulationNDM # Model.
                     ,parms=exampleRates) # Rates.
  
  # Removal of time points (first column).
  exampleData <- c(t(exampleData[,-1]))
  
  # Set of colnames.
  names(exampleData) <- c(sapply(TauGlobals,function(i)paste0(c("chpn","chmn","npn","nmn","cyn","pyn"
                                                                ,"chpp","chmp","npp","nmp","cyp","pyp"),i)))
  
  # Steady state solution (initial conditions).
  exampleData <- c(tail(yini,(length(yini)/2)),exampleData)
  names(exampleData)[1:(length(yini)/2)] <- c("chpt","chmt","npt","nmt","cyt","pyt")
  
  # Cytoplasmic RNA is re-defined as the sum of polysomal RNA and cytoplasmic RNA from the model
  # because this is closer to the experimental data.
  # Here we could add a flag to avoid this step if this assumption is not true for a specific experiment.
  exampleData[sort(names(exampleData[grep("^cyn",names(exampleData))]))] <- exampleData[sort(names(exampleData[grep("^cyn",names(exampleData))]))]+exampleData[sort(names(exampleData[grep("^pyn",names(exampleData))]))] 
  exampleData[sort(names(exampleData[grep("^cyp",names(exampleData))]))]  <- exampleData[sort(names(exampleData[grep("^cyp",names(exampleData))]))]+exampleData[sort(names(exampleData[grep("^pyp",names(exampleData))]))] 
  exampleData["cyt"] <- exampleData["cyt"]+exampleData["pyt"]
  
  exampleData <- c(exampleData[grep("t$",names(exampleData))] # We keep the total fractions
                   ,unlist(lapply(c(TauFractions,TauTotal),function(i)exampleData[grepl(paste0(i,"$"),names(exampleData))
                                                                                  &!grepl(paste0("^pyn",i,"$"),names(exampleData))
                                                                                  &!grepl(paste0("^pyp",i,"$"),names(exampleData))])) # We keep all the data for the profiled fraction (no polysomal which can be profiled at different time points).
                   ,unlist(lapply(TauPoly,function(i)exampleData[grepl(paste0("^pyn",i,"$"),names(exampleData))
                                                                 |grepl(paste0("^pyp",i,"$"),names(exampleData))])))
  
  # If some fractionation time-points have been profiled only for total RNA we compute the total nascent and total
  # pre-existing, and we remove the corresponding fractions.
  if(!is.null(TauTotal))
  {
    exampleData <- c(exampleData
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("n",i,"$"),names(exampleData))
                                                                  &!grepl(paste0("^pyn",i,"$"),names(exampleData))])}) # For certain time-points we just have nascent total RNA.
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("p",i,"$"),names(exampleData))
                                                                  &!grepl(paste0("^pyp",i,"$"),names(exampleData))])})) # For certain time-points we just have pre-existing total RNA.
    
    names(exampleData)[names(exampleData)==""] <- c(paste0("n",TauTotal),paste0("p",TauTotal))
    
    # Removal of fractions for time-points profiled without cellular fractionation.
    exampleData <- exampleData[!(names(exampleData)%in%c(sapply(TauTotal,function(i)paste0(c("chpn","chmn","npn","nmn","cyn","chpp","chmp","npp","nmp","cyp"),i))))]
    
    # Removal of eventual duplicated conditions
    exampleData <- exampleData[!duplicated(names(exampleData))]
  }
  
  exampleData <- exampleData[!duplicated(names(exampleData))]
  exampleData <- exampleData[sort(names(exampleData))]
}

## Full model with nuclear premature decay
fullModelSimulationNDP <- function(t,y,parms)
{
  k1=parms[["k1"]] # Synthesis.
  k2=parms[["k2"]] # Co-transcriptional splicing.
  k3=parms[["k3"]] # Detachment of spliced chromatin associated RNA.
  k4=parms[["k4"]] # Detachment of unspliced chromatin associated RNA.
  k5=parms[["k5"]] # Post-transcriptional splicing.
  k6=parms[["k6"]] # Export.
  k7=parms[["k7"]] # Translation.
  k8=parms[["k8"]] # Cytoplasmic decay
  k10=parms[["k10"]] # Polysomal decay
  k11=parms[["k11"]] # Premature nuclear RNA decay.
  
  # Implementation of the ODE model.
  # System of equations for NASCENT RNA
  dy1=k1 - (k2+k4)*y[1]                     #y1=Chp
  dy2=k2*y[1] - k3*y[2]                     #y2=Chm
  dy3=k4*y[1] - (k5+k11)*y[3]               #y3=Np
  dy4=k3*y[2] + k5*y[3] - k6*y[4]           #y4=Nm
  dy5=k6*y[4] - (k7+k8)*y[5]                #y5=C
  dy6=k7*y[5] - k10*y[6]                    #y6=P             
  
  nascentTmp <- c(dy1,dy2,dy3,dy4,dy5,dy6)
  
  # System of equations for PRE-EXISTING RNA
  dy7=0 - (k2+k4)*y[7]                       #y7=Chp
  dy8=k2*y[7] - k3*y[8]                      #y8=Chm
  dy9=k4*y[7] - (k5+k11)*y[9]                #y9=Np
  dy10=k3*y[8] + k5*y[9] - k6*y[10]          #y10=Nm
  dy11=k6*y[10] - (k7+k8)*y[11]              #y11=C
  dy12=k7*y[11] - k10*y[12]                   
  
  preExistingTmp <- c(dy7,dy8,dy9,dy10,dy11,dy12)
  
  # Output.
  list(c(nascentTmp,preExistingTmp))
}

dataGenerationFullNDP <- function(exampleRates # Set of rate.
                                  ,TauFractions # Time points with cellular fractionation.
                                  ,TauPoly # Time points with polysomal profiling.
                                  ,TauTotal) # Time points with total RNA profiling.
{
  k1=exampleRates[["k1"]] # Synthesis.
  k2=exampleRates[["k2"]] # Co-transcriptional splicing.
  k3=exampleRates[["k3"]] # Detachment of spliced chromatin associated RNA.
  k4=exampleRates[["k4"]] # Detachment of unspliced chromatin associated RNA.
  k5=exampleRates[["k5"]] # Post-transcriptional splicing.
  k6=exampleRates[["k6"]] # Export.
  k7=exampleRates[["k7"]] # Translation.
  k8=exampleRates[["k8"]] # Cytoplasmic decay
  k10=exampleRates[["k10"]] # Polysomal decay
  k11=exampleRates[["k11"]] # Premature nuclear RNA decay.
  
  # Initial conditions.
  yini <- c(y1=0,y2=0,y3=0,y4=0,y5=0,y6=0
            ,y7= k1/(k2+k4)
            ,y8= k1*k2/(k3*(k2+k4))
            ,y9= k4*k1/((k2+k4)*(k5+k11))
            ,y10= (k1/(k6*(k2+k4)))*(k2+(k4*k5)/(k5+k11))
            ,y11= (k1/((k7+k8)*(k2+k4)))*(k2+(k4*k5)/(k5+k11))
            ,y12= (k7/k10)*(k1/((k7+k8)*(k2+k4)))*(k2+(k4*k5)/(k5+k11)))
  
  # All time-points to be simulated
  TauGlobals <- unique(c(TauFractions,TauPoly,TauTotal))
  
  # Numerical solution of the ODE model.
  exampleData <- ode(times=TauGlobals # Time points.
                     ,y=yini # Initial conditions.
                     ,func=fullModelSimulationNDP # Model.
                     ,parms=exampleRates) # Rates.
  
  # Removal of time points (first column).
  exampleData <- c(t(exampleData[,-1]))
  
  # Set of colnames.
  names(exampleData) <- c(sapply(TauGlobals,function(i)paste0(c("chpn","chmn","npn","nmn","cyn","pyn"
                                                                ,"chpp","chmp","npp","nmp","cyp","pyp"),i)))
  
  # Steady state solution (initial conditions).
  exampleData <- c(tail(yini,(length(yini)/2)),exampleData)
  names(exampleData)[1:(length(yini)/2)] <- c("chpt","chmt","npt","nmt","cyt","pyt")
  
  # Cytoplasmic RNA is re-defined as the sum of polysomal RNA and cytoplasmic RNA from the model
  # because this is closer to the experimental data.
  # Here we could add a flag to avoid this step if this assumption is not true for a specific experiment.
  exampleData[sort(names(exampleData[grep("^cyn",names(exampleData))]))] <- exampleData[sort(names(exampleData[grep("^cyn",names(exampleData))]))]+exampleData[sort(names(exampleData[grep("^pyn",names(exampleData))]))] 
  exampleData[sort(names(exampleData[grep("^cyp",names(exampleData))]))]  <- exampleData[sort(names(exampleData[grep("^cyp",names(exampleData))]))]+exampleData[sort(names(exampleData[grep("^pyp",names(exampleData))]))] 
  exampleData["cyt"] <- exampleData["cyt"]+exampleData["pyt"]
  
  exampleData <- c(exampleData[grep("t$",names(exampleData))] # We keep the total fractions
                   ,unlist(lapply(c(TauFractions,TauTotal),function(i)exampleData[grepl(paste0(i,"$"),names(exampleData))
                                                                                  &!grepl(paste0("^pyn",i,"$"),names(exampleData))
                                                                                  &!grepl(paste0("^pyp",i,"$"),names(exampleData))])) # We keep all the data for the profiled fraction (no polysomal which can be profiled at different time points).
                   ,unlist(lapply(TauPoly,function(i)exampleData[grepl(paste0("^pyn",i,"$"),names(exampleData))
                                                                 |grepl(paste0("^pyp",i,"$"),names(exampleData))])))
  
  # If some fractionation time-points have been profiled only for total RNA we compute the total nascent and total
  # pre-existing, and we remove the corresponding fractions.
  if(!is.null(TauTotal))
  {
    exampleData <- c(exampleData
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("n",i,"$"),names(exampleData))
                                                                  &!grepl(paste0("^pyn",i,"$"),names(exampleData))])}) # For certain time-points we just have nascent total RNA.
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("p",i,"$"),names(exampleData))
                                                                  &!grepl(paste0("^pyp",i,"$"),names(exampleData))])})) # For certain time-points we just have pre-existing total RNA.
    
    names(exampleData)[names(exampleData)==""] <- c(paste0("n",TauTotal),paste0("p",TauTotal))
    
    # Removal of fractions for time-points profiled without cellular fractionation.
    exampleData <- exampleData[!(names(exampleData)%in%c(sapply(TauTotal,function(i)paste0(c("chpn","chmn","npn","nmn","cyn","chpp","chmp","npp","nmp","cyp"),i))))]
    
    # Removal of eventual duplicated conditions
    exampleData <- exampleData[!duplicated(names(exampleData))]
  }
  
  exampleData <- exampleData[!duplicated(names(exampleData))]
  exampleData <- exampleData[sort(names(exampleData))]
}

## Full model
fullModelSimulation <- function(t,y,parms)
{
  k1=parms[["k1"]] # Synthesis
  k2=parms[["k2"]]# Co-transcriptional splicing
  k3=parms[["k3"]] # Detachment of spliced chromatin
  k4=parms[["k4"]] # Detachment of unspliced chromatin
  k5=parms[["k5"]] # Post-transcriptional splicing
  k6=parms[["k6"]] # Export
  k7=parms[["k7"]] # Translation
  k8=parms[["k8"]] # Cytoplasmic decay
  k10=parms[["k10"]] # Polysomal decay
  
  # System of equations for NASCENT RNA
  dy1=k1 - (k2+k4)*y[1]                     #y1=Chp
  dy2=k2*y[1] - k3*y[2]                     #y2=Chm
  dy3=k4*y[1] - k5*y[3]                     #y3=Np
  dy4=k3*y[2] + k5*y[3] - k6*y[4]           #y4=Nm
  dy5=k6*y[4] - (k7+k8)*y[5]                #y5=C
  dy6=k7*y[5] - k10*y[6]                     #y6=P             
  
  nascentTmp <- c(dy1,dy2,dy3,dy4,dy5,dy6)
  
  # System of equations for PRE-EXISTING RNA
  dy7=0 - (k2+k4)*y[7]                       #y7=Chp
  dy8=k2*y[7] - k3*y[8]                      #y8=Chm
  dy9=k4*y[7] - k5*y[9]                      #y9=Np
  dy10=k3*y[8] + k5*y[9] - k6*y[10]          #y10=Nm
  dy11=k6*y[10] - (k7+k8)*y[11]              #y11=C
  dy12=k7*y[11] - k10*y[12]                   
  
  preExistingTmp <- c(dy7,dy8,dy9,dy10,dy11,dy12)
  
  list(c(nascentTmp,preExistingTmp))
}

dataGenerationFull <- function(exampleRates # Set of rate.
                               ,TauFractions # Time points with cellular fractionation.
                               ,TauPoly # Time points with polysomal profiling.
                               ,TauTotal) # Time points with total RNA profiling.
{
  k1=exampleRates[["k1"]] # Synthesis.
  k2=exampleRates[["k2"]] # Co-transcriptional splicing.
  k3=exampleRates[["k3"]] # Detachment of spliced chromatin associated RNA.
  k4=exampleRates[["k4"]] # Detachment of unspliced chromatin associated RNA.
  k5=exampleRates[["k5"]] # Post-transcriptional splicing.
  k6=exampleRates[["k6"]] # Export.
  k7=exampleRates[["k7"]] # Translation.
  k8=exampleRates[["k8"]] # Cytoplasmic decay.
  k10=exampleRates[["k10"]] #Polysomal decay.
  
  # Initial conditions.
  yini <- c(y1=0,y2=0,y3=0,y4=0,y5=0,y6=0
            ,y7= k1/(k2+k4)
            ,y8= k1*k2/(k3*(k2+k4))
            ,y9= k4*k1/(k5*(k2+k4))
            ,y10= k1/k6
            ,y11= k1/(k7+k8)
            ,y12= k7*k1/(k10*(k7+k8)))
  
  # All time-points to be simulated
  TauGlobals <- unique(c(TauFractions,TauPoly,TauTotal))
  
  # Numerical solution of the ODE model.
  exampleData <- ode(times=TauGlobals # Time points.
                     ,y=yini # Initial conditions.
                     ,func=fullModelSimulation # Model.
                     ,parms=exampleRates) # Rates.
  
  # Removal of time points (first column).
  exampleData <- c(t(exampleData[,-1]))
  
  # Set of colnames.
  names(exampleData) <- c(sapply(TauGlobals,function(i)paste0(c("chpn","chmn","npn","nmn","cyn","pyn"
                                                                ,"chpp","chmp","npp","nmp","cyp","pyp"),i)))
  
  # Steady state solution (initial conditions).
  exampleData <- c(tail(yini,(length(yini)/2)),exampleData)
  names(exampleData)[1:(length(yini)/2)] <- c("chpt","chmt","npt","nmt","cyt","pyt")
  
  # Cytoplasmic RNA is re-defined as the sum of polysomal RNA and cytoplasmic RNA from the model
  # because this is closer to the experimental data.
  # Here we could add a flag to avoid this step if this assumption is not true for a specific experiment.
  exampleData[sort(names(exampleData[grep("^cyn",names(exampleData))]))] <- exampleData[sort(names(exampleData[grep("^cyn",names(exampleData))]))]+exampleData[sort(names(exampleData[grep("^pyn",names(exampleData))]))] 
  exampleData[sort(names(exampleData[grep("^cyp",names(exampleData))]))]  <- exampleData[sort(names(exampleData[grep("^cyp",names(exampleData))]))]+exampleData[sort(names(exampleData[grep("^pyp",names(exampleData))]))] 
  exampleData["cyt"] <- exampleData["cyt"]+exampleData["pyt"]
  
  exampleData <- c(exampleData[grep("t$",names(exampleData))] # We keep the total fractions
                   ,unlist(lapply(c(TauFractions,TauTotal),function(i)exampleData[grepl(paste0(i,"$"),names(exampleData))
                                                                                  &!grepl(paste0("^pyn",i,"$"),names(exampleData))
                                                                                  &!grepl(paste0("^pyp",i,"$"),names(exampleData))]))
                   ,unlist(lapply(TauPoly,function(i)exampleData[grepl(paste0("^pyn",i,"$"),names(exampleData))
                                                                 |grepl(paste0("^pyp",i,"$"),names(exampleData))])))
  
  # If some fractionation time-points have been profiled only for total RNA we compute the total nascent and total
  # pre-existing, and we remove the corresponding fractions.
  if(!is.null(TauTotal))
  {
    exampleData <- c(exampleData
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("n",i,"$"),names(exampleData))
                                                                  &!grepl(paste0("^pyn",i,"$"),names(exampleData))])}) # For certain time-points we just have nascent total RNA.
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("p",i,"$"),names(exampleData))
                                                                  &!grepl(paste0("^pyp",i,"$"),names(exampleData))])})) # For certain time-points we just have pre-existing total RNA.
    
    names(exampleData)[names(exampleData)==""] <- c(paste0("n",TauTotal),paste0("p",TauTotal))
    
    # Removal of fractions for time-points profiled without cellular fractionation.
    exampleData <- exampleData[!(names(exampleData)%in%c(sapply(TauTotal,function(i)paste0(c("chpn","chmn","npn","nmn","cyn","chpp","chmp","npp","nmp","cyp"),i))))]
    
    # Removal of eventual duplicated conditions
    exampleData <- exampleData[!duplicated(names(exampleData))]
  }
  
  exampleData <- exampleData[!duplicated(names(exampleData))]
  exampleData <- exampleData[sort(names(exampleData))]
}

## Full model - No Poly
modelSimulationNoPoly <- function(t
                                  ,y
                                  ,parms)
{
  k1=parms[["k1"]] # Synthesis
  k2=parms[["k2"]]# Co-transcriptional splicing
  k3=parms[["k3"]] # Detachment of spliced chromatin
  k4=parms[["k4"]] # Detachment of unspliced chromatin
  k5=parms[["k5"]] # Post-transcriptional splicing
  k6=parms[["k6"]] # Export
  k8=parms[["k8"]] # Decay
  
  # System of equations for NASCENT RNA
  dy1=k1 - (k2+k4)*y[1]                    #y1=Chp
  dy2=k2*y[1] - k3*y[2]                   #y2=Chm
  dy3=k4*y[1] - k5*y[3]                     #y3=Np
  dy4=k3*y[2] + k5*y[3] - k6*y[4]          #y4=Nm
  dy5=k6*y[4] - k8*y[5]                     #y5=C
  
  nascentTmp <- c(dy1,dy2,dy3,dy4,dy5)
  
  # System of equations for PRE-EXISTING RNA
  dy6=0 - (k2+k4)*y[6]                          #y7=Chp
  dy7=k2*y[6] - k3*y[7]                        #y8=Chm
  dy8=k4*y[6] - k5*y[8]                          #y9=Np
  dy9=k3*y[7] + k5*y[8] - k6*y[9]             #y10=Nm
  dy10=k6*y[9] - k8*y[10]                       #y11=C
  
  preExistingTmp <- c(dy6,dy7,dy8,dy9,dy10)
  
  list(c(nascentTmp,preExistingTmp))
}

dataGenerationNoPoly <- function(exampleRates # Set of rate.
                                 ,TauFractions # Time points with cellular fractionation.
                                 ,TauPoly # Time points with polysomal profiling (kept only to use the same simulateData function).
                                 ,TauTotal) # Time points with polysomal profiling.
{
  k1=exampleRates[["k1"]] # Synthesis.
  k2=exampleRates[["k2"]] # Co-transcriptional splicing.
  k3=exampleRates[["k3"]] # Detachment of spliced chromatin associated RNA.
  k4=exampleRates[["k4"]] # Detachment of unspliced chromatin associated RNA.
  k5=exampleRates[["k5"]] # Post-transcriptional splicing.
  k6=exampleRates[["k6"]] # Export.
  k8=exampleRates[["k8"]] # Decay.
  
  # Initial conditions.
  yini <- c(y1=0,y2=0,y3=0,y4=0,y5=0
            ,y6= k1/(k2+k4)
            ,y7= k1*k2/(k3*(k2+k4))
            ,y8= k4*k1/(k5*(k2+k4))
            ,y9= k1/k6
            ,y10= k1/k8)
  
  # All time-points to be simulated
  TauGlobals <- unique(c(TauFractions,TauTotal))
  
  # Numerical solution of the ODE model.
  exampleData <- ode(times=TauGlobals # Time points.
                     ,y=yini # Initial conditions.
                     ,func=modelSimulationNoPoly # Model.
                     ,parms=exampleRates) # Rates.
  
  # Removal of time points (first column).
  exampleData <- c(t(exampleData[,-1]))
  
  # Set of colnames.
  names(exampleData) <- c(sapply(TauGlobals,function(i)paste0(c("chpn","chmn","npn","nmn","cyn"
                                                                ,"chpp","chmp","npp","nmp","cyp"),i)))
  
  # Steady state solution
  exampleData <- c(tail(yini,(length(yini)/2)),exampleData)
  names(exampleData)[1:(length(yini)/2)] <- c("chpt","chmt","npt","nmt","cyt")
  
  # If some fractionation time-points have been profiled only for total RNA we compute the total nascent and total
  # pre-existing, and we remove the corresponding fractions.
  if(!is.null(TauTotal))
  {
    exampleData <- c(exampleData
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("n",i,"$"),names(exampleData))])}) # For certain time-points we just have nascent total RNA.
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("p",i,"$"),names(exampleData))])})) # For certain time-points we just have pre-existing total RNA.
    
    names(exampleData)[names(exampleData)==""] <- c(paste0("n",TauTotal),paste0("p",TauTotal))
    
    # Removal of fractions for time-points profiled without cellular fractionation.
    exampleData <- exampleData[!(names(exampleData)%in%c(sapply(TauTotal,function(i)paste0(c("chpn","chmn","npn","nmn","cyn","chpp","chmp","npp","nmp","cyp"),i))))]
    
    # Removal of eventual duplicated conditions
    exampleData <- exampleData[!duplicated(names(exampleData))]
  }
  
  exampleData <- exampleData[!duplicated(names(exampleData))]
  exampleData <- exampleData[sort(names(exampleData))]
}

## Chr. P Only - Yes Poly
modelSimulationNoNucP <- function(t
                                  ,y
                                  ,parms)
{
  k1=parms[["k1"]] # Synthesis
  k2=parms[["k2"]] # Co-transcriptional splicing
  k3=parms[["k3"]] # Detachment of unspliced chromatin
  k6=parms[["k6"]] # Export
  k7=parms[["k7"]] #Translation
  k8=parms[["k8"]] #Cytoplasmic decay
  k10=parms[["k10"]] #Polysomal decay
  
  # System of equations for NASCENT RNA
  dy1=k1 - k2*y[1]                          #y1=Chp
  dy2=k2*y[1] - k3*y[2]                     #y2=Chm
  dy3=k3*y[2] - k6*y[3]                     #y3=N
  dy4=k6*y[3] - (k7+k8)*y[4]                #y4=C
  dy5=k7*y[4] - k10*y[5]                     #y5=P
  
  nascentTmp <- c(dy1,dy2,dy3,dy4,dy5)
  
  # System of equations for PRE-EXISTING RNA
  dy6=0 - k2*y[6]                           #y6=Chp
  dy7=k2*y[6] - k3*y[7]                     #y7=Chm
  dy8=k3*y[7] - k6*y[8]                     #y8=Nm
  dy9=k6*y[8] - (k7+k8)*y[9]                #y9=C
  dy10=k7*y[9] - k10*y[10]                   #y10=P
  
  
  preExistingTmp <- c(dy6,dy7,dy8,dy9,dy10)
  
  list(c(nascentTmp,preExistingTmp))
}

dataGenerationNoNucP <- function(exampleRates # Set of rate.
                                 ,TauFractions # Time points with cellular fractionation.
                                 ,TauPoly # Time points with polysomal profiling.
                                 ,TauTotal) # Time points with total RNA profiling.
{
  k1=exampleRates[["k1"]] # Synthesis.
  k2=exampleRates[["k2"]] # Co-transcriptional splicing.
  k3=exampleRates[["k3"]] # Detachment of spliced chromatin associated RNA.
  k6=exampleRates[["k6"]] # Export.
  k7=exampleRates[["k7"]] # Translation.
  k8=exampleRates[["k8"]] # Cytoplasmic decay.
  k10=exampleRates[["k10"]] #Polysomal decay.
  
  # Initial conditions.
  yini <- c(y1=0,y2=0,y3=0,y4=0,y5=0
            ,y6= k1/k2
            ,y7= k1/k3
            ,y8= k1/k6
            ,y9= k1/(k7+k8)
            ,y10=k1*k7/(k10*(k7+k8)))
  
  # All time-points to be simulated
  TauGlobals <- unique(c(TauFractions,TauPoly,TauTotal))
  
  # Numerical solution of the ODE model.
  exampleData <- ode(times=TauGlobals # Time points.
                     ,y=yini # Initial conditions.
                     ,func=modelSimulationNoNucP # Model.
                     ,parms=exampleRates) # Rates.
  
  # Removal of time points (first column).
  exampleData <- c(t(exampleData[,-1]))
  
  # Set of colnames.
  names(exampleData) <- c(sapply(TauGlobals,function(i)paste0(c("chpn","chmn","nn","cyn","pyn"
                                                                ,"chpp","chmp","np","cyp","pyp"),i)))
  
  # Steady state solution (initial conditions).
  exampleData <- c(tail(yini,(length(yini)/2)),exampleData)
  names(exampleData)[1:(length(yini)/2)] <- c("chpt","chmt","nt","cyt","pyt")
  
  # Cytoplasmic RNA is re-defined as the sum of polysomal RNA and cytoplasmic RNA from the model
  # because this is closer to the experimental data.
  # Here we could add a flag to avoid this step if this assumption is not true for a specific experiment.
  exampleData[sort(names(exampleData[grep("^cyn",names(exampleData))]))] <- exampleData[sort(names(exampleData[grep("^cyn",names(exampleData))]))]+exampleData[sort(names(exampleData[grep("^pyn",names(exampleData))]))] 
  exampleData[sort(names(exampleData[grep("^cyp",names(exampleData))]))]  <- exampleData[sort(names(exampleData[grep("^cyp",names(exampleData))]))]+exampleData[sort(names(exampleData[grep("^pyp",names(exampleData))]))] 
  exampleData["cyt"] <- exampleData["cyt"]+exampleData["pyt"]
  
  exampleData <- c(exampleData[grep("t$",names(exampleData))] # We keep the total fractions
                   ,unlist(lapply(c(TauFractions,TauTotal),function(i)exampleData[grepl(paste0(i,"$"),names(exampleData))
                                                                                  &!grepl(paste0("^pyn",i,"$"),names(exampleData))
                                                                                  &!grepl(paste0("^pyp",i,"$"),names(exampleData))]))
                   ,unlist(lapply(TauPoly,function(i)exampleData[grepl(paste0("^pyn",i,"$"),names(exampleData))
                                                                 |grepl(paste0("^pyp",i,"$"),names(exampleData))])))
  
  # If some fractionation time-points have been profiled only for total RNA we compute the total nascent and total
  # pre-existing, and we remove the corresponding fractions.
  if(!is.null(TauTotal))
  {
    exampleData <- c(exampleData
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("n",i,"$"),names(exampleData))
                                                                  &!grepl(paste0("^pyn",i,"$"),names(exampleData))])}) # For certain time-points we just have nascent total RNA.
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("p",i,"$"),names(exampleData))
                                                                  &!grepl(paste0("^pyp",i,"$"),names(exampleData))])})) # For certain time-points we just have pre-existing total RNA.
    
    names(exampleData)[names(exampleData)==""] <- c(paste0("n",TauTotal),paste0("p",TauTotal))
    
    # Removal of fractions for time-points profiled without cellular fractionation.
    exampleData <- exampleData[!(names(exampleData)%in%c(sapply(TauTotal,function(i)paste0(c("chpn","chmn","nn","cyn","chpp","chmp","np","cyp"),i))))]
    
    # Removal of eventual duplicated conditions
    exampleData <- exampleData[!duplicated(names(exampleData))]
  }
  
  exampleData <- exampleData[!duplicated(names(exampleData))]
  exampleData <- exampleData[sort(names(exampleData))]
}

## Chr. P Only - No Poly
modelSimulationNoPolyNoNucP <- function(t
                                        ,y
                                        ,parms)
{
  k1=parms[["k1"]] # Synthesis
  k2=parms[["k2"]] # Co-transcriptional splicing
  k3=parms[["k3"]] # Detachment of unspliced chromatin
  k6=parms[["k6"]] # Export
  k8=parms[["k8"]] # Decay
  
  # System of equations for NASCENT RNA
  dy1=k1 - k2*y[1]                          #y1=Chp
  dy2=k2*y[1] - k3*y[2]                     #y2=Chm
  dy3=k3*y[2] - k6*y[3]                     #y3=N
  dy4=k6*y[3] - k8*y[4]                     #y4=C
  
  nascentTmp <- c(dy1,dy2,dy3,dy4)
  
  # System of equations for PRE-EXISTING RNA
  dy5=0 - k2*y[5]                               #y5=Chp
  dy6=k2*y[5] - k3*y[6]                         #y6=Chm
  dy7=k3*y[6] - k6*y[7]                         #y7=Np
  dy8=k6*y[7] - k8*y[8]                         #y8=C
  
  
  preExistingTmp <- c(dy5,dy6,dy7,dy8)
  
  list(c(nascentTmp,preExistingTmp))
}

dataGenerationNoPolyNoNucP <- function(exampleRates # Set of rate.
                                       ,TauFractions # Time points with cellular fractionation.
                                       ,TauPoly # Time points with polysomal profiling (kept only to use the same simulateData function).
                                       ,TauTotal) # Time points with total RNA profiling.
{
  k1=exampleRates[["k1"]] # Synthesis.
  k2=exampleRates[["k2"]] # Co-transcriptional splicing.
  k3=exampleRates[["k3"]] # Detachment of spliced chromatin associated RNA.
  k6=exampleRates[["k6"]] # Export.
  k8=exampleRates[["k8"]] # Decay.
  
  # Initial conditions.
  yini <- c(y1=0,y2=0,y3=0,y4=0
            ,y5= k1/k2
            ,y6= k1/k3
            ,y7= k1/k6
            ,y8= k1/k8)
  
  # All time-points to be simulated
  TauGlobals <- unique(c(TauFractions,TauTotal))
  
  # Numerical solution of the ODE model.
  exampleData <- ode(times=TauGlobals # Time points.
                     ,y=yini # Initial conditions.
                     ,func=modelSimulationNoPolyNoNucP # Model.
                     ,parms=exampleRates) # Rates.
  
  # Removal of time points (first column).
  exampleData <- c(t(exampleData[,-1]))
  
  # Set of colnames.
  names(exampleData) <- c(sapply(TauGlobals,function(i)paste0(c("chpn","chmn","nn","cyn"
                                                                ,"chpp","chmp","np","cyp"),i)))
  
  # Steady state solution
  exampleData <- c(tail(yini,(length(yini)/2)),exampleData)
  names(exampleData)[1:(length(yini)/2)] <- c("chpt","chmt","nt","cyt")
  
  # If some fractionation time-points have been profiled only for total RNA we compute the total nascent and total
  # pre-existing, and we remove the corresponding fractions.
  if(!is.null(TauTotal))
  {
    exampleData <- c(exampleData
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("n",i,"$"),names(exampleData))])}) # For certain time-points we just have nascent total RNA.
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("p",i,"$"),names(exampleData))])})) # For certain time-points we just have pre-existing total RNA.
    
    names(exampleData)[names(exampleData)==""] <- c(paste0("n",TauTotal),paste0("p",TauTotal))
    
    # Removal of fractions for time-points profiled without cellular fractionation.
    exampleData <- exampleData[!(names(exampleData)%in%c(sapply(TauTotal,function(i)paste0(c("chpn","chmn","nn","cyn","chpp","chmp","np","cyp"),i))))]
    
    # Removal of eventual duplicated conditions
    exampleData <- exampleData[!duplicated(names(exampleData))]
  }
  
  exampleData <- exampleData[!duplicated(names(exampleData))]
  exampleData <- exampleData[sort(names(exampleData))]
}

## No Premature - Yes Poly.
modelSimulationNoP <- function(t
                               ,y
                               ,parms)
{
  k1=parms[["k1"]] # Synthesis
  k3=parms[["k3"]] # Detachment of spliced chromatin
  k6=parms[["k6"]] # Export
  k7=parms[["k7"]] # Translation
  k8=parms[["k8"]] # Cytoplasmic degradation
  k10=parms[["k10"]] # Polysomal degradation
  
  # System of equations for NASCENT RNA
  dy1=k1 - k3*y[1]                               #y1=Chm
  dy2=k3*y[1] - k6*y[2]                          #y2=Nm
  dy3=k6*y[2] - (k7+k8)*y[3]                     #y3=C
  dy4=k7*y[3] - k10*y[4]                          #y4=P
  
  nascentTmp <- c(dy1,dy2,dy3,dy4)
  
  # System of equations for PRE-EXISTING RNA
  dy5=0 - k3*y[5]                                #y5=Chm
  dy6=k3*y[5] - k6*y[6]                          #y6=Nm
  dy7=k6*y[6] - (k7+k8)*y[7]                     #y7=C
  dy8=k7*y[7] - k10*y[8]                          #y8=P
  
  preExistingTmp <- c(dy5,dy6,dy7,dy8)
  
  list(c(nascentTmp,preExistingTmp))
}

dataGenerationNoP <- function(exampleRates # Set of rate.
                              ,TauFractions # Time points with cellular fractionation.
                              ,TauPoly # Time points with polysomal profiling.
                              ,TauTotal) # Time points with total RNA profiling.
{
  k1=exampleRates[["k1"]] # Synthesis.
  k3=exampleRates[["k3"]] # Processing.
  k6=exampleRates[["k6"]] # Export.
  k7=exampleRates[["k7"]] # Translation.
  k8=exampleRates[["k8"]] # Cytoplasmic decay.
  k10=exampleRates[["k10"]] # Polysomal decay.
  
  # Initial conditions.
  yini <- c(y1=0,y2=0,y3=0,y4=0
            ,y5= k1/k3
            ,y6= k1/k6
            ,y7= k1/(k7+k8)
            ,y8=k7*k1/(k10*(k7+k8)))
  
  # All time-points to be simulated
  TauGlobals <- unique(c(TauFractions,TauPoly,TauTotal))
  
  # Numerical solution of the ODE model.
  exampleData <- ode(times=TauGlobals # Time points.
                     ,y=yini # Initial conditions.
                     ,func=modelSimulationNoP # Model.
                     ,parms=exampleRates) # Rates.
  
  # Removal of time points (first column).
  exampleData <- c(t(exampleData[,-1]))
  
  # Set of colnames.
  names(exampleData) <- c(sapply(TauGlobals,function(i)paste0(c("chn","nn","cyn","pyn"
                                                                ,"chp","np","cyp","pyp"),i)))
  
  # Steady state solution (initial conditions).
  exampleData <- c(tail(yini,(length(yini)/2)),exampleData)
  names(exampleData)[1:(length(yini)/2)] <- c("cht","nt","cyt","pyt")
  
  # Cytoplasmic RNA is re-defined as the sum of polysomal RNA and cytoplasmic RNA from the model
  # because this is closer to the experimental data.
  # Here we could add a flag to avoid this step if this assumption is not true for a specific experiment.
  exampleData[sort(names(exampleData[grep("^cyn",names(exampleData))]))] <- exampleData[sort(names(exampleData[grep("^cyn",names(exampleData))]))]+exampleData[sort(names(exampleData[grep("^pyn",names(exampleData))]))] 
  exampleData[sort(names(exampleData[grep("^cyp",names(exampleData))]))]  <- exampleData[sort(names(exampleData[grep("^cyp",names(exampleData))]))]+exampleData[sort(names(exampleData[grep("^pyp",names(exampleData))]))] 
  exampleData["cyt"] <- exampleData["cyt"]+exampleData["pyt"]
  
  exampleData <- c(exampleData[grep("t$",names(exampleData))] # We keep the total fractions
                   ,unlist(lapply(c(TauFractions,TauTotal),function(i)exampleData[grepl(paste0(i,"$"),names(exampleData))
                                                                                  &!grepl(paste0("^pyn",i,"$"),names(exampleData))
                                                                                  &!grepl(paste0("^pyp",i,"$"),names(exampleData))]))
                   ,unlist(lapply(TauPoly,function(i)exampleData[grepl(paste0("^pyn",i,"$"),names(exampleData))
                                                                 |grepl(paste0("^pyp",i,"$"),names(exampleData))])))
  
  # If some fractionation time-points have been profiled only for total RNA we compute the total nascent and total
  # pre-existing, and we remove the corresponding fractions.
  if(!is.null(TauTotal))
  {
    exampleData <- c(exampleData
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("n",i,"$"),names(exampleData))
                                                                  &!grepl(paste0("^pyn",i,"$"),names(exampleData))])}) # For certain time-points we just have nascent total RNA.
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("p",i,"$"),names(exampleData))
                                                                  &!grepl(paste0("^pyp",i,"$"),names(exampleData))])})) # For certain time-points we just have pre-existing total RNA.
    
    names(exampleData)[names(exampleData)==""] <- c(paste0("n",TauTotal),paste0("p",TauTotal))
    
    # Removal of fractions for time-points profiled without cellular fractionation.
    exampleData <- exampleData[!(names(exampleData)%in%c(sapply(TauTotal,function(i)paste0(c("chn","nn","cyn","chp","np","cyp"),i))))]
    
    # Removal of eventual duplicated conditions
    exampleData <- exampleData[!duplicated(names(exampleData))]
  }
  
  exampleData <- exampleData[!duplicated(names(exampleData))]
  exampleData <- exampleData[sort(names(exampleData))]
}

## No Premature - No Poly.
modelSimulationNoPolyNoP <- function(t
                                     ,y
                                     ,parms)
{
  k1=parms[["k1"]] # Synthesis.
  k3=parms[["k3"]] # Processing.
  k6=parms[["k6"]] # Export.
  k8=parms[["k8"]] # Degradation.
  
  # System of equations for NASCENT RNA
  dy1=k1 - k3*y[1]                               #y1=Chm
  dy2=k3*y[1] - k6*y[2]                          #y2=Nm
  dy3=k6*y[2] - k8*y[3]                          #y3=C
  
  nascentTmp <- c(dy1,dy2,dy3)
  
  # System of equations for PRE-EXISTING RNA
  dy4=0 - k3*y[4]                                #y5=Chm
  dy5=k3*y[4] - k6*y[5]                          #y8=Nm
  dy6=k6*y[5] - k8*y[6]                          #y9=C
  
  preExistingTmp <- c(dy4,dy5,dy6)
  
  list(c(nascentTmp,preExistingTmp))
}

dataGenerationNoPolyNoP <- function(exampleRates # Set of rate.
                                    ,TauFractions # Time points with cellular fractionation.
                                    ,TauPoly # Time points with polysomal profiling (kept only to use the same simulateData function).
                                    ,TauTotal) # Time points with total RNA profiling.
{
  k1=exampleRates[["k1"]] # Synthesis
  k3=exampleRates[["k3"]] # Processing
  k6=exampleRates[["k6"]] # Export
  k8=exampleRates[["k8"]] # Degradation
  
  # Initial conditions.
  yini <- c(y1=0,y2=0,y3=0
            ,y4= k1/k3
            ,y5= k1/k6
            ,y6= k1/k8)
  
  # All time-points to be simulated
  TauGlobals <- unique(c(TauFractions,TauTotal))
  
  # Numerical solution of the ODE model.
  exampleData <- ode(times=TauGlobals # Time points.
                     ,y=yini # Initial conditions.
                     ,func=modelSimulationNoPolyNoP # Model.
                     ,parms=exampleRates) # Rates.
  
  # Removal of time points (first column).
  exampleData <- c(t(exampleData[,-1]))
  
  # Set of colnames.
  names(exampleData) <- c(sapply(TauGlobals,function(i)paste0(c("chn","nn","cyn"
                                                                ,"chp","np","cyp"),i)))
  
  # Steady state solution
  exampleData <- c(tail(yini,(length(yini)/2)),exampleData)
  names(exampleData)[1:(length(yini)/2)] <- c("cht","nt","cyt")
  
  # If some fractionation time-points have been profiled only for total RNA we compute the total nascent and total
  # pre-existing, and we remove the corresponding fractions.
  if(!is.null(TauTotal))
  {
    exampleData <- c(exampleData
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("n",i,"$"),names(exampleData))])}) # For certain time-points we just have nascent total RNA.
                     ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("p",i,"$"),names(exampleData))])})) # For certain time-points we just have pre-existing total RNA.
    
    names(exampleData)[names(exampleData)==""] <- c(paste0("n",TauTotal),paste0("p",TauTotal))
    
    # Removal of fractions for time-points profiled without cellular fractionation.
    exampleData <- exampleData[!(names(exampleData)%in%c(sapply(TauTotal,function(i)paste0(c("chn","nn","cyn","chp","np","cyp"),i))))]
    
    # Removal of eventual duplicated conditions
    exampleData <- exampleData[!duplicated(names(exampleData))]
  }
  
  exampleData <- exampleData[!duplicated(names(exampleData))]
  exampleData <- exampleData[sort(names(exampleData))]
}

### Data simulations
##  - Function to mimic replicates,
##  - Data simulation function.

dataAverageFunction <- function(v
                                ,nRep
                                ,CVs)
{
  species <- c("ch","^n","cy","py")
  
  generated_gaussian <- lapply(seq_along(species), function(x){
    V <- v[grep(species[x],names(v))]
    qwe <- lapply(names(V), function(i){
      rnorm(nRep,m=V[i],sd=abs(CVs[[x]][i]*V[i]))})
    names(qwe) <- names(V)
    qwe
  })
  generated_gaussian <- t(do.call(rbind,unlist(generated_gaussian,recursive=FALSE)))
  
  #generated_gaussian <- sapply(v, function(x){rnorm(nRep,m=x,sd=abs(CV*x))})
  #these rows handle the case n=1
  if(is.null(nrow(generated_gaussian)))
  {
    namesTmp <- names(generated_gaussian)
    generated_gaussian <- matrix(generated_gaussian,nrow=1)
    colnames(generated_gaussian) <- namesTmp
  }
  if(nRep!=1){
    medie <- apply(generated_gaussian, 2, function(x){mean(x)})
    deviazioni <- apply(generated_gaussian, 2, function(x){sd(x)})
  }else{
    medie <- generated_gaussian[1,]
    deviazioni <- medie*unlist(unname(CVs))[sort(names(unlist(unname(CVs))))]
  }
  
  list(medie,deviazioni,generated_gaussian)
}

simulateData <- function(exampleRates # Mean rates to sample.
                         ,TauFractions # Time points with cellular fractionation.
                         ,TauPoly # Time points with polysomal profiling.
                         ,TauTotal # Time points with total nascent and pre-existing RNA profiling.
                         ,noise # TRUE to add noise to data.
                         ,CVs # Variation Coefficients for noise.
                         ,Reps # Number of replicates.
                         ,nGenes # Number of genes to be simulated.
                         ,seed # Seed for reproducibility.
                         ,ZeroThresh # Minimum expression value.
                         ,MultFact # Number of nGenes to be modeled.
                         ,cpus=1) # Number of CPUs.
{
  if(length(intersect(TauFractions,TauTotal))>0)
  {
    print("Define independent Total and Fractions RNA temporal profiles.")
    stop()
  }
  
  # All time-points to be simulated
  TauGlobals <- unique(c(TauFractions,TauTotal,TauPoly))
  
  # Selection of the configuration to be simulated according to the provided example rates.
  # This means definition of the modelSimulation function (ode implementation) and of dataGeneration function (ode solution).
  
  if("k9"%in%names(exampleRates))
  { #Full Model with nuclear mature RNA decay
    modelSimulation <- fullModelSimulationNDM
    dataGeneration <- dataGenerationFullNDM
    print("Full Model with nuclear mature RNA decay")
  }else if("k11"%in%names(exampleRates))
  { #Full Model with nuclear premature RNA decay
    modelSimulation <- fullModelSimulationNDP
    dataGeneration <- dataGenerationFullNDP
    print("Full Model with nuclear premature RNA decay")
  }else if(all(c("k1","k2","k3","k4","k5","k6","k7","k8","k10")%in%names(exampleRates)))
  { #Full Model
    modelSimulation <- fullModelSimulation
    dataGeneration <- dataGenerationFull
    print("Full Model")
  }else if(all(c("k1","k2","k3","k4","k5","k6","k8")%in%names(exampleRates)))
  { #Full Model without translation
    modelSimulation <- modelSimulationNoPoly
    dataGeneration <- dataGenerationNoPoly
    print("Full Model without translation")
  }else if(all(c("k1","k2","k3","k6","k7","k8","k10")%in%names(exampleRates)))
  { #No Nucleoplasmic Premature - yes tranlation
    modelSimulation <- modelSimulationNoNucP
    dataGeneration <- dataGenerationNoNucP
    print("No Nucleoplasmic Premature - yes tranlation")
  }else if(all(c("k1","k2","k3","k6","k8")%in%names(exampleRates)))
  { #No Nucleoplasmic Premature - no tranlation
    modelSimulation <- modelSimulationNoPolyNoNucP
    dataGeneration <- dataGenerationNoPolyNoNucP
    print("No Nucleoplasmic Premature - no tranlation")
  }else if(all(c("k1","k3","k6","k7","k8","k10")%in%names(exampleRates)))
  { #No premature - yes translation
    modelSimulation <- modelSimulationNoP
    dataGeneration <- dataGenerationNoP
    print("No premature - yes translation")
  }else if(all(c("k1","k3","k6","k8")%in%names(exampleRates)))
  { #No premature - no translation
    modelSimulation <- modelSimulationNoPolyNoP
    dataGeneration <- dataGenerationNoPolyNoP
    print("No premature - no translation")
  }else{
    print("Error in configuration!")
  }
  
  # Definition of a specific name for the configuration for output purposes.
  elements <- c("Noise","Reps","TauFractions","TauPoly","TauTotal",names(exampleRates))
  values <- c(noise,Reps,paste0(TauFractions,collapse="_"),paste0(TauPoly,collapse = "_"),paste0(TauTotal,collapse = "_"),unname(exampleRates))
  nameTmp <- character(0)
  for (i in 1:length(elements)){
    name <- paste0(c(elements[i],values[i]),collapse = "_")
    nameTmp <- paste0(c(nameTmp,name),collapse ="_")
  }
  
  # Sampling of the rates.
  set.seed(seed)
  multipleExampleRates <- mcsapply(exampleRates,function(x)rnorm(nGenes*MultFact,mean=x,sd=5*x),mc.cores=cpus)
  
  #Selection of positive rate-sets.
  multipleExampleRates <- multipleExampleRates[rowSums(multipleExampleRates<0)==0,]
  
  # Inclusion of the original rates set.
  multipleExampleRates <- rbind(exampleRates,multipleExampleRates)
  rownames(multipleExampleRates) <- 1:nrow(multipleExampleRates)
  
  # CVs names check
  exampleDataTmp <- dataGeneration(exampleRates=multipleExampleRates[1,] # Rates.
                                   ,TauFractions=TauFractions # Time points with cellular fractionation.
                                   ,TauPoly=TauPoly # Time points with polysomal profiling.
                                   ,TauTotal=TauTotal) # Time points with total nascent and pre-existing RNA profiling.
  
  if(!(all(names(exampleDataTmp)%in%names(unlist(unname(CVs))))))
  {
    print("Missing CVs.")
  }else{
    CVs <- lapply(CVs,function(i)i[names(i)%in%names(exampleDataTmp)])
  }
  
  # Data generation and addition of noise at the expression level.
  set.seed(seed)
  Obj <- mclapply(1:nrow(multipleExampleRates),function(i)
  {
    par <- multipleExampleRates[i,]
    # Data generation.
    exampleDataTmp <- dataGeneration(exampleRates=par # Rates.
                                     ,TauFractions=TauFractions # Time points with cellular fractionation.
                                     ,TauPoly=TauPoly # Time points with polysomal profiling.
                                     ,TauTotal=TauTotal) # Time points with total nascent and pre-existing RNA profiling.
    
    if(noise)
    {
      # Noise addition and averages computation.
      ListDataTmp <- dataAverageFunction(v=exampleDataTmp # Expression levels.
                                         ,nRep=Reps # Number of replicates.
                                         ,CVs=CVs) # Variation Coefficient.
      return(ListDataTmp)
    }else{
      return(exampleDataTmp)
    }
  },mc.cores=cpus)
  names(Obj) <- rownames(multipleExampleRates)
  
  Obj <- Obj[mcsapply(Obj,function(i){all(apply(i[[3]],2,min)[!grepl("n0$",colnames(i[[3]]))]>1e-10)},mc.cores=cpus)]
  
  # Extraction of the quantities of interest.
  exampleData <- t(sapply(Obj,"[[",1))           #deviations
  DevDataTmp <- t(sapply(Obj,"[[",2))           #deviations
  
  # Conditions with low coverage are set to the minimum coverage value.
  exampleData[abs(exampleData)<ZeroThresh] <- ZeroThresh
  
  # Selection of suitable examples without mean negative expression levels.
  # exampleData[which(exampleData<0)] <- NaN
  # exampleData <- exampleData[apply(exampleData,1,function(r)all(is.finite(r))),]
  
  # Selection of at last nGenes examples.
  if(nrow(exampleData)<nGenes){print(paste0("Warning: less than ",nGenes," genes provided."))}
  exampleData <- exampleData[1:min(nrow(exampleData),nGenes),]
  multipleExampleRates <- multipleExampleRates[rownames(exampleData),]
  DevDataTmp <- DevDataTmp[rownames(exampleData),]
  
  singleReplicateExpressionData <- lapply(Obj[rownames(exampleData)],"[[",3)
  
  rownames(exampleData) <- 1:nrow(exampleData)
  rownames(DevDataTmp) <- 1:nrow(DevDataTmp)
  rownames(multipleExampleRates) <- 1:nrow(multipleExampleRates)
  names(singleReplicateExpressionData) <- seq_along(singleReplicateExpressionData)
  
  singleReplicateExpressionData <- lapply(1:nrow(singleReplicateExpressionData[[1]]),function(i)sapply(singleReplicateExpressionData,function(j)j[i,]))
  singleReplicateExpressionData <- lapply(singleReplicateExpressionData,function(i){i[i<ZeroThresh]<-ZeroThresh;t(i)})
  
  # Output.
  return(list("exampleData"=exampleData,"DevDataTmp" = DevDataTmp
              ,"exampleRates"=multipleExampleRates,"name"=nameTmp
              ,"TauGlobals"=TauGlobals,"TauFractions"=TauFractions
              ,"TauPoly"=TauPoly,"TauTotal"=TauTotal
              ,"modelSimulation"=modelSimulation,"dataGeneration"=dataGeneration,"singleReplicateExpressionData"=singleReplicateExpressionData))
}

### Plot functions
RatesDistributions <- function(object
                               ,width
                               ,height
                               ,obj_name
                               ,xlimTmp=c(NA,NA)
                               ,ratesOrder=c("k1","k2","k3","k4","k5","k11","k9","k6","k8","k7","k10"))
{
  #Creating Palette
  myPalette <- RColorBrewer::brewer.pal(9,"Set3")
  
  inferedRates <- as.data.frame(object$inferedRates)
  ActualRates <- colnames(inferedRates) 
  # AllRates <- c("k1","k2","k3","k4","k5","k6","k7","k8")
  AllRates <- colnames(inferedRates)
  MissingRates <- setdiff(AllRates,ActualRates)
  
  RateNames <- c("Synthesis"
                 ,"Co-trascriptional\nsplicing"
                 ,"Detachment of\nmature chromatin"
                 ,"Detachment of\npremature chromatin"
                 ,"Post-transcriptional\nsplicing"
                 ,"Export"
                 ,"Polysomal association"
                 ,"Cytoplasmatic\ndegradation"
                 ,"Nucleoplasmic mature\ndegradation"
                 ,"Polysomal degradation"
                 ,"Nucleoplasmic premature\ndegradation")
  
  names(RateNames) <- paste0("k",1:11)
  
  RateNames <- RateNames[c("k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11")%in%colnames(inferedRates)]
  myPalette <- myPalette[c("k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11")%in%colnames(inferedRates)]
  plot_list <- list()
  
  names(AllRates) <- AllRates
  AllRates <- AllRates[ratesOrder[ratesOrder%in%AllRates]]
  
  inferedRates <- inferedRates[,AllRates]
  
  xLabs <- c("Log10(fg/(M.cells*h))",rep("Log10(1/h)",10))
  names(xLabs) <- paste0("k",1:11)
  xLabs <- xLabs[AllRates]
  
  for (i in 1:length(AllRates)){
    if(length(grep(AllRates[i], MissingRates))==0){
      inferedRates_curr <- log10(data.frame(rate = inferedRates[, AllRates[i]]))
      rownames(inferedRates_curr) <- c()
      
      plot_list[[i]] <- ggplot(inferedRates_curr, aes(x = rate)) +
        geom_histogram(aes(y = ..density..), fill=myPalette[i], bins=100) +
        geom_density(alpha=.4, fill="grey") +
        theme_classic() +
        labs(x=xLabs[[i]], y="density") +
        ggtitle(paste0(RateNames[AllRates[[i]]])) +
        xlim(xlimTmp[[1]],xlimTmp[[2]]) +
        theme(plot.title = element_text(size = 15,hjust = 0.5))
    }
  }
  
  w <- sapply(seq_along(plot_list),function(i){is.null(plot_list[[i]])})
  img <- do.call("grid.arrange", c(plot_list[!w],ncol = length(ActualRates))) 
  ggsave(filename=paste0("RatesDistributions",obj_name,".pdf"),plot = img,width=width,height=height)
  
}

DifferentialRatesDistributions <- function(objects
                                           ,width
                                           ,height
                                           ,obj_name
                                           ,allRates=FALSE
                                           ,ratesOrder=c("k1","k2","k3","k4","k5","k11","k9","k6","k8","k7","k10"))
{
  namesTmp <- names(objects)
  inferedRates <- lapply(objects,function(object)object$inferedRates)
  
  ActualRates <- names(which(table(unlist(sapply(inferedRates,colnames)))==length(objects)))
  
  RateNames <- c("Synthesis\n[fg/(M. cells*h)]"
                 ,"Co-trascriptional\nsplicing\n[h^-1]"
                 ,"Detachment of\nmature chromatin\n[h^-1]"
                 ,"Detachment of\npremature chromatin\n[h^-1]"
                 ,"Post-transcriptional\nsplicing\n[h^-1]"
                 ,"Export\n[h^-1]"
                 ,"Polysomal association\n[h^-1]"
                 ,"Cytoplasmatic\ndegradation\n[h^-1]"
                 ,"Nucleoplasmic mature\ndegradation\n[h^-1]"
                 ,"Polysomal degradation\n[h^-1]"
                 ,"Nucleoplasmic premature\ndegradation\n[h^-1]")
  names(RateNames) <- paste0("k",1:11)
  
  if(!allRates)
  {
    inferedRates <- lapply(inferedRates,function(i)i[,ActualRates])
    RateNames <- RateNames[names(RateNames)%in%ActualRates]
  }
  myPalette <- RColorBrewer::brewer.pal(10,"Set3")
  names(myPalette) <- c("k1","k2","k3","k4","k5","k6","k7","k8","k9","k10")
  myPalette <- myPalette[names(RateNames)]
  
  sortedRates <- names(RateNames)
  names(sortedRates) <- sortedRates
  
  sortedRates <- sortedRates[ratesOrder[ratesOrder%in%sortedRates]]
  
  pdf(paste0("DifferentialRatesDistributions_",obj_name,".pdf"),width=width,height=height)
  {
    par(mfrow=c(1,length(RateNames)))
    for(i in sortedRates)
    {
      plotTmp <- lapply(inferedRates,function(j)tryCatch(j[,i],error=function(e){NA}))
      if(any(unlist(sapply(plotTmp,is.finite))))
      {
        boxplot(plotTmp,ylab="",main=RateNames[[i]],varwidth=TRUE,outline=FALSE,las=2)        
      }else{print(paste0("Rate ",i," is never estimated."))}
    }      
  }
  dev.off()
}

differentialDistributionPlot <- function(object1=object1
                                         ,object2=object2
                                         ,name=""
                                         ,width=6
                                         ,height=4
                                         ,breaks=25
                                         ,col1="palevioletred"
                                         ,col2="lightblue"
                                         ,name1=""
                                         ,name2=""
                                         ,pchMed=NA)
{
  ratesTmp1 <- object1$inferedRates
  ratesTmp2 <- object2$inferedRates
  
  if(ncol(ratesTmp1)==9){
    ratesTmp1 <- ratesTmp1[,c(1:6,8,7,9)]
    colnames(ratesTmp1) <- paste0("k",1:9)
    
    ratesTmp2 <- ratesTmp2[,c(1:6,8,7,9)]
    colnames(ratesTmp2) <- paste0("k",1:9)
  }else{
    colnames(ratesTmp1) <- paste0("k",1:7)
    colnames(ratesTmp2) <- paste0("k",1:7)
  }
  
  df1 <- as.data.frame(log10(c(ratesTmp1)))
  names(df1) <- "values"
  df1$rate <- as.factor(c(sapply(colnames(ratesTmp1),function(i)rep(i,nrow(ratesTmp1)))))
  
  df2 <- as.data.frame(log10(c(ratesTmp2)))
  names(df2) <- "values"
  df2$rate <- as.factor(c(sapply(colnames(ratesTmp2),function(i)rep(i,nrow(ratesTmp2)))))
  
  pdf(paste0("differentialDistributions",name,".pdf"),width=width,height=height)
  histoplot(values~rate,breaks=breaks,data=df1, col = col1,border = NA, side = "left",colMed="black",pchMed = pchMed, plotCentre = "line")
  histoplot(values~rate,breaks=breaks,data=df2, col = col2,border = NA, side = "right", add = T,colMed="black",pchMed = pchMed, plotCentre = "line")
  title(xlab = "", ylab = "")
  legend("top", fill = c(col1, col2), legend = c(name1,name2), title = "",box.lwd=-1,ncol=2)
  dev.off()
}

HeatmapsPlot <- function(object
                         ,obj_name
                         ,n_clust
                         ,ratesWeight=1
                         ,expressionFlag
                         ,ratesOrder=c("k1","k2","k3","k4","k5","k11","k9","k6","k8","k7","k10"))
{
  if(is.numeric(expressionFlag))
  {
    if(length(object$expressionData)<expressionFlag){print("Error: incorrect option selected")}
    ExampleData <- object$expressionData[[expressionFlag]]
    inferedRates <- object$inferedRates
  }else if(expressionFlag=="MEAN")
  {
    if(length(object$expressionData)<2){print("Error: incorrect option selected")}
    ExampleData <- lapply(object$expressionData,function(i){
      i[i==1e-10] <- NA
      i
    })
    ExampleData <- sapply(colnames(ExampleData[[1]]),function(i)apply(cbind(ExampleData[[1]][,i],ExampleData[[2]][,i]),1,mean,na.rm=TRUE))
    inferedRates <- object$inferedRates
  }else if(expressionFlag=="ALL")
  {
    if(length(object$expressionData)<2){print("Error: incorrect option selected")}
    ExampleData <- Reduce(rbind,object$expressionData)
    rownames(ExampleData) <- paste0(unlist(lapply(object$expressionData,rownames)),"_",sapply(seq_along(object$expressionData),function(i)rep(i,nrow(object$expressionData[[i]]))))
    
    inferedRates <- object$inferedRates[sapply(strsplit(rownames(ExampleData),"_"),"[[",1),]
    rownames(inferedRates) <- rownames(ExampleData)
  }
  
  polyFlag <- ("pyt"%in%colnames(ExampleData))
  
  if(polyFlag)
  {
    ExampleData <- ExampleData[,c("chpt","chmt","npt","nmt","cyt","pyt")]
  }else{
    ExampleData <- ExampleData[,c("chpt","chmt","npt","nmt","cyt")]
  }
  
  colnames(ExampleData) <- c("Chp","Chm","Np","Nm","C","P")[1:ncol(ExampleData)]
  
  sortedRates <- colnames(inferedRates)
  names(sortedRates) <- sortedRates
  
  sortedRates <- sortedRates[ratesOrder[ratesOrder%in%sortedRates]]
  
  inferedRates <- inferedRates[,sortedRates]
  
  # RateNames <- c("Syn","CoSpl","ChmDet","ChpDet","PostSpl","Exp","Tran","UnTranDeg","NucDeg","TranDeg")
  RateNames <- c("k1","k2","k3","k4","k5","k6","k8","k7","k10","k9")
  names(RateNames) <- paste0("k",1:10)
  
  colnames(inferedRates) <- RateNames[colnames(inferedRates)]
  
  pheatmap(cor(inferedRates,method="s")
           ,cluster_rows = FALSE,cluster_cols = FALSE
           ,display_numbers = TRUE,color = colorRampPalette(c("blue","white","red"))(101),breaks=seq(-1,1,length=100),
           ,filename = paste0("correlationsHeatmap_",obj_name,".pdf"),width=4,height=4)
  
  ExampleData <- log10(ExampleData)
  inferedRates <- log10(inferedRates)
  
  ExampleDataBinned <- apply(ExampleData,2,function(i)
  {
    i[i>quantile(i,0.99)] <- quantile(i,0.99)
    i[i<quantile(i,0.01)] <- quantile(i,0.01)
    
    if(all(i>0))
    {
      as.numeric(cut(i, seq(min(i)*0.99,max(i)*1.01,length.out=101)))  
    }else if(all(i<0)){
      as.numeric(cut(i, seq(min(i)*1.01,max(i)*0.99,length.out=101)))  
    }else{
      as.numeric(cut(i, seq(min(i)*1.01,max(i)*1.01,length.out=101)))  
    }
    
  })
  rownames(ExampleDataBinned) <- rownames(ExampleData)
  colnames(ExampleDataBinned) <- colnames(ExampleData)
  
  inferedRatesBinned <- apply(inferedRates,2,function(i)
  {
    i[i>quantile(i,0.99)] <- quantile(i,0.99)
    i[i<quantile(i,0.01)] <- quantile(i,0.01)
    
    if(all(i>0))
    {
      as.numeric(cut(i, seq(min(i)*0.99,max(i)*1.01,length.out=101)))  
    }else if(all(i<0)){
      as.numeric(cut(i, seq(min(i)*1.01,max(i)*0.99,length.out=101)))  
    }else{
      as.numeric(cut(i, seq(min(i)*1.01,max(i)*1.01,length.out=101)))  
    }
  })
  rownames(inferedRatesBinned) <- rownames(ExampleData)
  
  foe <- pheatmap(cbind(ExampleDataBinned,do.call("cbind",lapply(1:ratesWeight,function(i)inferedRatesBinned)))
                  ,cluster_rows = TRUE,cluster_cols = FALSE,
                  ,show_rownames = FALSE,color = colorRampPalette(c("blue","white","red"))(101)
                  ,silent=TRUE)
  
  ord <- foe$tree_row$order
  ord2 <- cutree(foe$tree_row,n_clust)[ord]
  
  ord2 <- split(ord2,ord2)[unique(ord2)]
  
  ord2 <- unlist(sapply(seq_along(ord2),function(i){rep(LETTERS[i],length(ord2[[i]]))}))
  names(ord2) <- rownames(ExampleData)[ord]
  
  rowAnnotation <- as.data.frame(list("Cluster"=ord2))
  
  rowAnnotationColors <- carto_pal(length(unique(rowAnnotation[,"Cluster"])), "Safe")
  names(rowAnnotationColors) <- unique(rowAnnotation[,"Cluster"])
  
  pheatmap(cbind(ExampleDataBinned,inferedRatesBinned)[ord,]
           ,cluster_rows = FALSE,cluster_cols = FALSE,annotation_row=rowAnnotation
           ,show_rownames = FALSE,color = colorRampPalette(c("white","pink","red"))(101)
           ,filename = paste0("expressionAndRatesHeatmap_",obj_name,".pdf"),width=5,height=4
           ,annotation_colors=list("Cluster"=rowAnnotationColors))
  
  saveRDS(ord2,file=paste0("Clustering_",obj_name,".rds"))
  return(ord2)
}

DifferentialHeatmapsPlot <- function(refObject
                                     ,object
                                     ,name
                                     ,n_clust
                                     ,expressionFlag
                                     ,ylim=c(-5,5)
                                     ,ratesWeight=1
                                     ,ratesOrder=c("k1","k2","k3","k4","k5","k11","k9","k6","k8","k7","k10")
                                     ,show_rownames=FALSE
                                     ,width,height
                                     ,externalAnnotationRow=NULL
                                     ,externalAnnotationCol=NULL
                                     ,externalAnnotationColors=NULL)
{
  if(is.numeric(expressionFlag))
  {
    if(length(object$expressionData)<expressionFlag){print("Error: incorrect option selected")}
    ExampleData1 <- refObject$expressionData[[expressionFlag]]
    ExampleData2 <- object$expressionData[[expressionFlag]]
    
    inferedRates1 <- refObject$inferedRates
    inferedRates2 <- object$inferedRates 
    
  }else if(expressionFlag=="MEAN")
  {
    if(length(object$expressionData)<2){print("Error: incorrect option selected")}
    ExampleData1 <- lapply(refObject$expressionData,function(i){
      i[i==1e-10] <- NA
      i
    })
    ExampleData1 <- sapply(colnames(ExampleData1[[1]]),function(i)apply(cbind(ExampleData1[[1]][,i],ExampleData1[[2]][,i]),1,mean,na.rm=TRUE))
    ExampleData2 <- lapply(object$expressionData,function(i){
      i[i==1e-10] <- NA
      i
    })
    ExampleData2 <- sapply(colnames(ExampleData2[[1]]),function(i)apply(cbind(ExampleData2[[1]][,i],ExampleData2[[2]][,i]),1,mean,na.rm=TRUE))
    
    inferedRates1 <- refObject$inferedRates
    inferedRates2 <- object$inferedRates 
    
  }else if(expressionFlag=="ALL")
  {
    if(length(object$expressionData)<2){print("Error: incorrect option selected")}
    ExampleData1 <- Reduce(rbind,refObject$expressionData)
    ExampleData2 <- Reduce(rbind,object$expressionData)
    
    rownames(ExampleData1) <- paste0(unlist(lapply(refObject$expressionData,rownames)),"_",sapply(seq_along(refObject$expressionData),function(i)rep(i,nrow(refObject$expressionData[[i]]))))
    rownames(ExampleData2) <- paste0(unlist(lapply(object$expressionData,rownames)),"_",sapply(seq_along(object$expressionData),function(i)rep(i,nrow(object$expressionData[[i]]))))
    
    inferedRates1 <- refObject$inferedRates[sapply(strsplit(rownames(ExampleData1),"_"),"[[",1),]
    rownames(inferedRates1) <- rownames(ExampleData1)
    
    inferedRates2 <- object$inferedRates[sapply(strsplit(rownames(ExampleData2),"_"),"[[",1),]
    rownames(inferedRates2) <- rownames(ExampleData2)
    
  }
  
  polyFlag1 <- ("pyt"%in%colnames(ExampleData1))
  if(polyFlag1)
  {
    ExampleData1 <- ExampleData1[,c("chpt","chmt","npt","nmt","cyt","pyt")]
  }else{
    ExampleData1 <- ExampleData1[,c("chpt","chmt","npt","nmt","cyt")]
  }
  
  polyFlag2 <- ("pyt"%in%colnames(ExampleData2))
  if(polyFlag2)
  {
    ExampleData2 <- ExampleData2[,c("chpt","chmt","npt","nmt","cyt","pyt")]
  }else{
    ExampleData2 <- ExampleData2[,c("chpt","chmt","npt","nmt","cyt")]
  }
  
  commonGenes <- intersect(rownames(ExampleData1),rownames(ExampleData2))
  commonRates <- intersect(colnames(inferedRates1),colnames(inferedRates2))
  commonSpecies <- intersect(colnames(ExampleData1),colnames(ExampleData2))
  
  names(commonRates) <- commonRates
  
  commonRates <- commonRates[ratesOrder[ratesOrder%in%commonRates]]
  
  print(paste0(length(commonGenes)," common genes."))
  
  cor1 <- cor(log10(inferedRates1[,commonRates]),method="s")
  cor2 <- cor(log10(inferedRates2[,commonRates]),method="s")
  
  corA <- corB <- cor2-cor1
  for(i in 1:nrow(corA))
  {
    for(j in 1:ncol(corA))
    {
      if(j<i){corA[i,j]=sign(cor1[i,j])*sign(cor2[i,j])}
    }
  }
  
  for(i in 1:nrow(corB))
  {
    for(j in 1:ncol(corB))
    {
      if(j<i){corB[i,j]=cor2[i,j]}
    }
  }
  
  # pheatmap(corA
  # ,cluster_rows = FALSE,cluster_cols = FALSE
  # ,display_numbers = TRUE,color = colorRampPalette(c("blue","white","red"))(101),breaks=seq(-0.5,0.5,length=100),
  # ,filename = paste0("differentialCorrelationsHeatmap_",name,".pdf"),width=4,height=4)
  
  RateNames <- c("k1","k2","k3","k4","k5","k6","k8","k7","k10","k9")
  names(RateNames) <- paste0("k",1:10)
  
  colnames(corB) <- rownames(corB) <- RateNames[colnames(corB)]
  
  pheatmap(corB
           ,cluster_rows = FALSE,cluster_cols = FALSE
           ,display_numbers = TRUE,color = colorRampPalette(c("blue","white","red"))(101),breaks=seq(-1,1,length=100),
           ,filename = paste0("differentialCorrelationsHeatmap_",name,".pdf"),width=4,height=4)
  
  ExampleData1 <- ExampleData1[commonGenes,commonSpecies]
  inferedRates1 <- inferedRates1[commonGenes,commonRates]
  ExampleData2 <- ExampleData2[commonGenes,commonSpecies]
  inferedRates2 <- inferedRates2[commonGenes,commonRates]
  
  ExampleData <- log2(ExampleData2/ExampleData1)  
  inferedRates <- log2(inferedRates2/inferedRates1)  
  
  colnames(inferedRates) <- RateNames[colnames(inferedRates)]
  
  ExampleData[ExampleData<ylim[[1]]] <- ylim[[1]]
  inferedRates[inferedRates<ylim[[1]]] <- ylim[[1]]
  ExampleData[ExampleData>ylim[[2]]] <- ylim[[2]]
  inferedRates[inferedRates>ylim[[2]]] <- ylim[[2]]
  
  if(!is.null(externalAnnotationRow))
  {
    ExampleData <- ExampleData[rownames(ExampleData)%in%rownames(externalAnnotationRow),]
    externalAnnotationRow <- externalAnnotationRow[rownames(ExampleData),]
  }
  
  foe <- pheatmap(cbind(ExampleData,do.call("cbind",lapply(1:ratesWeight,function(i)inferedRates)))
                  ,cluster_rows = TRUE,cluster_cols = FALSE,
                  ,show_rownames = FALSE,color = colorRampPalette(c("blue","white","red"))(101)
                  ,silent=TRUE)
  
  ord <- foe$tree_row$order
  ord2 <- cutree(foe$tree_row,n_clust)[ord]
  
  ord2 <- split(ord2,ord2)[unique(ord2)]
  
  ord2 <- unlist(lapply(seq_along(ord2),function(i){rep(LETTERS[i],length(ord2[[i]]))}))
  names(ord2) <- rownames(ExampleData)[ord]
  
  rowAnnotation <- as.data.frame(list("Cluster"=ord2))
  
  rowAnnotationColors <- carto_pal(length(unique(rowAnnotation[,"Cluster"])), "Safe")
  names(rowAnnotationColors) <- unique(rowAnnotation[,"Cluster"])
  
  if(!is.null(externalAnnotationRow))
  {
    rowAnnotation <- cbind(rowAnnotation,externalAnnotationRow[rownames(rowAnnotation),])
  }
  
  if(!is.null(externalAnnotationCol))
  {
    externalAnnotationColTmp <- matrix(NA,nrow=ncol(cbind(ExampleData,inferedRates)),ncol=1)
    rownames(externalAnnotationColTmp) <- colnames(cbind(ExampleData,inferedRates))
    colnames(externalAnnotationColTmp) <- "DifferentialGenes"
    externalAnnotationColTmp[names(externalAnnotationCol),1] <- externalAnnotationCol
    externalAnnotationCol <- as.data.frame(externalAnnotationColTmp)
  }
  
  pheatmap(cbind(ExampleData,inferedRates)[ord,]
           ,cluster_rows = FALSE,cluster_cols = FALSE,annotation_row=rowAnnotation,annotation_col=externalAnnotationCol
           ,annotation_colors=append(list("Cluster"=rowAnnotationColors),externalAnnotationColors)
           ,show_rownames = show_rownames,color = colorRampPalette(c("blue","white","red"))(101)
           ,filename = paste0("differentialExpressionAndRatesHeatmap_",name,".pdf"),width=width,height=height)
  
  saveRDS(ord2,file=paste0("DifferentialClustering_",name,".rds"))
  return(ord2)
}

GenesCharacterization <- function(txdb,clustering=NULL,sortedRates=NULL,top=NULL,bottom=NULL,obj_name,cgFeatures=NULL,entropyFeatures=NULL,width=8,height=6,ids=LETTERS){
  
  require("TxDb.Hsapiens.UCSC.hg38.knownGene")
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  foe <- transcriptsBy(txdb,'gene')
  #extract n° of isoforms per transcripts
  n_isoforms <- unlist(as.list(table(names(unlist(foe)))))
  
  #extract number of exons and their length per transcripts
  exonsDB <- reduce(exonsBy(txdb ,'gene'))
  exonsDB <- exonsDB[elementNROWS(range(exonsDB))==1]
  n_exons <- unlist(as.list(table(names(unlist(exonsDB)))))
  exonsLength <- sapply(width(exonsDB),sum)
  
  #extract number of introns and their length per transcripts
  intronsDB <- psetdiff(unlist(range(exonsDB)),exonsDB)
  intronsDB <- intronsDB[elementNROWS(intronsDB)>0]
  intronsLength <- sapply(width(intronsDB), sum)
  
  tx2gene <- unlist(transcriptsBy(txdb,"gene"))
  tx2gene <- names(tx2gene)
  names(tx2gene) <- unlist(transcriptsBy(txdb,"gene"))$tx_name
  
  #Extract utr3 median lengths
  utr3 <-threeUTRsByTranscript(txdb,use.names=TRUE)
  utr3 <- unlist(width(utr3))
  names(utr3) <- tx2gene[names(utr3)]
  utr3max <- sapply(split(utr3,names(utr3)),max)
  utr3 <- sapply(split(utr3,names(utr3)),median)
  
  #Extract utr5 median lengths
  utr5 <-fiveUTRsByTranscript(txdb,use.names=TRUE)
  utr5 <- unlist(width(utr5))
  names(utr5) <- tx2gene[names(utr5)]
  utr5max <- sapply(split(utr5,names(utr5)),max)
  utr5 <- sapply(split(utr5,names(utr5)),median)
  
  if(is.matrix(sortedRates))
  {
    
    #Characterization of genes with fastest and slowest dynamics
    
    
    ratesNames <- c(rep(unlist(lapply(colnames(sortedRates),function(i){rep(i,top)})),2),
                    unlist(lapply(colnames(sortedRates),function(i){rep(i,nrow(sortedRates))})))
    df <- data.frame(cbind("Rank"=c(rep("TOP",top*ncol(sortedRates)),rep("BOT",bottom*ncol(sortedRates)),rep("ALL",nrow(sortedRates)*9)),
                           "Rate"=ratesNames,
                           "Genes"=c(matrix(sortedRates[1:top,],ncol=1),
                                     matrix(sortedRates[(nrow(sortedRates)-bottom+1):nrow(sortedRates),],ncol=1),
                                     matrix(sortedRates,ncol=1))))
    
    pdf(file=paste0("StructuralFeatures_TopBot",obj_name,".pdf"),width=width,height=height)
    if(!is.null(cgFeatures)|!is.null(entropyFeatures)){par(mfrow=c(3,4))}else{par(mfrow=c(2,4))}
    
    df$utr5 <- utr5[df$Genes]/1000
    df$utr5max <- utr5max[df$Genes]/1000
    df$utr3 <- utr3[df$Genes]/1000
    df$utr3max <- utr3max[df$Genes]/1000
    df$exonsLength <- exonsLength[df$Genes]/1000
    df$intronsLength <- intronsLength[df$Genes]/1000
    df$n_isoforms <- n_isoforms[df$Genes]
    df$n_exons <- n_exons[df$Genes]
    df$cg_3UTR <- cgFeatures[[1]][df$Genes]
    df$cg_5UTR <- cgFeatures[[2]][df$Genes]
    df$entropy_3UTR <- entropyFeatures[[1]][df$Genes]
    df$entropy_5UTR <- entropyFeatures[[2]][df$Genes]
    
    df$Rate <- factor(df$Rate, levels=colnames(sortedRates))
    df$Rank <- factor(df$Rank, levels=c("ALL","TOP","BOT",""))
    
    Wtest <- lapply(colnames(df)[4:15],function(j){
      
      feat <- lapply(levels(df$Rate),function(i){
        
        w1 = wilcox.test(df[df$Rank=="ALL"&df$Rate==i,j],df[df$Rank=="TOP"&df$Rate==i,j])$p.value
        w2 = wilcox.test(df[df$Rank=="ALL"&df$Rate==i,j],df[df$Rank=="BOT"&df$Rate==i,j])$p.value
        
        return(c("top"=w1,"bot"=w2))
      })
      names(feat) <- levels(df$Rate)
      return(unlist(feat))
    })
    names(Wtest) <- colnames(df)[4:15]
    
    MatCol <- matrix("grey",nrow=4*ncol(sortedRates),ncol=12)
    rownames(MatCol) <- sapply(c(paste0("k",1:8),"k10"),function(i)paste0(i,".",c("all","top","bot","")))
    colnames(MatCol) <- names(Wtest)
    for(i in names(Wtest)){
      
      MatCol[names(Wtest[[i]]),i] <- c("blue","red")[(Wtest[[i]]<1e-4)+1]
    }
    
    boxplot(utr5~Rank+Rate,data=df,outline=FALSE,las=2,ylab="[kbs]",main="Median Length\n5' UTR",varwidth=FALSE,col=MatCol[,"utr5"],xlab=NULL,cex.axis = 0.8)
    
    boxplot(utr5max~Rank+Rate,data=df,outline=FALSE,las=2,ylab="[kbs]",main="Max Length\n5' UTR",varwidth=FALSE,col=MatCol[,"utr5max"],xlab=NULL,cex.axis = 0.8)
    
    boxplot(utr3~Rank+Rate,data=df,outline=FALSE,las=2,ylab="[kbs]",main="Median Length\n3' UTR",varwidth=FALSE,col=MatCol[,"utr3"],xlab=NULL,cex.axis = 0.8)
    
    boxplot(utr3max~Rank+Rate,data=df,outline=FALSE,las=2,ylab="[kbs]",main="Max Length\n3' UTR",varwidth=FALSE,col=MatCol[,"utr3max"],xlab=NULL,cex.axis = 0.8)
    
    boxplot(exonsLength~Rank+Rate,data=df,outline=FALSE,las=2,ylab="[kbs]",main="Exons length",varwidth=FALSE,col=MatCol[,"exonsLength"],xlab=NULL,cex.axis = 0.8)
    
    boxplot(intronsLength~Rank+Rate,data=df,outline=FALSE,las=2,ylab="[kbs]",main="Introns length",varwidth=FALSE,col=MatCol[,"intronsLength"],xlab=NULL,cex.axis = 0.8)
    
    boxplot(n_isoforms~Rank+Rate,data=df,outline=FALSE,las=2,ylab="[#]",main="Isoforms number",varwidth=FALSE,col=MatCol[,"n_isoforms"],xlab=NULL,cex.axis = 0.8)
    
    boxplot(n_exons~Rank+Rate,data=df,outline=FALSE,las=2,ylab="[#]",main="Exons number",varwidth=FALSE,col=MatCol[,"n_exons"],xlab=NULL,cex.axis = 0.8)
    
    if(!is.null(cgFeatures))
    {
      boxplot(cg_3UTR~Rank+Rate,data=df,outline=FALSE,las=2,ylab="[%]",main="CG percentage\nMax 5' UTR",varwidth=FALSE,col=MatCol[,"cg_3UTR"],xlab=NULL,cex.axis = 0.8)
      boxplot(cg_5UTR~Rank+Rate,data=df,outline=FALSE,las=2,ylab="[%]",main="CG percentage\nMax 3' UTR",varwidth=FALSE,col=MatCol[,"cg_5UTR"],xlab=NULL,cex.axis = 0.8)
    }
    
    if(!is.null(entropyFeatures))
    {
      boxplot(entropy_3UTR~Rank+Rate,data=df,outline=FALSE,las=2,ylab="[bits]",main="Entropy\nMax 5' UTR",varwidth=FALSE,col=MatCol[,"entropy_3UTR"],xlab=NULL,cex.axis = 0.8)
      boxplot(entropy_5UTR~Rank+Rate,data=df,outline=FALSE,las=2,ylab="[bits]",main="Entropy\nMax 3' UTR",varwidth=FALSE,col=MatCol[,"entropy_5UTR"],xlab=NULL,cex.axis = 0.8)
    }
    
    dev.off()
    
  }else{
    
    #Gene characterization based on clustering
    clustering <- split(names(clustering),clustering)
    
    pdf(file=paste0("clustersStructuralFeatures_",obj_name,".pdf"),width=width,height=height)
    if(!is.null(cgFeatures)|!is.null(entropyFeatures)){par(mfrow=c(3,4))}else{par(mfrow=c(2,4))}
    
    boxplot(lapply(seq_along(clustering),function(i){utr5[names(utr5)%in%clustering[[i]]]/1000}),ylab="[kbs]",main="Median Length\n5' UTR",names=ids[seq_along(clustering)],outline=FALSE,varwidth=TRUE,las=2)
    boxplot(lapply(seq_along(clustering),function(i){utr5max[names(utr5max)%in%clustering[[i]]]/1000}),ylab="[kbs]",main="Max Length\n5' UTR",names=ids[seq_along(clustering)],outline=FALSE,varwidth=TRUE,las=2)
    boxplot(lapply(seq_along(clustering),function(i){utr3[names(utr3)%in%clustering[[i]]]/1000}),ylab="[kbs]",main="Median Length\n3' UTR",names=ids[seq_along(clustering)],outline=FALSE,varwidth=TRUE,las=2)
    boxplot(lapply(seq_along(clustering),function(i){utr3max[names(utr3max)%in%clustering[[i]]]/1000}),ylab="[kbs]",main="Max Length\n3' UTR",names=ids[seq_along(clustering)],outline=FALSE,varwidth=TRUE,las=2)
    boxplot(lapply(seq_along(clustering),function(i){exonsLength[names(exonsLength)%in%clustering[[i]]]/1000}),ylab="[kbs]",main="Exons length",names=ids[seq_along(clustering)],outline=FALSE,varwidth=TRUE,las=2)
    boxplot(lapply(seq_along(clustering),function(i){intronsLength[names(intronsLength)%in%clustering[[i]]]/1000}),ylab="[kbs]",main="Introns length",names=ids[seq_along(clustering)],outline=FALSE,varwidth=TRUE,las=2)
    boxplot(lapply(seq_along(clustering),function(i){n_isoforms[names(n_isoforms)%in%clustering[[i]]]}),ylab="[#]",main="Isoforms number",names=ids[seq_along(clustering)],outline=FALSE,varwidth=TRUE,las=2)
    boxplot(lapply(seq_along(clustering),function(i){n_exons[names(n_exons)%in%clustering[[i]]]}),ylab="[#]",main="Exons number",names=ids[seq_along(clustering)],outline=FALSE,varwidth=TRUE,las=2)
    
    if(!is.null(cgFeatures))
    {
      boxplot(lapply(seq_along(clustering),function(i){cgFeatures[[1]][names(cgFeatures[[1]])%in%clustering[[i]]]}),ylab="[%]",main="CG percentage\nMax 5' UTR",names=ids[seq_along(clustering)],outline=FALSE,varwidth=TRUE,las=2)
      boxplot(lapply(seq_along(clustering),function(i){cgFeatures[[2]][names(cgFeatures[[2]])%in%clustering[[i]]]}),ylab="[%]",main="CG percentage\nMax 3' UTR",names=ids[seq_along(clustering)],outline=FALSE,varwidth=TRUE,las=2)
    }
    
    if(!is.null(entropyFeatures))
    {
      boxplot(lapply(seq_along(clustering),function(i){entropyFeatures[[1]][names(entropyFeatures[[1]])%in%clustering[[i]]]}),ylab="[bits]",main="Entropy\nMax 5' UTR",names=ids[seq_along(clustering)],outline=FALSE,varwidth=TRUE,las=2)
      boxplot(lapply(seq_along(clustering),function(i){entropyFeatures[[2]][names(entropyFeatures[[2]])%in%clustering[[i]]]}),ylab="[bits]",main="Entropy\nMax 3' UTR",names=ids[seq_along(clustering)],outline=FALSE,varwidth=TRUE,las=2)
    }
    
    dev.off()
  }
}


TailsCharacterization <- function(txdb,clustering=NULL,obj_name,bam_untreated=NULL,bam_treated,tailsPaths_untreatetd=NULL,
                                  tailsPaths_treated,minoverlap_I,minoverlap_E,mergeReplicates,width,height)
{
  if(mergeReplicates)
  {
    if(!is.null(bam_untreated))names(bam_untreated) <- gsub("1$","",gsub("2$","",names(bam_untreated)))
    names(bam_treated) <- gsub("1$","",gsub("2$","",names(bam_treated)))
    if(!is.null(tailsPaths_untreatetd))names(tailsPaths_untreatetd) <- gsub("1$","",gsub("2$","",names(tailsPaths_untreatetd)))
    names(tailsPaths_treated) <- gsub("1$","",gsub("2$","",names(tailsPaths_treated)))
  }
  
  bams <- list(list("BAM"=bam_untreated,"TAIL"=tailsPaths_untreatetd),list("BAM"=bam_treated,"TAIL"=tailsPaths_treated))
  
  Out <- lapply(bams,function(j)
  {
    if(!is.null(j[[1]]))
    {
      outTmp <- lapply(sort(unique(names(j$BAM))),function(i)
      {
        print(i)
        
        bamTmp <- j$BAM[names(j$BAM)==i]
        tailTmp <- j$TAIL[names(j$TAIL)==i]
        
        # BAM loading
        if(length(bamTmp>1))
        {
          bamTmp2 <- get(load(bamTmp[[1]]))
          for(i in bamTmp[-1])
          {
            bamTmp2 <- c(bamTmp2,get(load(i)))
          }
          bamTmp <- bamTmp2
        }else{bamTmp <- get(load(bamTmp))}
        
        # Definition of regions classified as exonic in at-least one isoform
        exonsDB <- reduce(exonsBy(txdb ,'gene'))
        exonsDB <- exonsDB[elementNROWS(range(exonsDB))==1]
        
        # Definition of gaps as intronic regions
        intronsDB <- psetdiff(unlist(range(exonsDB)),exonsDB)
        intronsDB <- intronsDB[elementNROWS(intronsDB)>0]
        
        # Correct chromosomes names
        seqlevelsStyle(exonsDB) <- "ENSEMBL"
        seqlevelsStyle(intronsDB) <- "ENSEMBL"
        
        # Overlap between reads and exonic regions
        geneOverlaps <- findOverlaps(grglist(bamTmp,drop.D.ranges=TRUE),unlist(range(exonsDB)),minoverlap=minoverlap_E)
        
        # Unique overlaps
        geneOverlaps <- geneOverlaps[isUnique(queryHits(geneOverlaps)),]
        bamTmp <- bamTmp[names(bamTmp)[queryHits(geneOverlaps)]]
        
        # Overlap between reads and intronic regions
        intronicOverlaps <- findOverlaps(grglist(bamTmp,drop.D.ranges=TRUE),intronsDB,minoverlap=minoverlap_I)
        intronicOverlaps <- intronicOverlaps[isUnique(queryHits(intronicOverlaps)),]
        intronicReads <- names(bamTmp)[queryHits(intronicOverlaps)]
        names(intronicReads) <- names(intronsDB)[subjectHits(intronicOverlaps)]
        
        bamTmp <- bamTmp[!(names(bamTmp)%in%intronicReads)]
        
        # Annotation of no-intronic reads
        exonicOverlaps <- findOverlaps(grglist(bamTmp,drop.D.ranges=TRUE),exonsDB,minoverlap=minoverlap_E)
        exonicOverlaps <- exonicOverlaps[isUnique(queryHits(exonicOverlaps)),]
        exonicReads <- names(bamTmp)[queryHits(exonicOverlaps)]
        names(exonicReads) <- names(exonsDB)[subjectHits(exonicOverlaps)]
        
        # Tails loading
        tailTmp <- lapply(tailTmp,function(i)read.table(i,header=TRUE,sep="\t"))
        if(is.list(tailTmp)){tailTmp <- do.call("rbind",tailTmp)}
        tailTmp <- tailTmp[tailTmp[,"qc_tag"]=="PASS",]
        tailLengthTmp <- tailTmp[,"polya_length"]
        names(tailLengthTmp) <- tailTmp[,1]
        
        intronicTailLengths <- tailLengthTmp[names(tailLengthTmp)%in%intronicReads]
        exonicTailLengths <- tailLengthTmp[names(tailLengthTmp)%in%exonicReads]
        
        readToGene <- c(names(exonicReads),names(intronicReads))
        names(readToGene) <- c(exonicReads,intronicReads)
        
        allTailLengths <- split(c(intronicTailLengths,exonicTailLengths),readToGene[c(names(intronicTailLengths),names(exonicTailLengths))])
        
        intronicTailLengths <- split(intronicTailLengths,readToGene[names(intronicTailLengths)])
        exonicTailLengths <- split(exonicTailLengths,readToGene[names(exonicTailLengths)])
        
        list(intronic=intronicTailLengths,exonic=exonicTailLengths,all=allTailLengths)
      })
      names(outTmp) <- sort(unique(names(j$BAM)))
      return(outTmp)
    }else{list(list(intronic=NULL,exonic=NULL,all=NULL)
               ,list(intronic=NULL,exonic=NULL,all=NULL)
               ,list(intronic=NULL,exonic=NULL,all=NULL)
               ,list(intronic=NULL,exonic=NULL,all=NULL))}
  })
  names(Out) <- c("Untreated","Treated")
  
  tailsTmp <- lapply(Out,function(j){unlist(lapply(j,function(i)list(sapply(i[[1]],median),sapply(i[[2]],median),sapply(i[[3]],median))),recursive=FALSE)})
  
  tailsTmp <- lapply(tailsTmp,function(j)
  {
    names(j) <- sapply(names(j),function(i){gsub("1$","_In",i)})
    names(j) <- sapply(names(j),function(i){gsub("2$","_Ex",i)})
    names(j) <- sapply(names(j),function(i){gsub("3$","",i)})
    return(j)
  })
  
  tailsPlot <- lapply(tailsTmp,function(i){i[!(grepl("_In",names(i))|grepl("_Ex",names(i)))]})
  
  l <- unlist(unname(lapply(tailsPlot,function(i){sapply(i,length)})))
  s <- unlist(sapply(seq_along(l),function(i){rep(names(l[i]),l[i])}))
  
  cnd <- sapply(strsplit(s,"_"),'[',2)
  df <- as.data.frame(cbind("Species"=s,
                            "Condition"=cnd,
                            "Values"=unlist(tailsPlot)))
  df$Species <- sapply(strsplit(df$Species,"_"),'[',1)
  rownames(df) <-  NULL
  df$Values <- as.numeric(df$Values)
  df$GenesNames <- unname(unlist(lapply(tailsPlot,function(i){lapply(i,function(j){names(j)})})))
  
  df$Species <- factor(df$Species,levels=c("Chr","Nuc","Cyt","Poly"))
  
  pdf(paste0("polyATailLength_",obj_name,".pdf"),width=2.5,height=4)
  boxplot(Values~Condition+Species,data=df,outline=FALSE,ylab="Median PolyA tail\nlength [bs]",xlab=NULL,las=2)
  dev.off()
  
  if(!is.null(clustering))
  {
    df2 <- df[df$GenesNames%in%names(clustering),]
    df2$Clustering <- sapply(df2$GenesNames,function(i){clustering[i]})
    
    pdf(paste0("polyATailLengthClusters_",obj_name,".pdf"),width=width,height=height)
    par(mfrow=c(2,length(unique(clustering))/2))
    
    for(k in unique(clustering))
    {
      dfTmp <- df2[df2$Clustering==k,]
      dfTmp$Species <- factor(dfTmp$Species,levels=c("Chr","Nuc","Cyt","Poly"))
      boxplot(Values~Condition+Species,data=dfTmp,outline=FALSE,ylab="Median PolyA tail length [bs]",xlab=NULL,las=2,main=paste0("Cluster ",k))
      pl <- ggplot(dfTmp,aes(x=Species,y=Values,fill=Condition)) + 
        geom_boxplot(outlier.shape = NA) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.title.x=element_blank()) +
        ylab("Median PolyA tail\nlength [bs]")
    }
    dev.off()
  }
  return(Out)
}

cgContentFunction <- function(fastaPath)
{
  fastaTmp <- readLines(fastaPath)
  namesTmp <- fastaTmp[grepl(">",fastaTmp)]
  namesTmp <- gsub(">","",namesTmp)
  namesTmp <- sapply(strsplit(namesTmp,"[(]"),"[[",1)
  
  fastaTmp <- fastaTmp[!grepl(">",fastaTmp)]
  
  CGTmp <- sapply(fastaTmp,function(i)sum(strsplit(i,"")[[1]]=="C"|strsplit(i,"")[[1]]=="G")/length(strsplit(i,""
  )[[1]]))
  names(CGTmp) <- namesTmp
  CGTmp
}

entropyFunction <- function(fastaPath)
{
  fastaTmp <- readLines(fastaPath)
  namesTmp <- fastaTmp[grepl(">",fastaTmp)]
  namesTmp <- gsub(">","",namesTmp)
  namesTmp <- sapply(strsplit(namesTmp,"[(]"),"[[",1)
  
  fastaTmp <- fastaTmp[!grepl(">",fastaTmp)]
  
  CGTmp <- sapply(fastaTmp,function(i)Entropy(table(strsplit(i,"")[[1]])))
  names(CGTmp) <- namesTmp
  CGTmp
}

############### CORRELATION BETWEEN DATA ##################
## Goodness of fit
GoodnessOfFit <- function(inferedRates,expressionFlag,width=10,height=15,name="",lowSat=0,upSat=1)
{
  inferedDataTmp <- inferedRates$inferedData
  
  if(is.numeric(expressionFlag))
  {
    if(length(inferedRates$expressionData)<expressionFlag){print("Error: incorrect option selected")}
    expressionDataTmp <- inferedRates$expressionData[[expressionFlag]]
  }else if(expressionFlag=="MEAN")
  {
    if(length(inferedRates$expressionData)<2){print("Error: incorrect option selected")}
    expressionDataTmp <- lapply(inferedRates$expressionData,function(i){
      i[i==1e-10] <- NA
      i
    })
    expressionDataTmp <- sapply(colnames(expressionDataTmp[[1]]),function(i)apply(cbind(expressionDataTmp[[1]][,i],expressionDataTmp[[2]][,i]),1,mean,na.rm=TRUE))
  }else if(expressionFlag=="ALL")
  {
    if(length(inferedRates$expressionData)<2){print("Error: incorrect option selected")}
    expressionDataTmp <- Reduce(rbind,inferedRates$expressionData)
    inferedDataTmp <- inferedDataTmp[rownames(expressionDataTmp),]
    
    expressionDataTmp[expressionDataTmp==1e-10] <- NaN
  }
  
  colnames(expressionDataTmp) <- gsub("chpp","Chp_PreEx_",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("chpn","Chp_Nas_",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("chpt","Chp_Tot",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("chmp","Chm_PreEx_",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("chmn","Chm_Nas_",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("chmt","Chm_Tot",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("npp","Np_PreEx_",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("npn","Np_Nas_",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("npt","Np_Tot",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("nmp","Nm_PreEx_",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("nmn","Nm_Nas_",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("nmt","Nm_Tot",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("cyp","C_PreEx_",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("cyn","C_Nas_",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("cyt","C_Tot",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("pyp","P_PreEx_",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("pyn","P_Nas_",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("pyt","P_Tot",colnames(expressionDataTmp))
  colnames(expressionDataTmp) <- gsub("0.33","20min",colnames(expressionDataTmp))
  
  colnames(inferedDataTmp) <- gsub("chpp","Chp_PreEx_",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("chpn","Chp_Nas_",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("chpt","Chp_Tot",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("chmp","Chm_PreEx_",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("chmn","Chm_Nas_",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("chmt","Chm_Tot",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("npp","Np_PreEx_",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("npn","Np_Nas_",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("npt","Np_Tot",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("nmp","Nm_PreEx_",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("nmn","Nm_Nas_",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("nmt","Nm_Tot",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("cyp","C_PreEx_",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("cyn","C_Nas_",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("cyt","C_Tot",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("pyp","P_PreEx_",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("pyn","P_Nas_",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("pyt","P_Tot",colnames(inferedDataTmp))
  colnames(inferedDataTmp) <- gsub("0.33","20min",colnames(inferedDataTmp))
  
  name <- paste0("DataCorrelations_",name,".pdf")
  pdf(name,width=width,height=height)
  par(mfrow=c(6,3))
  
  outList <- list()
  
  for(i in colnames(inferedDataTmp)[!grepl("0$",colnames(inferedDataTmp))])
  {
    corTmp <- tryCatch(cor(inferedDataTmp[,i],expressionDataTmp[,i],method="s",use="c"),error=function(e)NaN)
    
    x <- log10(inferedDataTmp[,i])
    y <- log10(expressionDataTmp[,i])
    
    if(any(!is.finite(x))|any(!is.finite(y)))
    {
      print("Not finite data-points.")
      y <- y[is.finite(x)]
      x <- x[is.finite(x)]
      
      x <- x[is.finite(y)]
      y <- y[is.finite(y)]
    }
    
    xSat <- x
    ySat <- y
    
    xSatBounds <- quantile(xSat,c(lowSat,upSat))
    ySatBounds <- quantile(ySat,c(lowSat,upSat))
    
    xSat <- xSat[xSat<xSatBounds[[1]]|xSat>xSatBounds[[2]]|ySat<ySatBounds[[1]]|ySat>ySatBounds[[2]]]
    ySat <- ySat[names(xSat)]
    
    xSat[xSat<xSatBounds[[1]]] <- xSatBounds[[1]]
    xSat[xSat>xSatBounds[[2]]] <- xSatBounds[[2]]
    
    ySat[ySat<ySatBounds[[1]]] <- ySatBounds[[1]]
    ySat[ySat>ySatBounds[[2]]] <- ySatBounds[[2]]
    
    outList <- append(outList,corTmp)
    
    if(is.finite(corTmp))
    {
      par(mar = c(5,4,4,7) + .1)
      smoothScatter(x,y,postPlotHook = fudgeit
                    ,xlab="Inferred",ylab="Expected",main=i,xlim=xSatBounds,ylim=ySatBounds)
      abline(0,1,col=2,lwd=2)
      points(xSat,ySat,pch=16,col=2)
      tryCatch(points(lowess(x,y),col=3,lwd=2,type="l"),error=function(e){return(NULL)})
      legend("bottomright",legend=paste0("Cor.=",round(corTmp,2)),bty = "n",text.font=2)
    }
  }
  dev.off()
  
  names(outList) <- colnames(inferedDataTmp)[!grepl("0$",colnames(inferedDataTmp))]
  outList
}


replicatesDataCorrelations <- function(object1,object2,expressionFlag,width=7,height=12,name="",lowSat=0,upSat=1)
{
  
  if(is.numeric(expressionFlag))
  {
    if(length(object1$expressionData)<expressionFlag|length(object2$expressionData)<expressionFlag){print("Error: incorrect option selected")}
    expressionLevels1 <- object1$expressionData[[expressionFlag]]
    expressionLevels2 <- object2$expressionData[[expressionFlag]]
  }else if(expressionFlag=="MEAN")
  {
    if(length(object1$expressionData)<2|length(object2$expressionData)<2){print("Error: incorrect option selected")}
    
    expressionLevels1 <- lapply(object1$expressionData,function(i){
      i[i==1e-10] <- NA
      i
    })
    expressionLevels1 <- sapply(colnames(expressionLevels1[[1]]),function(i)apply(cbind(expressionLevels1[[1]][,i],expressionLevels1[[2]][,i]),1,mean,na.rm=TRUE))
    
    expressionLevels2 <- lapply(object1$expressionData,function(i){
      i[i==1e-10] <- NA
      i
    })
    expressionLevels2 <- sapply(colnames(expressionLevels2[[1]]),function(i)apply(cbind(expressionLevels2[[1]][,i],expressionLevels2[[2]][,i]),1,mean,na.rm=TRUE))
    
  }else if(expressionFlag=="ALL")
  {
    if(length(object1$expressionData)<2|length(object2$expressionData)<2){print("Error: incorrect option selected")}
    expressionLevels1 <- Reduce(rbind,object1$expressionData) 
    expressionLevels2 <- Reduce(rbind,object2$expressionData)
  }
  
  colnames(expressionLevels1) <- gsub("chpp","Chp_PreEx_",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("chpn","Chp_Nas_",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("chpt","Chp_Tot",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("chmp","Chm_PreEx_",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("chmn","Chm_Nas_",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("chmt","Chm_Tot",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("npp","Np_PreEx_",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("npn","Np_Nas_",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("npt","Np_Tot",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("nmp","Nm_PreEx_",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("nmn","Nm_Nas_",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("nmt","Nm_Tot",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("cyp","C_PreEx_",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("cyn","C_Nas_",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("cyt","C_Tot",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("pyp","P_PreEx_",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("pyn","P_Nas_",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("pyt","P_Tot",colnames(expressionLevels1))
  colnames(expressionLevels1) <- gsub("0.33","20min",colnames(expressionLevels1))
  
  colnames(expressionLevels2) <- gsub("chpp","Chp_PreEx_",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("chpn","Chp_Nas_",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("chpt","Chp_Tot",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("chmp","Chm_PreEx_",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("chmn","Chm_Nas_",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("chmt","Chm_Tot",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("npp","Np_PreEx_",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("npn","Np_Nas_",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("npt","Np_Tot",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("nmp","Nm_PreEx_",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("nmn","Nm_Nas_",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("nmt","Nm_Tot",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("cyp","C_PreEx_",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("cyn","C_Nas_",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("cyt","C_Tot",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("pyp","P_PreEx_",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("pyn","P_Nas_",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("pyt","P_Tot",colnames(expressionLevels2))
  colnames(expressionLevels2) <- gsub("0.33","20min",colnames(expressionLevels2))
  
  commonGenes <- intersect(rownames(object1$inferedRates),rownames(object2$inferedRates))
  
  outList <- list()
  
  pdf(paste0("rawDataCorrelations_",name,".pdf"),width=10,height=16)
  par(mfrow=c(6,3))
  for(i in colnames(expressionLevels1)[!grepl("0$",colnames(expressionLevels1))])
  {
    corTmp <- tryCatch(cor(expressionLevels1[commonGenes,i],expressionLevels2[commonGenes,i],method="s",use="c"),error=function(e)NaN)
    
    x <- log10(expressionLevels1[commonGenes,i])
    y <- log10(expressionLevels2[commonGenes,i])
    
    xSat <- x
    ySat <- y
    
    xSatBounds <- quantile(xSat,c(lowSat,upSat))
    ySatBounds <- quantile(ySat,c(lowSat,upSat))
    
    xSat <- xSat[xSat<xSatBounds[[1]]|xSat>xSatBounds[[2]]|ySat<ySatBounds[[1]]|ySat>ySatBounds[[2]]]
    ySat <- ySat[names(xSat)]
    
    xSat[xSat<xSatBounds[[1]]] <- xSatBounds[[1]]
    xSat[xSat>xSatBounds[[2]]] <- xSatBounds[[2]]
    
    ySat[ySat<ySatBounds[[1]]] <- ySatBounds[[1]]
    ySat[ySat>ySatBounds[[2]]] <- ySatBounds[[2]]
    
    commonGenes <- intersect(rownames(expressionLevels1),rownames(expressionLevels2))
    # corTmp <- tryCatch(cor(x,y,method="s",use="c"),error=function(e)NaN)
    
    outList <- append(outList,corTmp)
    
    if(is.finite(corTmp))
    {
      par(mar = c(5,5,5,5) + .1)
      smoothScatter(x,y,postPlotHook = fudgeit
                    ,xlab="Replicate 1",ylab="Replicate 2",main=i,xlim=xSatBounds,ylim=ySatBounds)
      abline(0,1,col=2,lwd=2)
      points(xSat,ySat,pch=16,col=2)
      tryCatch(points(lowess(x,y),col=3,lwd=2,type="l"),error=function(e){return(NULL)})
      legend("bottomright",legend=paste0("Cor.=",round(corTmp,2)),bty = "n",text.font=2)
    }
  }
  dev.off()
  names(outList) <- colnames(expressionLevels1)[!grepl("0$",colnames(expressionLevels1))]
  outList
}

replicatesRatesCorrelations <- function(object1,object2,width=7,height=12,name="",lowSat=0,upSat=1,ratesOrder=c("k1","k2","k3","k4","k5","k11","k9","k6","k8","k7","k10"))
{
  inferredRates1 <- object1$inferedRates
  inferredRates2 <- object2$inferedRates
  
  commonGenes <- intersect(rownames(object1$inferedRates),rownames(object2$inferedRates))
  
  outList <- list()
  
  RateNames <- c("Synthesis\n[fg/(M. cells*h)]"
                 ,"Co-trascriptional\nsplicing\n[h^-1]"
                 ,"Detachment of\nmature chromatin\n[h^-1]"
                 ,"Detachment of\npremature chromatin\n[h^-1]"
                 ,"Post-transcriptional\nsplicing\n[h^-1]"
                 ,"Export\n[h^-1]"
                 ,"Polysomal association\n[h^-1]"
                 ,"Cytoplasmatic\ndegradation\n[h^-1]"
                 ,"Nucleoplasmic mature\ndegradation\n[h^-1]"
                 ,"Polysomal degradation\n[h^-1]"
                 ,"Nucleoplasmic premature\ndegradation\n[h^-1]")
  names(RateNames) <- paste0("k",1:11)
  
  ratesSorted <- colnames(inferredRates1)
  names(ratesSorted) <- ratesSorted
  
  ratesSorted <- ratesSorted[ratesOrder[ratesOrder%in%ratesSorted]]
  
  pdf(paste0("ratesCorrelations_",name,".pdf"),width=10,height=13)
  par(mfrow=c(4,3))
  for(i in ratesSorted)
  {
    corTmp <- tryCatch(cor(inferredRates1[commonGenes,i],inferredRates2[commonGenes,i],method="s",use="c"),error=function(e)NaN)
    
    x <- log10(inferredRates1[commonGenes,i])
    y <- log10(inferredRates2[commonGenes,i])
    
    xSat <- x
    ySat <- y
    
    xSatBounds <- quantile(xSat,c(lowSat,upSat))
    ySatBounds <- quantile(ySat,c(lowSat,upSat))
    print(c(xSatBounds,ySatBounds))
    
    xSat <- xSat[xSat<xSatBounds[[1]]|xSat>xSatBounds[[2]]|ySat<ySatBounds[[1]]|ySat>ySatBounds[[2]]]
    ySat <- ySat[names(xSat)]
    
    xSat[xSat<xSatBounds[[1]]] <- xSatBounds[[1]]
    xSat[xSat>xSatBounds[[2]]] <- xSatBounds[[2]]
    
    ySat[ySat<ySatBounds[[1]]] <- ySatBounds[[1]]
    ySat[ySat>ySatBounds[[2]]] <- ySatBounds[[2]]
    
    commonGenes <- intersect(rownames(inferredRates1),rownames(inferredRates2))
    # corTmp <- tryCatch(cor(x,y,method="s",use="c"),error=function(e)NaN)
    
    outList <- append(outList,corTmp)
    
    xSat5 <- x
    ySat5 <- y
    
    xSat5Bounds <- quantile(xSat5,c(0.1,0.9))
    ySat5Bounds <- quantile(ySat5,c(0.1,0.9))
    print(c(xSat5Bounds,ySat5Bounds))
    
    xSat5 <- xSat5[xSat5<xSat5Bounds[[1]]|xSat5>xSat5Bounds[[2]]|ySat5<ySat5Bounds[[1]]|ySat5>ySat5Bounds[[2]]]
    ySat5 <- ySat5[names(xSat5)]
    
    lw_x <- x[setdiff(names(x),names(xSat5))]
    lw_y <- y[setdiff(names(y),names(ySat5))]
    
    if(is.finite(corTmp))
    {
      par(mar = c(5,5,5,7) + .1)
      smoothScatter(x,y,postPlotHook = fudgeit
                    ,xlab="Replicate 1",ylab="Replicate 2",main=RateNames[[i]],xlim=xSatBounds,ylim=ySatBounds)
      abline(0,1,col=2,lwd=2)
      points(xSat,ySat,pch=16,col=2)
      tryCatch(points(lowess(lw_x,lw_y),col=3,lwd=3,type="l"),error=function(e){return(NULL)})
      legend("bottomright",legend=paste0("Cor.=",round(corTmp,2)),bty = "n",text.font=2)
    }
  }
  dev.off()
  names(outList) <- colnames(inferredRates1)[!grepl("0$",colnames(inferredRates1))]
  outList
}

simulatedRatesCorrelations <- function(object1,simObject,width=20,height=12,imageStructure=NULL,name=""
                                       ,lowSat=0,upSat=1,ratesOrder=c("k1","k2","k3","k4","k5","k11","k9","k6","k8","k7","k10"))
{
  inferredRates1 <- object1$inferedRates
  inferredRates2 <- simObject$exampleRates
  
  commonGenes <- intersect(rownames(inferredRates1),rownames(inferredRates2))
  
  outList <- list()
  
  RateNames <- c("Synthesis\n[fg/(M. cells*h)]"
                 ,"Co-trascriptional\nsplicing\n[h^-1]"
                 ,"Detachment of\nmature chromatin\n[h^-1]"
                 ,"Detachment of\npremature chromatin\n[h^-1]"
                 ,"Post-transcriptional\nsplicing\n[h^-1]"
                 ,"Export\n[h^-1]"
                 ,"Polysomal association\n[h^-1]"
                 ,"Cytoplasmatic\ndegradation\n[h^-1]"
                 ,"Nucleoplasmic mature\ndegradation\n[h^-1]"
                 ,"Polysomal degradation\n[h^-1]"
                 ,"Nucleoplasmic premature\ndegradation\n[h^-1]")
  names(RateNames) <- paste0("k",1:11)
  
  ratesSorted <- colnames(inferredRates1)
  names(ratesSorted) <- ratesSorted
  
  ratesSorted <- ratesSorted[ratesOrder[ratesOrder%in%ratesSorted]]
  
  pdf(paste0("simRatesCorrelations_",name,".pdf"),width=width,height=height)
  if(is.null(imageStructure)){par(mfrow=c(1,length(ratesSorted)))}else{par(mfrow=imageStructure)}
  for(i in ratesSorted)
  {
    corTmp <- tryCatch(cor(inferredRates1[commonGenes,i],inferredRates2[commonGenes,i],method="s",use="c"),error=function(e)NaN)
    
    x <- log10(inferredRates1[commonGenes,i])
    y <- log10(inferredRates2[commonGenes,i])
    
    xSat <- x
    ySat <- y
    
    xSatBounds <- quantile(xSat,c(lowSat,upSat))
    ySatBounds <- quantile(ySat,c(lowSat,upSat))
    
    xSat <- xSat[xSat<xSatBounds[[1]]|xSat>xSatBounds[[2]]|ySat<ySatBounds[[1]]|ySat>ySatBounds[[2]]]
    ySat <- ySat[names(xSat)]
    
    xSat[xSat<xSatBounds[[1]]] <- xSatBounds[[1]]
    xSat[xSat>xSatBounds[[2]]] <- xSatBounds[[2]]
    
    ySat[ySat<ySatBounds[[1]]] <- ySatBounds[[1]]
    ySat[ySat>ySatBounds[[2]]] <- ySatBounds[[2]]
    
    commonGenes <- intersect(rownames(inferredRates1),rownames(inferredRates2))
    # corTmp <- tryCatch(cor(x,y,method="s",use="c"),error=function(e)NaN)
    
    outList <- append(outList,corTmp)
    
    if(is.finite(corTmp))
    {
      par(mar = c(5,4,4,5) + .1,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
      smoothScatter(x,y,postPlotHook = fudgeit
                    ,xlab="Inferred",ylab="Expected",main=RateNames[[i]],xlim=xSatBounds,ylim=ySatBounds)
      abline(0,1,col=2,lwd=2)
      points(xSat,ySat,pch=16,col=2)
      tryCatch(points(lowess(x,y),col=3,lwd=2,type="l"),error=function(e){return(NULL)})
      legend("bottomright",legend=paste0("Cor.=",round(corTmp,2)),bty = "n",text.font=2,cex=1.5)
    }
  }
  dev.off()
  names(outList) <- ratesSorted
  outList
}

differentialSpecies <- function(CIControl,CITreatment,nameTmp=NULL,disjointed=FALSE)
{
  commonGenesTmp <- intersect(rownames(CIControl$mean),rownames(CITreatment$mean))
  
  CIControlBottom <- CIControl$bottom
  CIControlUp <- CIControl$up
  
  CITreatmentBottom <- CITreatment$bottom
  CITreatmentUp <- CITreatment$up
  
  CIControlBottom[CIControl$mean==1e-10] <- NaN
  CIControlUp[CIControl$mean==1e-10] <- NaN
  
  CITreatmentBottom[CITreatment$mean==1e-10] <- NaN
  CITreatmentUp[CITreatment$mean==1e-10] <- NaN
  
  CIControl$mean[CIControl$mean==1e-10] <- NaN
  CITreatment$mean[CITreatment$mean==1e-10] <- NaN
  
  CIControlMean <- CIControl$mean[commonGenesTmp,]
  CITreatmentMean <- CITreatment$mean[commonGenesTmp,]
  
  CIControlBottom <- CIControlBottom[commonGenesTmp,]
  CIControlUp <- CIControlUp[commonGenesTmp,]
  
  CITreatmentBottom <- CITreatmentBottom[commonGenesTmp,]
  CITreatmentUp <- CITreatmentUp[commonGenesTmp,]
  
  overlapsTmp <- matrix(0,nrow=nrow(CIControlMean),ncol=ncol(CIControlMean))
  rownames(overlapsTmp) <- rownames(CIControlMean)
  colnames(overlapsTmp) <- colnames(CIControlMean)
  
  CIControlBottom[!is.finite(CIControlMean)] <- NaN
  CIControlUp[!is.finite(CIControlMean)] <- NaN
  
  CITreatmentBottom[!is.finite(CITreatmentMean)] <- NaN
  CITreatmentUp[!is.finite(CITreatmentMean)] <- NaN
  
  if(disjointed)
  {
    overlapsTmp[CIControlBottom>CITreatmentUp] <- (-1)
    overlapsTmp[CIControlUp<CITreatmentBottom] <- 1
  }else{
    overlapsTmp <- (((CITreatmentMean<CIControlBottom)|(CITreatmentMean>CIControlUp))
                    &((CIControlMean<CITreatmentBottom)|(CIControlMean>CITreatmentUp)))
    
    overlapsTmp[(CITreatmentMean<CIControlBottom)&(CIControlMean>CITreatmentUp)] <- (-1)
    overlapsTmp[(CIControlMean<CITreatmentBottom)&(CITreatmentMean>CIControlUp)] <- 1
  }
  
  overlapsTmp <- overlapsTmp[,!grepl("n0$",colnames(overlapsTmp))]
  # overlapsTmp <- overlapsTmp[,!grepl("t$",colnames(overlapsTmp))]
  
  if(!is.null(nameTmp))
  {
    barplot(apply(abs(overlapsTmp),2,function(i)sum(is.finite(i))),ylim=c(0,2000),las=2,main=paste0("Number of finite genes\n",nameTmp),ylab="")
    barplot(apply(abs(overlapsTmp),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i))),ylim=c(0,1),las=2,main=paste0("Fraction of differential genes\n",nameTmp),ylab="")
    x <- table(apply(abs(overlapsTmp),1,sum,na.rm=TRUE))
    y <- rep(0,length(0:max(as.numeric(names(x)))))
    names(y) <- 0:max(as.numeric(names(x)))
    y[names(x)] <- x
    barplot(y,las=2,main=paste0("Number of differential spices\n",nameTmp),ylab="")    
  }
  
  overlapsTmp
}

## Plots relative to couplings

correlationSignificance <- function(matTmp,refMat,method="s")
{ 
  refMat <- refMat$inferedRates
  matTmp <- matTmp$inferedRates
  
  if(ncol(matTmp)==9){
    colnames(matTmp) <- colnames(refMat) <- c("k1","k2","k3","k4","k5","k6","k8","k7","k9")
    matTmp <- matTmp[,sort(colnames(matTmp))]
    refMat <- refMat[,sort(colnames(refMat))]
  }else{
    colnames(matTmp) <- colnames(refMat) <- c("k1","k2","k3","k4","k5","k6","k7")
    matTmp <- matTmp[,sort(colnames(matTmp))]
    refMat <- refMat[,sort(colnames(refMat))]
  }
  
  common_genes <- intersect(rownames(refMat),rownames(matTmp))
  matTmp <- log2(matTmp[common_genes,]/refMat[common_genes,colnames(matTmp)])
  combi <- apply(combn(colnames(matTmp),2),2,paste,collapse='_')
  
  CorTmp <- cor(matTmp,method=method)
  RatesCorTmp <- CorTmp[t(upper.tri(CorTmp, diag = FALSE))]; names(RatesCorTmp) <- combi
  
  RatesCorSigTmp <- sapply(names(RatesCorTmp),function(i)
  {
    x <- strsplit(i,"_")[[1]]
    cor.test(matTmp[,x[[1]]],matTmp[,x[[2]]],method=method)$p.value
  })
  
  list(cor=RatesCorTmp,sig=RatesCorSigTmp,significant=names(RatesCorSigTmp[RatesCorSigTmp<(1e-4)]))
}

Mapping <- function(mat,refMat,thr=NULL,name,nclust,plot,annotation=NULL)
{
  
  cor <- correlationSignificance(mat,refMat,method='s')
  
  #setting names and taking common genes
  refMat <- refMat$inferedRates
  mat <- mat$inferedRates
  
  if(ncol(mat)==9){
    colnames(mat) <- colnames(refMat) <- c("k1","k2","k3","k4","k5","k6","k8","k7","k9")
    mat <- mat[,sort(colnames(mat))]
    refMat <- refMat[,sort(colnames(refMat))]
  }else{
    colnames(mat) <- colnames(refMat) <- c("k1","k2","k3","k4","k5","k6","k7")
    mat <- mat[,sort(colnames(mat))]
    refMat <- refMat[,sort(colnames(refMat))]
  }
  
  common_genes <- intersect(rownames(refMat),rownames(mat))
  #FC matrix
  mat <- log2(mat[common_genes,]/refMat[common_genes,colnames(mat)])
  combi <- apply(combn(colnames(mat),2),2,paste,collapse='_')
  
  matTmp <- mat 
  
  if(!is.null(thr)){
    matTmp[matTmp<thr&matTmp>(-thr)] = 0
    matTmp[matTmp>0] = 2
    matTmp[matTmp<0] = 1
    out <- sapply(combi,function(i){
      x <- strsplit(i,"_")[[1]]
      p <- matTmp[,x[1]]*matTmp[,x[2]]
    })
  }else{
    matTmp[matTmp>0] = 2
    matTmp[matTmp<0] = 1
    out <- sapply(combi,function(i){
      x <- strsplit(i,"_")[[1]]
      p <- matTmp[,x[1]]*matTmp[,x[2]]
    })
  }
  
  out <- out[,names(which(cor$sig<(1e-4)))]
  out <- out[apply(out!=0,1,any),]
  
  out[out==2] <- (-1)
  out[out==4] <- (1.25)
  out[out==1] <- (0.75)
  
  rowDist <- hclust(as.dist(hamming.distance(out)),method="complete")
  colDist <- hclust(as.dist(hamming.distance(t(out))),method="complete")
  
  foe <- pheatmap(out,cluster_cols=colDist,
                  cluster_rows=rowDist,
                  show_rownames=FALSE,
                  color=c("navyblue","beige","darkolivegreen3","gold"),
                  breaks=c(-1.5,-0.5,0.5,1,1.5),
                  cutree_rows=nclust,silent=TRUE)
  
  ord <- foe$tree_row$order
  ord2 <- cutree(foe$tree_row,nclust)[ord]
  
  ord2 <- split(ord2,ord2)[unique(ord2)]
  
  ord2 <- unlist(sapply(seq_along(ord2),function(i){rep(LETTERS[i],length(ord2[[i]]))}))
  names(ord2) <- rownames(out)[ord]
  
  rowAnnotation <- as.data.frame(list("Cluster"=ord2))
  rowAnnotationColors <- carto_pal(length(unique(rowAnnotation[,"Cluster"])), "Safe")
  names(rowAnnotationColors) <- unique(rowAnnotation[,"Cluster"])
  
  if(!is.null(annotation)){
    rowAnnotation <- cbind(rowAnnotation,annotation)
    annCol <- list(Hamming = c("-"="blue","+"="red","0"="white"))
  }
  
  
  
  if(plot){
    pdf(file=paste0("HeatmapBinary_",name,".pdf"))
    foe <- pheatmap(out[ord,],cluster_cols=FALSE,
                    cluster_rows=FALSE,
                    show_rownames=FALSE,
                    color=c("navyblue","beige","darkolivegreen3","gold"),
                    breaks=c(-1.5,-0.5,0.5,1,1.5),
                    cutree_rows=nclust,
                    annotation_row=rowAnnotation,
                    annotation_colors=list("Cluster"=rowAnnotationColors))
    dev.off()
  }
  return(list(out,mat,cor,ord2))
}

PlotNetwork <- function(matTmp,refMat,name,plot)
{
  if(is.list(matTmp[[1]]))
  {
    mat_corr <- lapply(matTmp,function(i)
    {
      outTmp <- correlationSignificance(i,refMat,method="s")
      outTmp$cor[outTmp$sig>1e-4] <- 0
      outTmp
    })
 
    commonEdgesTmp <- names(which(table(unlist(lapply(mat_corr,function(i)i$significant)))==length(matTmp)))
    corTmp <- apply(sapply(mat_corr,function(i)i$cor[commonEdgesTmp]),1,mean)
    sigTmp <- apply(sapply(mat_corr,function(i)i$sig[commonEdgesTmp]),1,mean)
 
    mat_corr  <- list(cor=corTmp,sig=sigTmp,significant=commonEdgesTmp)
 
  }else{
    mat_corr <- correlationSignificance(matTmp,refMat,method="s")
    mat_corr$cor[mat_corr$sig>1e-4] <- 0    
  }
 
  if(plot){
    pdf(file=paste0(name,"_network.pdf"),width=4,height=4)
    qwe <- graph.data.frame(data.frame("from"=sapply(strsplit(names(mat_corr$cor),"_"),"[[",1)
                                            ,"to"=sapply(strsplit(names(mat_corr$cor),"_"),"[[",2))
                                  ,directed=FALSE,vertices=paste0("k",1:9))
          igraph.options(plot.layout=layout.circle)
          colTmp <- c("blue","white","red")[sign(mat_corr$cor)+2]
          plot.igraph(qwe,edge.width=abs(mat_corr$cor)^3*15,edge.color=colTmp,vertex.size=35,main=name)
    dev.off()
  }
 
  return(mat_corr)
}

PlotBar <- function(mat,name){
  
  pdf(file=paste0("Barplot_",name,".pdf"),width=15)
  barplot(mat,col=c("beige","gold","darkolivegreen3","navyblue") , 
          border="white", 
          space=0.04, 
          font.axis=2, 
          xlab="edge",
          las=2)
  dev.off()
}

EnrichFunc <- function(cluster,binaryMat,name)
{
  
  GO <- lapply(cluster,function(i){
    
    enrichGO(gene = names(i),
             OrgDb = org.Hs.eg.db,
             keyType = 'ENTREZID',
             ont = "ALL",
             pAdjustMethod = "BH",
             universe = rownames(binaryMat),
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.1,
             readable = TRUE,
             maxGSSize = 3000)
    
  })
  GO <- GO[sapply(GO,function(i){nrow(i@result)!=0})]
  GO <- lapply(GO,function(i){cbind(geneID=i$geneID,pAdjusted=signif(i$p.adjust,digits=2),Description=i$Description)})
  saveRDS(file=paste0(name,"_GO.rds"),GO)
}

NullModelHD <- function(mat1,mat2,nShuffle,method)
{
  if(method==1){
    
    qwe <- sapply(1:nShuffle,function(k){
      set.seed(k)
      mat1 <- mat1[sample(nrow(mat1)),]
      mat2 <- mat2[sample(nrow(mat2)),]
      foe <- sapply(1:nrow(mat1),function(i){hamming.distance(mat1[i,],mat2[i,])})
    })
    
  }else{
    
    rownames(mat1) <- 1:nrow(mat1)
    rownames(mat2) <- (nrow(mat2)+1):(nrow(mat2)*2)
    
    newMat1 <- rbind(mat1[sample(floor(nrow(mat1)/2)),],mat2[sample(floor(nrow(mat1)/2)),])
    newMat2 <- rbind(mat1,mat2)
    newMat2 <- newMat2[setdiff(rownames(newMat2),rownames(newMat1)),]
    
    qwe <- sapply(1:nShuffle,function(k){
      set.seed(k)
      newMat1 <- newMat1[sample(nrow(newMat1)),]
      newMat2 <- newMat2[sample(nrow(newMat2)),]
      foe <- sapply(1:nrow(newMat1),function(i){hamming.distance(newMat1[i,],newMat2[i,])})
    })
  }  
}

FoldChanges <- function(matrix)
{
  combi <- apply(combn(colnames(matrix),2),2,paste,collapse='_')
  fc_edges <- sapply(combi,function(i){
    x <- strsplit(i,"_")[[1]]
    matrix[,x[1]]*matrix[,x[2]]
  })
  fc_edges
}

PlotBars <- function(rbpMat,tfMat,name,significant)
{  
  rbpMat <- lapply(rbpMat,function(i){cbind("ID"=i$ID,"setSize"=i$setSize,"pVal"=i$p.adjust,"type"=rep("skyblue",length(i$ID)))})
  tfMat <- lapply(tfMat,function(i){cbind("ID"=i$ID,"setSize"=as.numeric(i$setSize),"pVal"=i$p.adjust,"type"=rep("gray",length(i$ID)))})
  rbpMat[intersect(names(tfMat),names(rbpMat))] <- lapply(intersect(names(tfMat),names(rbpMat)),function(k)
  {
    qwe <- rbind(rbpMat[[k]],tfMat[[k]])
    qwe[order(qwe[,"pVal"],decreasing=FALSE),]
  })
  dfMat <- c(rbpMat,tfMat[!names(tfMat)%in%names(rbpMat)])
  dfMat <- dfMat[names(dfMat)%in%significant]
  dfMat <- lapply(names(dfMat),function(i){
    dfMat[[i]][,"ID"] <- sapply(strsplit(dfMat[[i]][,1],"-"),'[',1)
    dfMat[[i]][,"ID"] <- sapply(dfMat[[i]][,"ID"],function(k){paste0(k,"-",i)})
    dfMat[[i]][,"pVal"] <- as.numeric(signif(as.numeric(dfMat[[i]][,"pVal"]),digits=2))
    if(nrow(dfMat[[i]])>5){
      dfMat[[i]] <- dfMat[[i]][1:5,]
    }
    dfMat[[i]]
  })
  dfMat <- do.call(rbind,dfMat)
  
  pdf(file=paste0("Barplot_",name,".pdf"),width=20,height=5)
  par(mar = c(7, 5, 5, 7),xpd=TRUE)
  
  foe <- barplot(as.numeric(dfMat[,"setSize"]),
                 names.arg=dfMat[,"ID"],
                 las=2,
                 col=dfMat[,"type"],
                 ylab="Set size",
                 xlab=NULL,
                 cex.lab=1.3,
                 space=1)
  text(foe,as.numeric(dfMat[,"setSize"])+40,labels=signif(as.numeric(dfMat[,"pVal"]),2),xpd=TRUE,srt=270,cex=0.7)
  legend("topright",legend=c("RBP","TF"),fill=c("skyblue","gray"),horiz=FALSE,inset=c(-0.2, 0))
  
  dev.off()
}

GSEA_table <- function(rbpMat,tfMat,name)
{
  rbpMat <- lapply(rbpMat,function(i){i$ID})
  rbpMat <- lapply(rbpMat,function(i){sapply(strsplit(i,"-"),'[',1)})
  tfMat <- lapply(tfMat,function(i){i$ID})
  tfMat <- lapply(tfMat,function(i){sapply(strsplit(i,"-"),'[',1)})
  df <- data.frame(cbind("RBPs"=rep(0,length(union(names(rbpMat),names(tfMat)))),"TFs"=rep(0,length(union(names(rbpMat),names(tfMat))))))
  rownames(df) <- sort(union(names(rbpMat),names(tfMat)))
  
  RBP <- sapply(rownames(df),function(i){
    
    if(i%in%names(rbpMat)){
      df[i,"RBPs"] <- list(paste(rbpMat[[i]], collapse=","))
    }else{
      df[i,"RBPs"] <- ""
    }
  })
  TF <- sapply(rownames(df),function(i){
    
    if(i%in%names(tfMat)){
      df[i,"TFs"] <- list(paste(tfMat[[i]], collapse=","))
    }else{
      df[i,"TFs"] <- ""
    }
  })
  
  df$RBPs <- RBP
  df$TFs <- TF
  saveRDS(file=paste0("df",name,".rds"),df)
}

NumberOfCouplings <- function(mat,name){

    zero <- apply(mat,1,function(i){length(i[i!=0])})
    zero <- table(zero)
    vec <- rep(0,ncol(mat))
    names(vec) <- 0:(ncol(mat)-1)
    vec[names(zero)] <- zero
    
    pdf(paste0("UnregulatedEdges_",name,".pdf"))
    barplot(vec,xlab="Edges",ylab="Genes",las=2,space=0.3)
    dev.off()
}


PlotGO <- function(file,name){

  df <- lapply(names(file),function(i)
  {
    if(nrow(file[[i]])>1){
      foe <- cbind(file[[i]][,2:3],"Cluster"=rep(i,nrow(file[[i]])))
      rownames(foe) <- foe[,"Description"]
    foe
    }else{
      foe <- t(as.matrix(c(file[[i]][,2:3],"Cluster"=i),nrow=1))
      rownames(foe) <-  foe[2]
      foe
    }
    
  })
  names(df) <- names(file)

  lapply(names(df),function(k){

    pdf(paste0("GO",name,"_",k,".pdf"),width=8,height=5)
    par(mar=c(5,20,2,2))
      barplot(as.numeric(df[[k]][,1]),horiz=TRUE,xlab="Adjusted p-value",names.arg = rownames(df[[k]]),las=2)
    dev.off()
  })
}
