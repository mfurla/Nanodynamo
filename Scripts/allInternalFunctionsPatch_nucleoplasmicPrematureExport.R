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
      # print("Counts threshold impact:")
      # print(table(selectionTmp))
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
      # print("Counts threshold impact:")
      # print(table(selectionTmp))
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
                                     ,expDesign # Experimental design.
                                     ,tbl # Matrix of sequencing statistics.
                                     ,labelingTime # Labeling time of the fractionation samples.
                                     ,labelingTimePoly=NULL # Labeling time of the polysomal samples.
                                     ,txdb # Annotation object.
                                     ,countsTh = TRUE # Select genes with at least on read in all the fractions of interest.
                                     ,cytPTh # DEPRECATED - TO BE REMOVED.
                                     ,minoverlap_I=10 # Minimum overlap for intronic reads.
                                     ,minoverlap_E=10 # Minimum overlap for exonic reads.
                                     ,spikeInsConcentrations=NULL # DEPRECATED - TO BE REMOVED.
                                     ,yeastCounts=NULL # Yeast counts for samples not aligned with ENO2 sequence; TO BE REMOVED.
                                     ,mergeSamples=FALSE) # If TRUE merge replicates counts instead of mediate them
{
  if(!identical(sort(names(expDesign)),sort(names(bamPaths)))|
     !identical(sort(names(expDesign)),sort(names(nascentPaths))))
  {
    print("Issue in the experimental expDesign.")
    return(NULL)
  }
  
  # Replicates
  samples <- paste0("Sample_",unique(sapply(strsplit(names(expDesign),"_"),"[[",2)))
  
  print("Detected samples:")
  print(samples)
  
  # Giving names to samples
  chrSamples <- names(expDesign)[grepl("chr",expDesign)]
  nucSamples <- names(expDesign)[grepl("nuc",expDesign)]
  cytSamples <- names(expDesign)[grepl("cyt",expDesign)]
  polySamples <- names(expDesign)[grepl("poly",expDesign)]
  
  if(length(polySamples)>0&!is.null(labelingTimePoly))
  {
    print("Polysomal RNA mode.")
  }else if(length(polySamples)==0&is.null(labelingTimePoly)){
    print("No Polysomal RNA mode.")
    polysomal = NULL
    polysomalP = NULL
    normalizedCountsPolysomal = NULL
  }else{
    print("Check polysomal labeling time and/or experimental design.")
    stop()
  }
  
  # Estimate counts for each chromatin species
  chromatin <- lapply(chrSamples,function(i)
  {
    bamPathTmp <- bamPaths[[i]]
    nascentPathTmp <- nascentPaths[[i]]
    
    genesCountsProfiling(bamPath=bamPathTmp
                         ,nascentPath=nascentPathTmp
                         ,label="ch"
                         ,labelingTime=labelingTime
                         ,prematureFlag=TRUE
                         ,countsTh=0
                         ,txdb=txdb
                         ,minoverlap_I=minoverlap_I
                         ,minoverlap_E=minoverlap_E)
  })
  names(chromatin) <- chrSamples
  
  # Estimate counts for each nucleoplasmic species
  nucleoplasmic <- lapply(nucSamples,function(i)
  {
    bamPathTmp <- bamPaths[[i]]
    nascentPathTmp <- nascentPaths[[i]]
    
    genesCountsProfiling(bamPath=bamPathTmp
                         ,nascentPath=nascentPathTmp
                         ,label="n"
                         ,labelingTime=labelingTime
                         ,prematureFlag=TRUE
                         ,countsTh=0
                         ,txdb=txdb
                         ,minoverlap_I=minoverlap_I
                         ,minoverlap_E=minoverlap_E)
  })
  names(nucleoplasmic) <- nucSamples
  
  # Estimate counts for each cytoplasmic species - not considering premature cytoplasmic
  cytoplasmic <- lapply(cytSamples,function(i)
  {
    bamPathTmp <- bamPaths[[i]]
    nascentPathTmp <- nascentPaths[[i]]
    
    genesCountsProfiling(bamPath=bamPathTmp
                         ,nascentPath=nascentPathTmp
                         ,label="cy"
                         ,labelingTime=labelingTime
                         ,prematureFlag=FALSE           #cytoplasmic without premature
                         ,countsTh=0
                         ,txdb=txdb
                         ,minoverlap_I=minoverlap_I
                         ,minoverlap_E=minoverlap_E)
  })
  names(cytoplasmic) <- cytSamples
  
  # Estimate counts for each cytoplasmic species - considering premature cytoplasmic
  cytoplasmicP <- lapply(cytSamples,function(i)
  {
    bamPathTmp <- bamPaths[[i]]
    nascentPathTmp <- nascentPaths[[i]]
    
    genesCountsProfiling(bamPath=bamPathTmp
                         ,nascentPath=nascentPathTmp
                         ,label="cy"
                         ,labelingTime=labelingTime
                         ,prematureFlag=TRUE
                         ,countsTh=0
                         ,txdb=txdb
                         ,minoverlap_I=minoverlap_I
                         ,minoverlap_E=minoverlap_E)
  })
  names(cytoplasmicP) <- cytSamples
  
  if(length(polySamples)>0){
    polysomal <- lapply(polySamples,function(i)
    {
      bamPathTmp <- bamPaths[[i]]
      nascentPathTmp <- nascentPaths[[i]]
      
      genesCountsProfiling(bamPath=bamPathTmp
                           ,nascentPath=nascentPathTmp
                           ,label="py"
                           ,labelingTime=labelingTimePoly
                           ,prematureFlag=FALSE           #polysomal without premature
                           ,countsTh=0
                           ,txdb=txdb
                           ,minoverlap_I=minoverlap_I
                           ,minoverlap_E=minoverlap_E)
    })
    names(polysomal) <- polySamples
    
    polysomalP <- lapply(polySamples,function(i)
    {
      bamPathTmp <- bamPaths[[i]]
      nascentPathTmp <- nascentPaths[[i]]
      
      genesCountsProfiling(bamPath=bamPathTmp
                           ,nascentPath=nascentPathTmp
                           ,label="py"
                           ,labelingTime=labelingTimePoly
                           ,prematureFlag=TRUE           #polysomal with premature
                           ,countsTh=0
                           ,txdb=txdb
                           ,minoverlap_I=minoverlap_I
                           ,minoverlap_E=minoverlap_E)
    })
    names(polysomalP) <- polySamples
  }
  
  if(mergeSamples)
  {
    allGenesChromatin <- unique(unlist(sapply(chromatin,rownames)))
    mergedChromatin <- list()
    mergedChromatin[[1]] <- data.frame(matrix(0,nrow=length(allGenesChromatin),ncol=ncol(chromatin[[1]])))
    rownames(mergedChromatin[[1]]) <- allGenesChromatin
    colnames(mergedChromatin[[1]]) <- colnames(chromatin[[1]])
    
    for(i in chromatin)
    {
      mergedChromatin[[1]][rownames(i),colnames(i)] <- mergedChromatin[[1]][rownames(i),colnames(i)]+i
    }
    
    allGenesNucleoplasmic <- unique(unlist(sapply(nucleoplasmic,rownames)))
    mergedNucleoplasmic <- list()
    mergedNucleoplasmic[[1]] <- data.frame(matrix(0,nrow=length(allGenesNucleoplasmic),ncol=ncol(nucleoplasmic[[1]])))
    rownames(mergedNucleoplasmic[[1]]) <- allGenesNucleoplasmic
    colnames(mergedNucleoplasmic[[1]]) <- colnames(nucleoplasmic[[1]])
    
    for(i in nucleoplasmic)
    {
      mergedNucleoplasmic[[1]][rownames(i),colnames(i)] <- mergedNucleoplasmic[[1]][rownames(i),colnames(i)]+i
    }
    
    allGenesCytoplasmic <- unique(unlist(sapply(cytoplasmic,rownames)))
    mergedCytoplasmic <- list()
    mergedCytoplasmic[[1]] <- data.frame(matrix(0,nrow=length(allGenesCytoplasmic),ncol=ncol(cytoplasmic[[1]])))
    rownames(mergedCytoplasmic[[1]]) <- allGenesCytoplasmic
    colnames(mergedCytoplasmic[[1]]) <- colnames(cytoplasmic[[1]])
    
    for(i in cytoplasmic)
    {
      mergedCytoplasmic[[1]][rownames(i),colnames(i)] <- mergedCytoplasmic[[1]][rownames(i),colnames(i)]+i
    }
    
    allGenesCytoplasmicP <- unique(unlist(sapply(cytoplasmicP,rownames)))
    mergedCytoplasmicP <- list()
    mergedCytoplasmicP[[1]] <- data.frame(matrix(0,nrow=length(allGenesCytoplasmicP),ncol=ncol(cytoplasmicP[[1]])))
    rownames(mergedCytoplasmicP[[1]]) <- allGenesCytoplasmicP
    colnames(mergedCytoplasmicP[[1]]) <- colnames(cytoplasmicP[[1]])
    
    for(i in cytoplasmicP)
    {
      mergedCytoplasmicP[[1]][rownames(i),colnames(i)] <- mergedCytoplasmicP[[1]][rownames(i),colnames(i)]+i
    }
    
    names(mergedChromatin) <- "Chr_Merged"
    names(mergedNucleoplasmic) <- "Nuc_Merged"
    names(mergedCytoplasmic) <- "Cyt_Merged"
    names(mergedCytoplasmicP) <- "Cyt_Merged"
    
    chromatin <- mergedChromatin
    nucleoplasmic <- mergedNucleoplasmic
    cytoplasmic <- mergedCytoplasmic
    cytoplasmicP <- mergedCytoplasmicP
    
    tblMerged <- rbind(apply(tbl[chrSamples,],2,sum)
                       ,apply(tbl[nucSamples,],2,sum)
                       ,apply(tbl[cytSamples,],2,sum))
    
    rownames(tblMerged) <- c("Chr_Merged","Nuc_Merged","Cyt_Merged")
    
    if(length(polySamples)>0){
      
      allGenesPolysomal <- unique(unlist(sapply(polysomal,rownames)))
      mergedPolysomal <- list()
      mergedPolysomal[[1]] <- data.frame(matrix(0,nrow=length(allGenesPolysomal),ncol=ncol(polysomal[[1]])))
      rownames(mergedPolysomal[[1]]) <- allGenesPolysomal
      colnames(mergedPolysomal[[1]]) <- colnames(polysomal[[1]])
      
      for(i in polysomal)
      {
        mergedPolysomal[[1]][rownames(i),colnames(i)] <- mergedPolysomal[[1]][rownames(i),colnames(i)]+i
      }
      
      allGenesPolysomalP <- unique(unlist(sapply(polysomalP,rownames)))
      mergedPolysomalP <- list()
      mergedPolysomalP[[1]] <- data.frame(matrix(0,nrow=length(allGenesPolysomalP),ncol=ncol(polysomalP[[1]])))
      rownames(mergedPolysomalP[[1]]) <- allGenesPolysomalP
      colnames(mergedPolysomalP[[1]]) <- colnames(polysomalP[[1]])
      
      for(i in polysomalP)
      {
        mergedPolysomalP[[1]][rownames(i),colnames(i)] <- mergedPolysomalP[[1]][rownames(i),colnames(i)]+i
      }
      
      names(mergedPolysomal) <- "Poly_Merged"
      names(mergedPolysomalP) <- "Poly_Merged"
      
      polysomal <- mergedPolysomal
      polysomalP <- mergedPolysomalP
      
      tblMerged <- rbind(tblMerged,apply(tbl[polySamples,],2,sum))
      rownames(tblMerged)[4] <- "Poly_Merged"
    }
    
    tbl <- tblMerged
    
    samples <- "Sample_Merged"
  }
  
  chromatinExpressedGenes <- lapply(chromatin,rownames)
  nucleoplasmicExpressedGenes <- lapply(nucleoplasmic,rownames)
  cytoplasmicExpressedGenes <- lapply(cytoplasmicP,rownames)
  polysomalExpressedGenes <- lapply(polysomalP,rownames)

# Genes with at least one read in all the fractions except for the polysomal
expressedGenes <- lapply(seq_along(cytoplasmicExpressedGenes),function(i)
{
  expressedGenes <- intersect(intersect(chromatinExpressedGenes[[i]]
                                        ,nucleoplasmicExpressedGenes[[i]])
                              ,cytoplasmicExpressedGenes[[i]])
  print(paste0("Expressed genes condition ",i,":",length(expressedGenes)))
  return(expressedGenes)      
})

# Expressed genes with chromatin premature RNA
if(countsTh)
{
  expressedGenesYesChp <- lapply(seq_along(chromatin),function(j)
  {
    i <- chromatin[[j]]
    expressedGenes[[j]][i[expressedGenes[[j]],paste0("chpn",labelingTime)]>0
                        &i[expressedGenes[[j]],paste0("chmn",labelingTime)]>0
                        &i[expressedGenes[[j]],paste0("chpp",labelingTime)]>0
                        &i[expressedGenes[[j]],paste0("chmp",labelingTime)]>0]
  })    
}else{
  expressedGenesYesChp <- lapply(seq_along(chromatin),function(j)
  {
    i <- chromatin[[j]]
    expressedGenes[[j]][i[expressedGenes[[j]],"chpt"]>0
                        &i[expressedGenes[[j]],"chmt"]>0]
  })
}

# Expressed genes with nucleoplasmic premature RNA
if(countsTh)
{
  expressedGenesYesNp <- lapply(seq_along(nucleoplasmic),function(j)
  {
    i <- nucleoplasmic[[j]]
    expressedGenes[[j]][i[expressedGenes[[j]],paste0("npn",labelingTime)]>0
                        &i[expressedGenes[[j]],paste0("nmn",labelingTime)]>0
                        &i[expressedGenes[[j]],paste0("npp",labelingTime)]>0
                        &i[expressedGenes[[j]],paste0("nmp",labelingTime)]>0]
  })
}else{
  expressedGenesYesNp <- lapply(seq_along(nucleoplasmic),function(j)
  {
    i <- nucleoplasmic[[j]]
    expressedGenes[[j]][i[expressedGenes[[j]],"npt"]>0
                        &i[expressedGenes[[j]],"nmt"]>0]
  })
}

if(length(polySamples)>0)
{
  if(countsTh)
  {
    expressedGenesYesP <- lapply(seq_along(polysomal),function(j)
    {
      i <- polysomal[[j]]
      rownames(i)[i[,paste0("pyn",labelingTimePoly)]>0&i[,paste0("pyp",labelingTimePoly)]>0]
    })
  }else{
    expressedGenesYesP <- lapply(seq_along(polysomal),function(j)
    {
      i <- polysomal[[j]]
      rownames(i)[i[,"pyt"]>0]
    })    
  }
}else{
  expressedGenesYesP <- NULL
}

# Expressed genes with chromatin premature RNA, nucleoplasmic premature RNA, and polysomal RNA
expressedGenesYesChpNpP <- lapply(seq_along(expressedGenesYesNp),function(j)
{
  expressedGenesYesChpNpP <- intersect(intersect(expressedGenesYesChp[[j]]
                                                 ,expressedGenesYesNp[[j]])
                                       ,expressedGenesYesP[[j]])
  intersect(expressedGenes[[j]],expressedGenesYesChpNpP)
})

# Expressed genes with chromatin premature RNA, nucleoplasmic premature RNA, and without polysomal RNA
expressedGenesYesChpNpNoP <- lapply(seq_along(expressedGenesYesNp),function(j)
{
  expressedGenesYesChpNp <- intersect(expressedGenesYesChp[[j]]
                                      ,expressedGenesYesNp[[j]])
  expressedGenesYesChpNpNoP <- setdiff(expressedGenesYesChpNp
                                       ,expressedGenesYesChpNpP[[j]])
  intersect(expressedGenes[[j]],expressedGenesYesChpNpNoP)
})

# Expressed genes with chromatin premature RNA, polysomal RNA, and without nucleoplasmic premature RNA
expressedGenesYesChpPNoNp <- lapply(seq_along(expressedGenesYesNp),function(j)
{
  expressedGenesYesChpP <- intersect(expressedGenesYesChp[[j]]
                                     ,expressedGenesYesP[[j]])
  expressedGenesYesChpPNoNp <- setdiff(expressedGenesYesChpP,expressedGenesYesNp[[j]])
  intersect(expressedGenes[[j]],expressedGenesYesChpPNoNp)
})

# Expressed genes with chromatin premature RNA, and without nucleoplasmic premature RNA and polysomal RNA
expressedGenesYesChpNoNpP <- lapply(seq_along(expressedGenesYesNp),function(j)
{
  expressedGenesYesChpNoP <- setdiff(expressedGenesYesChp[[j]],expressedGenesYesP[[j]])
  expressedGenesYesChpNoNp <- setdiff(expressedGenesYesChp[[j]],expressedGenesYesNp[[j]])
  
  expressedGenesYesChpNoNpP <- intersect(expressedGenesYesChpNoP,expressedGenesYesChpNoNp)
  intersect(expressedGenes[[j]],expressedGenesYesChpNoNpP)
})

# Expressed genes with polysomal RNA, and without chromatin premature RNA
expressedGenesYesPNoChp <- lapply(seq_along(expressedGenesYesNp),function(j)
{
  expressedGenesYesPNoChp <- setdiff(expressedGenesYesP[[j]],expressedGenesYesChp[[j]])
  intersect(expressedGenes[[j]],expressedGenesYesPNoChp)
})

# Expressed genes without polysomal RNA, and chromatin premature RNA
expressedGenesNoChpP <- lapply(seq_along(expressedGenesYesNp),function(j)
{
  expressedGenesNoChp <- setdiff(expressedGenes[[j]],expressedGenesYesChp[[j]])
  expressedGenesNoP <- setdiff(expressedGenes[[j]],expressedGenesYesP[[j]])
  expressedGenesNoChpP <- intersect(expressedGenesNoChp,expressedGenesNoP)
  intersect(expressedGenes[[j]],expressedGenesNoChpP)
})

# First set of outputs
outTmp1 <- list(chromatin=chromatin
                ,nucleoplasmic=nucleoplasmic
                ,cytoplasmic=cytoplasmic
                ,cytoplasmicP=cytoplasmicP
                ,polysomal = polysomal
                ,polysomalP = polysomalP
                ,expressedGenesYesChpNpP = expressedGenesYesChpNpP
                ,expressedGenesYesChpNpNoP = expressedGenesYesChpNpNoP
                ,expressedGenesYesChpPNoNp = expressedGenesYesChpPNoNp
                ,expressedGenesYesChpNoNpP = expressedGenesYesChpNoNpP
                ,expressedGenesYesPNoChp = expressedGenesYesPNoChp
                ,expressedGenesNoChpP = expressedGenesNoChpP)

# Spike-ins and Enolase2 quantification
chromatinSpikes <- lapply(chrSamples,function(i)
{
  print(i)
  bamPathTmp <- bamPaths[[i]]
  spikesQuantification(bamPathTmp,spikes=c("ERCC","ENO2","SIRV"))
})
names(chromatinSpikes) <- chrSamples

nucleoplasmicSpikes <- lapply(nucSamples,function(i)
{
  print(i)
  bamPathTmp <- bamPaths[[i]]
  spikesQuantification(bamPathTmp,spikes=c("ERCC","ENO2","SIRV"))
})
names(nucleoplasmicSpikes) <- nucSamples

cytoplasmicSpikes <- lapply(cytSamples,function(i)
{
  print(i)
  bamPathTmp <- bamPaths[[i]]
  spikesQuantification(bamPathTmp,spikes=c("ERCC","ENO2","SIRV"))
})
names(cytoplasmicSpikes) <- cytSamples

if(length(polySamples)>0)
{
  polysomalSpikes <- lapply(polySamples,function(i)
  {
    print(i)
    bamPathTmp <- bamPaths[[i]]
    spikesQuantification(bamPathTmp,spikes=c("ERCC","ENO2","SIRV"))
  })
  names(polysomalSpikes) <- polySamples
  
  allSpikes <- unique(c(unlist(sapply(chromatinSpikes,names))
                        ,unlist(sapply(nucleoplasmicSpikes,names))
                        ,unlist(sapply(cytoplasmicSpikes,names))
                        ,unlist(sapply(polysomalSpikes,names))))
}else{
  allSpikes <- unique(c(unlist(sapply(chromatinSpikes,names))
                        ,unlist(sapply(nucleoplasmicSpikes,names))
                        ,unlist(sapply(cytoplasmicSpikes,names))))    
}

# Spike-ins and total counts matrix
spikesCountsTmp <- matrix(0,nrow=length(allSpikes),ncol=(length(chrSamples)+length(nucSamples)+length(cytSamples)))

rownames(spikesCountsTmp) <- allSpikes
colnames(spikesCountsTmp) <- c(chrSamples,nucSamples,cytSamples)

for(i in chrSamples)
{
  spikesTmp <- chromatinSpikes[[i]]
  spikesCountsTmp[names(spikesTmp),i] <- spikesTmp
}
for(i in nucSamples)
{
  spikesTmp <- nucleoplasmicSpikes[[i]]
  spikesCountsTmp[names(spikesTmp),i] <- spikesTmp
}
for(i in cytSamples)
{
  spikesTmp <- cytoplasmicSpikes[[i]]
  spikesCountsTmp[names(spikesTmp),i] <- spikesTmp
}

if(mergeSamples)
{
  spikesCountsTmpMerged <- cbind(Chr_Merged=apply(spikesCountsTmp[,chrSamples],1,sum)
                                 ,Nuc_Merged=apply(spikesCountsTmp[,nucSamples],1,sum)
                                 ,Cyt_Merged=apply(spikesCountsTmp[,cytSamples],1,sum))    
}

if(length(polySamples)>0)
{
  spikesCountsTmp <- cbind(spikesCountsTmp,matrix(0,nrow=length(allSpikes),ncol=length(polySamples)))
  colnames(spikesCountsTmp) <- c(chrSamples,nucSamples,cytSamples,polySamples)
  
  for(i in polySamples)
  {
    spikesTmp <- polysomalSpikes[[i]]
    spikesCountsTmp[names(spikesTmp),i] <- spikesTmp
  } 
  
  if(mergeSamples)
  {
    spikesCountsTmpMerged <- cbind(spikesCountsTmpMerged,"Poly_Merged"=apply(spikesCountsTmp[,polySamples],1,sum))
  }
  
  expressedGenes <- unique(unlist(sapply(c(outTmp1[["chromatin"]],outTmp1[["nucleoplasmic"]],outTmp1[["cytoplasmic"]],outTmp1[["polysomal"]]),rownames)))
  
  totalCountsTmp <- cbind(sapply(outTmp1[["chromatin"]],function(i)
  {
    x <- i[expressedGenes,"chmt"];x[is.na(x)] <- 0
    y <- i[expressedGenes,"chpt"];y[is.na(y)] <- 0
    
    x+y
  })
  ,sapply(outTmp1[["nucleoplasmic"]],function(i)
  {
    x <- i[expressedGenes,"nmt"];x[is.na(x)] <- 0
    y <- i[expressedGenes,"npt"];y[is.na(y)] <- 0
    
    x+y
  })
  ,sapply(outTmp1[["cytoplasmic"]],function(i)
  {
    x <- i[expressedGenes,"cyt"];x[is.na(x)] <- 0
    x
  })
  ,sapply(outTmp1[["polysomal"]],function(i)
  {
    x <- i[expressedGenes,"pyt"];x[is.na(x)] <- 0
    x
  })
  )
  
}else{
  expressedGenes <- unique(unlist(sapply(c(outTmp1[["chromatin"]],outTmp1[["nucleoplasmic"]],outTmp1[["cytoplasmic"]]),rownames)))    
  
  totalCountsTmp <- cbind(sapply(outTmp1[["chromatin"]],function(i)
  {
    x <- i[expressedGenes,"chmt"];x[is.na(x)] <- 0
    y <- i[expressedGenes,"chpt"];y[is.na(y)] <- 0
    
    x+y
  })
  ,sapply(outTmp1[["nucleoplasmic"]],function(i)
  {
    x <- i[expressedGenes,"nmt"];x[is.na(x)] <- 0
    y <- i[expressedGenes,"npt"];y[is.na(y)] <- 0
    
    x+y
  })
  ,sapply(outTmp1[["cytoplasmic"]],function(i)
  {
    x <- i[expressedGenes,"cyt"];x[is.na(x)] <- 0
    x
  }))
}
rownames(totalCountsTmp) <- expressedGenes

if(mergeSamples)
{
  spikesCountsTmp <- spikesCountsTmpMerged
}

# Normalization on library size
normCountsTmp <- apply(totalCountsTmp,2,function(i)1e6*(i/sum(i)))

# Human reads (all - ERCC - ENO2)
if(!is.null(yeastCounts)) # If some external yeast counts are provided
{
  if(identical(names(yeastCounts),colnames(spikesCountsTmp))) # If external yeast counts are provided for all the samples
  {
    spikesCountsTmp <- rbind(yeastCounts,spikesCountsTmp)
    rownames(spikesCountsTmp)[[1]] <- "ENO2"
  }else if(all(names(yeastCounts)%in%colnames(spikesCountsTmp))&"ENO2"%in%rownames(spikesCountsTmp)){ # If external yeast counts are provided for some of the samples
    spikesCountsTmp["ENO2",names(yeastCounts)] <- yeastCounts
  }else if(mergeSamples)
  {
    tmp <- spikesCountsTmp["ENO2",]
    names(tmp) <- sapply(strsplit(names(tmp),"_"),"[[",1)
    names(yeastCounts) <- sapply(strsplit(names(yeastCounts),"_"),"[[",1)
    if(all(names(yeastCounts)%in%names(tmp)))
    {
      tmp[names(yeastCounts)] <- tmp[names(yeastCounts)]+yeastCounts
      spikesCountsTmp["ENO2",] <- tmp
    }else{
      print("Issues with yeast counts.")
      stop()
    }
  }else{
    print("Issues with yeast counts.")
    stop()
  }
}
human_reads = tbl[colnames(spikesCountsTmp),"reads"] - colSums(spikesCountsTmp)

# All human aligned reads
sum_RC <- colSums(totalCountsTmp)

# Normalization factors
normFactor <- (tbl[colnames(spikesCountsTmp),"PolyA"]/tbl[colnames(spikesCountsTmp),"cells"]) * 
  sum_RC/human_reads * 
  1e-3/(colSums(spikesCountsTmp[!(rownames(spikesCountsTmp)%in%"ENO2"),])/sum_RC)

# Normalized counts
normalizedCounts <- sapply(seq_along(normFactor),function(i)normFactor[[i]]*normCountsTmp[,i])
colnames(normalizedCounts) <- colnames(normCountsTmp)

# Expression levels subdivision
if(mergeSamples)
{
  chrSamples <- "Chr_Merged"
  nucSamples <- "Nuc_Merged"
  cytSamples <- "Cyt_Merged"
  if(length(polySamples)>0){polySamples <- "Poly_Merged"}
}

normalizedCountsChromatin <- lapply(chrSamples,function(i)
{
  commonGenesTmp <- intersect(rownames(normalizedCounts),rownames(chromatin[[i]]))
  outTmp <- t(sapply(commonGenesTmp,function(j)
  {
    as.matrix(chromatin[[i]][j,])*as.numeric(normalizedCounts[j,i])/as.numeric(totalCountsTmp[j,i])
  }))
  colnames(outTmp) <- colnames(chromatin[[i]])
  outTmp
})
names(normalizedCountsChromatin) <- chrSamples

normalizedCountsNucleoplasmic <- lapply(nucSamples,function(i)
{
  commonGenesTmp <- intersect(rownames(normalizedCounts),rownames(nucleoplasmic[[i]]))
  outTmp <- t(sapply(commonGenesTmp,function(j)
  {
    as.matrix(nucleoplasmic[[i]][j,])*as.numeric(normalizedCounts[j,i])/as.numeric(totalCountsTmp[j,i])
  }))
  colnames(outTmp) <- colnames(nucleoplasmic[[i]])
  outTmp
})
names(normalizedCountsNucleoplasmic) <- nucSamples

normalizedCountsCytoplasmic <- lapply(cytSamples,function(i)
{
  commonGenesTmp <- intersect(rownames(normalizedCounts),rownames(cytoplasmic[[i]]))
  outTmp <- t(sapply(commonGenesTmp,function(j)
  {
    as.matrix(cytoplasmic[[i]][j,])*as.numeric(normalizedCounts[j,i])/as.numeric(totalCountsTmp[j,i])
  }))
  colnames(outTmp) <- colnames(cytoplasmic[[i]])
  outTmp
})
names(normalizedCountsCytoplasmic) <- cytSamples

if(length(polySamples)>0)
{
  normalizedCountsPolysomal <- lapply(polySamples,function(i)
  {
    commonGenesTmp <- intersect(rownames(normalizedCounts),rownames(polysomal[[i]]))
    outTmp <- t(sapply(commonGenesTmp,function(j)
    {
      as.matrix(polysomal[[i]][j,])*as.numeric(normalizedCounts[j,i])/as.numeric(totalCountsTmp[j,i])
    }))
    colnames(outTmp) <- colnames(polysomal[[i]])
    outTmp
  })
  names(normalizedCountsPolysomal) <- polySamples    
}

# Second set of outputs
outTmp2 <- list("spikesCounts"=spikesCountsTmp
                , "sum_RC"=sum_RC
                , "human_reads"=human_reads
                , "normFactor"=normFactor
                , "normalizedCounts"=normalizedCounts
                , "normalizedCountsChromatin"=normalizedCountsChromatin
                , "normalizedCountsNucleoplasmic"=normalizedCountsNucleoplasmic
                , "normalizedCountsCytoplasmic"=normalizedCountsCytoplasmic
                , "normalizedCountsPolysomal"=normalizedCountsPolysomal)

### Split of genes according to the expression levels of different RNA species
splitFunctionTmp <- function(x,genes)
{
  lapply(seq_along(x),function(j){i <- x[[j]];i[genes[[j]],]})
}

normalizedCountsYesChpNpP <- normalizedCountsYesChpPNoNp <- normalizedCountsYesPNoChp <- NULL

if(length(polySamples)>0)
{
  normalizedCountsYesChpNpP <- list("chromatin"=splitFunctionTmp(normalizedCountsChromatin,outTmp1$expressedGenesYesChpNpP)
                                    ,"nucleoplasmic"=splitFunctionTmp(normalizedCountsNucleoplasmic,outTmp1$expressedGenesYesChpNpP)
                                    ,"cytoplasmic"=splitFunctionTmp(normalizedCountsCytoplasmic,outTmp1$expressedGenesYesChpNpP)
                                    ,"polysomal"= splitFunctionTmp(normalizedCountsPolysomal,outTmp1$expressedGenesYesChpNpP))
  
  normalizedCountsYesChpPNoNp <- list("chromatin"=splitFunctionTmp(normalizedCountsChromatin,outTmp1$expressedGenesYesChpPNoNp)
                                      ,"nucleoplasmic"=splitFunctionTmp(normalizedCountsNucleoplasmic,outTmp1$expressedGenesYesChpPNoNp)
                                      ,"cytoplasmic"=splitFunctionTmp(normalizedCountsCytoplasmic,outTmp1$expressedGenesYesChpPNoNp)
                                      ,"polysomal"= splitFunctionTmp(normalizedCountsPolysomal,outTmp1$expressedGenesYesChpPNoNp))
  
  normalizedCountsYesPNoChp <- list("chromatin"=splitFunctionTmp(normalizedCountsChromatin,outTmp1$expressedGenesYesPNoChp)
                                    ,"nucleoplasmic"=splitFunctionTmp(normalizedCountsNucleoplasmic,outTmp1$expressedGenesYesPNoChp)
                                    ,"cytoplasmic"=splitFunctionTmp(normalizedCountsCytoplasmic,outTmp1$expressedGenesYesPNoChp)
                                    ,"polysomal"= splitFunctionTmp(normalizedCountsPolysomal,outTmp1$expressedGenesYesPNoChp))
  
  ### Sum of nucleoplasmic premature and mature RNA for genes without Np and/or Chp
  normalizedCountsYesChpPNoNp$nucleoplasmic <- lapply(normalizedCountsYesChpPNoNp$nucleoplasmic,function(i)
  {
    i[,"npn0"] <- i[,"npn0"]+i[,"nmn0"]
    i[,paste0("npn",labelingTime)] <- i[,paste0("npn",labelingTime)]+i[,paste0("nmn",labelingTime)]
    i[,"npp0"] <- i[,"npp0"]+i[,"nmp0"]
    i[,paste0("npp",labelingTime)] <- i[,paste0("npp",labelingTime)]+i[,paste0("nmp",labelingTime)]
    i[,"npt"] <- i[,"npt"]+i[,"nmt"]
    i <- i[,-grep("^nm",colnames(i))]
    colnames(i) <- gsub("^np","n",colnames(i))
    i
  })
  
  normalizedCountsYesPNoChp$nucleoplasmic <- lapply(normalizedCountsYesPNoChp$nucleoplasmic,function(i)
  {
    i[,"npn0"] <- i[,"npn0"]+i[,"nmn0"]
    i[,paste0("npn",labelingTime)] <- i[,paste0("npn",labelingTime)]+i[,paste0("nmn",labelingTime)]
    i[,"npp0"] <- i[,"npp0"]+i[,"nmp0"]
    i[,paste0("npp",labelingTime)] <- i[,paste0("npp",labelingTime)]+i[,paste0("nmp",labelingTime)]
    i[,"npt"] <- i[,"npt"]+i[,"nmt"]
    i <- i[,-grep("^nm",colnames(i))]
    colnames(i) <- gsub("^np","n",colnames(i))
    i
  })
  
  ### Sum of chromatin premature and mature RNA for genes without Chp
  normalizedCountsYesPNoChp$chromatin <- lapply(normalizedCountsYesPNoChp$chromatin,function(i)
  {
    i[,"chpn0"] <- i[,"chpn0"]+i[,"chmn0"]
    i[,paste0("chpn",labelingTime)] <- i[,paste0("chpn",labelingTime)]+i[,paste0("chmn",labelingTime)]
    i[,"chpp0"] <- i[,"chpp0"]+i[,"chmp0"]
    i[,paste0("chpp",labelingTime)] <- i[,paste0("chpp",labelingTime)]+i[,paste0("chmp",labelingTime)]
    i[,"chpt"] <- i[,"chpt"]+i[,"chmt"]
    i <- i[,-grep("^chm",colnames(i))]
    colnames(i) <- gsub("^chp","ch",colnames(i))
    i
  })
}

normalizedCountsYesChpNpNoP <- list("chromatin"=splitFunctionTmp(normalizedCountsChromatin,outTmp1$expressedGenesYesChpNpNoP)
                                    ,"nucleoplasmic"=splitFunctionTmp(normalizedCountsNucleoplasmic,outTmp1$expressedGenesYesChpNpNoP)
                                    ,"cytoplasmic"=splitFunctionTmp(normalizedCountsCytoplasmic,outTmp1$expressedGenesYesChpNpNoP))

normalizedCountsYesChpNoNpP <- list("chromatin"=splitFunctionTmp(normalizedCountsChromatin,outTmp1$expressedGenesYesChpNoNpP)
                                    ,"nucleoplasmic"=splitFunctionTmp(normalizedCountsNucleoplasmic,outTmp1$expressedGenesYesChpNoNpP)
                                    ,"cytoplasmic"=splitFunctionTmp(normalizedCountsCytoplasmic,outTmp1$expressedGenesYesChpNoNpP))

normalizedCountsNoChpP <- list("chromatin"=splitFunctionTmp(normalizedCountsChromatin,outTmp1$expressedGenesNoChpP)
                               ,"nucleoplasmic"=splitFunctionTmp(normalizedCountsNucleoplasmic,outTmp1$expressedGenesNoChpP)
                               ,"cytoplasmic"=splitFunctionTmp(normalizedCountsCytoplasmic,outTmp1$expressedGenesNoChpP))

### Sum of nucleoplasmic premature and mature RNA for genes without Np and/or Chp
normalizedCountsYesChpNoNpP$nucleoplasmic <- lapply(normalizedCountsYesChpNoNpP$nucleoplasmic,function(i)
{
  i[,"npn0"] <- i[,"npn0"]+i[,"nmn0"]
  i[,paste0("npn",labelingTime)] <- i[,paste0("npn",labelingTime)]+i[,paste0("nmn",labelingTime)]
  i[,"npp0"] <- i[,"npp0"]+i[,"nmp0"]
  i[,paste0("npp",labelingTime)] <- i[,paste0("npp",labelingTime)]+i[,paste0("nmp",labelingTime)]
  i[,"npt"] <- i[,"npt"]+i[,"nmt"]
  i <- i[,-grep("^nm",colnames(i))]
  colnames(i) <- gsub("^np","n",colnames(i))
  i
})

normalizedCountsNoChpP$nucleoplasmic <- lapply(normalizedCountsNoChpP$nucleoplasmic,function(i)
{
  i[,"npn0"] <- i[,"npn0"]+i[,"nmn0"]
  i[,paste0("npn",labelingTime)] <- i[,paste0("npn",labelingTime)]+i[,paste0("nmn",labelingTime)]
  i[,"npp0"] <- i[,"npp0"]+i[,"nmp0"]
  i[,paste0("npp",labelingTime)] <- i[,paste0("npp",labelingTime)]+i[,paste0("nmp",labelingTime)]
  i[,"npt"] <- i[,"npt"]+i[,"nmt"]
  i <- i[,-grep("^nm",colnames(i))]
  colnames(i) <- gsub("^np","n",colnames(i))
  i
})

### Sum of chromatin premature and mature RNA for genes without Chp
normalizedCountsNoChpP$chromatin <- lapply(normalizedCountsNoChpP$chromatin,function(i)
{
  i[,"chpn0"] <- i[,"chpn0"]+i[,"chmn0"]
  i[,paste0("chpn",labelingTime)] <- i[,paste0("chpn",labelingTime)]+i[,paste0("chmn",labelingTime)]
  i[,"chpp0"] <- i[,"chpp0"]+i[,"chmp0"]
  i[,paste0("chpp",labelingTime)] <- i[,paste0("chpp",labelingTime)]+i[,paste0("chmp",labelingTime)]
  i[,"chpt"] <- i[,"chpt"]+i[,"chmt"]
  i <- i[,-grep("^chm",colnames(i))]
  colnames(i) <- gsub("^chp","ch",colnames(i))
  i
})

if(length(polySamples)>0)
{
  normalizedCountsYesChpNpP <- lapply(seq_along(normalizedCountsYesChpNpP[[1]]),function(j)do.call(cbind,lapply(normalizedCountsYesChpNpP,function(i)i[[j]])))
  normalizedCountsYesChpPNoNp <- lapply(seq_along(normalizedCountsYesChpPNoNp[[1]]),function(j)do.call(cbind,lapply(normalizedCountsYesChpPNoNp,function(i)i[[j]])))
  normalizedCountsYesPNoChp <- lapply(seq_along(normalizedCountsYesPNoChp[[1]]),function(j)do.call(cbind,lapply(normalizedCountsYesPNoChp,function(i)i[[j]])))
  
  normalizedCountsYesChpNpP <- lapply(normalizedCountsYesChpNpP,function(i){i[abs(i)<1e-10] <- 1e-10;i})
  normalizedCountsYesChpPNoNp <- lapply(normalizedCountsYesChpPNoNp,function(i){i[abs(i)<1e-10] <- 1e-10;i})
  normalizedCountsYesPNoChp <- lapply(normalizedCountsYesPNoChp,function(i){i[abs(i)<1e-10] <- 1e-10;i})
  
  names(normalizedCountsYesChpNpP) <- 
    names(normalizedCountsYesChpPNoNp) <- 
    names(normalizedCountsYesPNoChp) <- samples
}

normalizedCountsYesChpNpNoP <- lapply(seq_along(normalizedCountsYesChpNpNoP[[1]]),function(j)do.call(cbind,lapply(normalizedCountsYesChpNpNoP,function(i)i[[j]])))
normalizedCountsYesChpNoNpP <- lapply(seq_along(normalizedCountsYesChpNoNpP[[1]]),function(j)do.call(cbind,lapply(normalizedCountsYesChpNoNpP,function(i)i[[j]])))
normalizedCountsNoChpP <- lapply(seq_along(normalizedCountsNoChpP[[1]]),function(j)do.call(cbind,lapply(normalizedCountsNoChpP,function(i)i[[j]])))

normalizedCountsYesChpNpNoP <- lapply(normalizedCountsYesChpNpNoP,function(i){i[abs(i)<1e-10] <- 1e-10;i})
normalizedCountsYesChpNoNpP <- lapply(normalizedCountsYesChpNoNpP,function(i){i[abs(i)<1e-10] <- 1e-10;i})
normalizedCountsNoChpP <- lapply(normalizedCountsNoChpP,function(i){i[abs(i)<1e-10] <- 1e-10;i})

names(normalizedCountsYesChpNpNoP) <- 
  names(normalizedCountsYesChpNoNpP) <- 
  names(normalizedCountsNoChpP) <- samples

### At this point we assume these are replicates
## Selection of the common genes 
if(length(samples)>1)
{
  normalizedCountsList <- lapply(list(normalizedCountsYesChpNpP
                                      ,normalizedCountsYesChpPNoNp
                                      ,normalizedCountsYesPNoChp
                                      ,normalizedCountsYesChpNpNoP
                                      ,normalizedCountsYesChpNoNpP
                                      ,normalizedCountsNoChpP)
                                 ,function(l){
                                   if(!is.null(l))
                                   {
                                     cmn_genes1 <- table(unlist(sapply(l,rownames)))
                                     cmn_genes1 <- names(cmn_genes1[cmn_genes1==length(l)])
                                     
                                     normalizedCountsTmp<-lapply(seq_along(l),function(i)
                                     { 
                                       j <- l[[i]]
                                       j <- j[cmn_genes1,]
                                       if(!is.matrix(j))
                                       {
                                         k <- j
                                         j <- matrix(j,nrow=1)
                                         colnames(j) <- names(k)
                                         rownames(j) <- cmn_genes1
                                         j
                                       }else{j}
                                     })
                                     names(normalizedCountsTmp) <- samples
                                     normalizedCountsTmp
                                   }else{NULL}
                                 })
  
  names(normalizedCountsList) <- c("normalizedCountsYesChpNpP"
                                   ,"normalizedCountsYesChpPNoNp"
                                   ,"normalizedCountsYesPNoChp"
                                   ,"normalizedCountsYesChpNpNoP"
                                   ,"normalizedCountsYesChpNoNpP"
                                   ,"normalizedCountsNoChpP")
  
  ### Final normalization step to account for differences in RNA production between experiments
  scalingFunction <- function(alpha,x,med)
  {
    abs(median(alpha*log10(x) - log10(med),na.rm=TRUE))
  }
  
  if(!is.null(normalizedCountsList$normalizedCountsYesPNoChp))
  {
    normalizedCounts <- normalizedCountsList$normalizedCountsYesPNoChp
  }else{
    normalizedCounts <- normalizedCountsList$normalizedCountsNoChpP
  }
  normalizedCountsMedian <- sapply(1:ncol(normalizedCounts[[1]]),function(i){apply(sapply(normalizedCounts,function(j)j[,i]),1,median)})
  normalizedCountsScaling <- sapply(normalizedCounts,function(i)optimize(scalingFunction,x=i,med=normalizedCountsMedian,lower=1e-6,upper=1e6)$minimum)    
  
  normalizedCountsListScaled <- lapply(normalizedCountsList,function(j)
  {
    jScales <- lapply(seq_along(j),function(i)j[[i]]**as.numeric(normalizedCountsScaling[[i]]))
    names(jScales) <- names(i)
    jScales
  })
  
  normalizedCountsListMean <- lapply(normalizedCountsListScaled,function(normalizedCountsTmp)
  {
    if(length(normalizedCountsTmp)!=0)
    {
      normalizedCountsTmpMean <- matrix(mapply(function(x,y) mean(c(x,y)),normalizedCountsTmp[[1]]
                                               ,normalizedCountsTmp[[2]]), ncol=ncol(normalizedCountsTmp[[1]]))
      rownames(normalizedCountsTmpMean) <- rownames(normalizedCountsTmp[[1]])
      colnames(normalizedCountsTmpMean) <- colnames(normalizedCountsTmp[[1]])
      normalizedCountsTmpMean
    }else{
      NULL
    }
  })
  
  normalizedCountsListSd <- lapply(normalizedCountsListScaled,function(normalizedCountsTmp)
  {
    if(length(normalizedCountsTmp)!=0)
    {
      normalizedCountsTmpSd <- matrix(mapply(function(x,y) sd(c(x,y)),normalizedCountsTmp[[1]]
                                             ,normalizedCountsTmp[[2]]), ncol=ncol(normalizedCountsTmp[[1]]))
      rownames(normalizedCountsTmpSd) <- rownames(normalizedCountsTmp[[1]])
      colnames(normalizedCountsTmpSd) <- colnames(normalizedCountsTmp[[1]])
      normalizedCountsTmpSd
    }else{
      NULL
    }
  })
}else{
  normalizedCountsScaling <- NULL
  normalizedCountsList <- NULL
  normalizedCountsListScaled <- NULL
  normalizedCountsListMean <- NULL
  normalizedCountsListSd <- NULL
}

outTmp3 <- list(normalizedCountsYesChpNpP=normalizedCountsYesChpNpP
                ,normalizedCountsYesChpPNoNp=normalizedCountsYesChpPNoNp
                ,normalizedCountsYesPNoChp=normalizedCountsYesPNoChp
                ,normalizedCountsYesChpNpNoP=normalizedCountsYesChpNpNoP
                ,normalizedCountsYesChpNoNpP=normalizedCountsYesChpNoNpP
                ,normalizedCountsNoChpP=normalizedCountsNoChpP
                ,normalizedCountsScaling=normalizedCountsScaling
                ,normalizedCountsList=normalizedCountsList
                ,normalizedCountsListScaled=normalizedCountsListScaled
                ,normalizedCountsListMean=normalizedCountsListMean
                ,normalizedCountsListSd=normalizedCountsListSd)

return(c(outTmp1,outTmp2,outTmp3))
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
  
  # Sorting of expression levels.
  colOrderTmp <- sort(colnames(expressionData))
  expressionData <- expressionData[,colOrderTmp]
  expressionDataDev <- expressionDataDev[,colOrderTmp]
  
  # If standard deviation are missing assume a CV of 1.
  if(is.null(expressionDataDev)){expressionDataDev<-expressionData}
  
  # If the optimization in the logarithmic space is required transform the initial rates.
  if(logOptim){initialRates <- lapply(initialRates,log)}
  
  # Selection of the configuration to be simulated according to the provided example rates or the RNA species.
  # This means definition of the modelSimulation function (ode implementation), of dataGeneration function (ode solution),
  # and in case of polysomal RNA profiling exclusion of nascent RNA from the species for the cost function optimization.
  
  if("k9"%in%names(unlist(initialRates))
     &any(grepl("^chpp",colnames(expressionData)))
     &any(grepl("^npp",colnames(expressionData)))
     &any(grepl("^py",colnames(expressionData))))
  {
    dataGeneration <- dataGenerationFullND
    # excludeSpecies=colnames(expressionData[,c(grep("^pyn",colnames(expressionData))
    #                                          ,grep("^pyp",colnames(expressionData)))])
    print("Inference mode: Nuclear decay.")
  }else if(any(grepl("^chpp",colnames(expressionData)))
           &any(grepl("^npp",colnames(expressionData)))
           &!any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: No polysomal RNA, Chromatin premature RNA, Nucleoplasmic premature RNA.")
    dataGeneration <- dataGenerationNoPoly
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k7","k9","k10")]))
  }else if(any(grepl("^chpp",colnames(expressionData)))
           &any(grepl("^npp",colnames(expressionData)))
           &any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: Full model")
    dataGeneration <- dataGenerationFullCP
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k9")]))
  }else if(any(grepl("^chpp",colnames(expressionData)))
           &!any(grepl("^npp",colnames(expressionData)))
           &!any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: No polysomal RNA, Chromatin premature RNA.")
    dataGeneration <- dataGenerationNoPolyNoNucP
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k4","k5","k7","k9","k10")]))
  }else if(any(grepl("^chpp",colnames(expressionData)))
           &!any(grepl("^npp",colnames(expressionData)))
           &any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: Yes polysomal RNA, Chromatin premature RNA.")
    dataGeneration <- dataGenerationNoNucP
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k4","k5","k9")]))
  }else if(!any(grepl("^chpp",colnames(expressionData)))
           &!any(grepl("^npp",colnames(expressionData)))
           &!any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: No polysomal RNA, No premature RNA.")
    dataGeneration <- dataGenerationNoPolyNoP
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k2","k4","k5","k7","k9","k10")]))
  }else if(!any(grepl("^chpp",colnames(expressionData)))
           &!any(grepl("^npp",colnames(expressionData)))
           &any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: Yes polysomal RNA, No premature RNA.")
    dataGeneration <- dataGenerationNoP
    # excludeSpecies=colnames(expressionData[,c(grep("^pyn",colnames(expressionData))
    #                                          ,grep("^pyp",colnames(expressionData)))])
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k2","k4","k5","k9")]))
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
                       ,control=list(maxit=1e9,reltol=1e-50)) # Optimization algorithm FIXED parameter.
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
  rownames(inferedRates) <- rownames(expressionData)
  
  # Generation of the expression levels for the optimal set of rates.
  inferedData <- t(sapply(1:nrow(inferedRates),function(i)dataGeneration(exampleRates=inferedRates[i,] # Rates.
                                                                         ,TauFractions=TauFractions # Time points with cellular fractionation.
                                                                         ,TauPoly=TauPoly # Time points with polysomal profiling.
                                                                         ,TauTotal=TauTotal))) # Time points with total nascent and pre-existing RNA profiling.
  rownames(inferedData) <- rownames(inferedRates)
  
  # Loop to estimate for each gene the model Log Likelihood as a specific cost function.
  LogL <- t(sapply(1:nrow(expressionData),function(x){costFunction(par=inferedRates[x,] # Set of rates.
                                                                   ,data=expressionData[x,] # Experimental data.
                                                                   ,dev=expressionDataDev[x,] # Experimental data standard deviation.
                                                                   ,TauFractions = TauFractions # Time points with cellular fractionation.
                                                                   ,TauPoly = TauPoly # Time points with polysomal profiling.
                                                                   ,TauTotal = TauTotal # Time points with total RNA profiling.
                                                                   ,dataGeneration=dataGeneration # Data generation function.
                                                                   ,parFixed=NULL # List of rates to be excluded from the optimization.
                                                                   ,logOptim=FALSE # Always FALSE because we previously converted (potential) logarithmic rates.
                                                                   ,excludeSpecies=NULL # List of species to be excluded from the cost function.
                                                                   ,lowB=1e-6 # Lower bound for rates.
                                                                   ,upB=1e4 # Upper bound for rates.
                                                                   ,FlagDev="LL" # We are interested in LogLikelihood.
                                                                   ,lambda=0.05) # Regularization strength.
  }))
  
  # AIC and BIC goodness of fit metrics.
  AIC <- 2*ncol(inferedRates) - 2*LogL
  BIC <- ncol(inferedRates)*log(ncol(expressionData))- 2*LogL
  
  # All goodness of fit metrics.
  metrics <- rbind(LL=LogL,AIC=AIC,BIC=BIC)
  rownames(metrics) <- c("LL","AIC","BIC")
  colnames(metrics) <- rownames(expressionData)
  
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

#Full model with cytoplasmic premature
fullModelSimulationCP <- function(t,y,parms)
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
  dy3=k4*y[1] - (k5+k6)*y[3]                     #y3=Np
  dy4=k3*y[2] + k5*y[3] - k6*y[4]           #y4=Nm
  dy5=k6*(y[3]+y[4]) - (k7+k8)*y[5]                #y5=C
  dy6=k7*y[5] - k10*y[6]                     #y6=P             
  
  nascentTmp <- c(dy1,dy2,dy3,dy4,dy5,dy6)
  
  # System of equations for PRE-EXISTING RNA
  dy7=0 - (k2+k4)*y[7]                       #y7=Chp
  dy8=k2*y[7] - k3*y[8]                      #y8=Chm
  dy9=k4*y[7] - (k5+k6)*y[9]                      #y9=Np
  dy10=k3*y[8] + k5*y[9] - k6*y[10]          #y10=Nm
  dy11=k6*(y[9]+y[10]) - (k7+k8)*y[11]              #y11=C
  dy12=k7*y[11] - k10*y[12]                   
  
  preExistingTmp <- c(dy7,dy8,dy9,dy10,dy11,dy12)
  
  list(c(nascentTmp,preExistingTmp))
}

dataGenerationFullCP <- function(exampleRates # Set of rate.
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
            ,y9= k1*k4/((k5+k6)*(k2+k4))
            ,y10= ((k1/(k2+k4))*(k2+(k4*k5)/(k5+k6)))/k6
            ,y11= k1/(k7+k8)
            ,y12= k7*k1/(k10*(k7+k8)))
  
  # All time-points to be simulated
  TauGlobals <- unique(c(TauFractions,TauPoly,TauTotal))
  
  # Numerical solution of the ODE model.
  exampleData <- ode(times=TauGlobals # Time points.
                     ,y=yini # Initial conditions.
                     ,func=fullModelSimulationCP # Model.
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