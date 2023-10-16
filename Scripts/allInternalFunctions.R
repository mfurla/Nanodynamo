### Loop functions
mcsapply <- function( X, FUN, ... ) do.call('cbind', mclapply( X, FUN, ... ))

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

  nucleoplasmicExpressedGenes <- lapply(nucleoplasmic,rownames)

  # Removal of cytoplasmic premature reads counts
  cytoplasmicNoPremature <- lapply(names(cytoplasmicP),function(i)
  {
    cytoplasmic[[i]][[paste0("cyp",labelingTime)]] <- cytoplasmic[[i]][[paste0("cyp",labelingTime)]] - cytoplasmicP[[i]][[paste0("cypp",labelingTime)]]
    cytoplasmic[[i]][[paste0("cyn",labelingTime)]] <- cytoplasmic[[i]][[paste0("cyn",labelingTime)]] - cytoplasmicP[[i]][[paste0("cypn",labelingTime)]]
    cytoplasmic[[i]]$cyn0 <- cytoplasmic[[i]]$cyn0 - cytoplasmicP[[i]]$cypn0
    cytoplasmic[[i]]$cyp0 <- cytoplasmic[[i]]$cyp0 - cytoplasmicP[[i]]$cypp0
    cytoplasmic[[i]]$cyt <- cytoplasmic[[i]]$cyt - cytoplasmicP[[i]]$cypt
    cytoplasmic[[i]]
  })
  names(cytoplasmicNoPremature) <- names(cytoplasmicP)
 
  cytoplasmicExpressedGenes <- lapply(cytoplasmicNoPremature,rownames)
  
  if(length(polySamples)>0)
  {
    # Removal of polysomal premature reads counts
    polysomalNoPremature <- lapply(names(polysomalP),function(i)
    {
      polysomal[[i]][[paste0("pyp",labelingTimePoly)]] <- polysomal[[i]][[paste0("pyp",labelingTimePoly)]] - polysomalP[[i]][[paste0("pypp",labelingTimePoly)]]
      polysomal[[i]][[paste0("pyn",labelingTimePoly)]] <- polysomal[[i]][[paste0("pyn",labelingTimePoly)]] - polysomalP[[i]][[paste0("pypn",labelingTimePoly)]]
      polysomal[[i]]$pyn0 <- polysomal[[i]]$pyn0 - polysomalP[[i]]$pypn0
      polysomal[[i]]$pyp0 <- polysomal[[i]]$pyp0 - polysomalP[[i]]$pypp0
      polysomal[[i]]$pyt <- polysomal[[i]]$pyt - polysomalP[[i]]$pypt
      polysomal[[i]]
    })
    names(polysomalNoPremature) <- names(polysomalP)
    
    polysomalExpressedGenes <- lapply(polysomalNoPremature,rownames)
    
  }

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
        data <- data[names(simulatedData)]
        dev <- dev[names(simulatedData)]
      }

      # Exclusion of configurations with significantly negative expression levels.
      simulatedData[abs(simulatedData)<1e-6] <- 1e-6
      if(any(simulatedData<0))
      {
        return(1e8)
      }else{
        # Selection of the cost function.
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

#Function for the inference of rates
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
    print("Inference mode: Nuclear decay. This is designed for simulated data only.")
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
    dataGeneration <- dataGenerationFull
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

## Full model with nuclear decay
fullModelSimulationND <- function(t,y,parms)
{
  k1=parms[["k1"]] # Synthesis.
  k2=parms[["k2"]] # Co-transcriptional splicing.
  k3=parms[["k3"]] # Detachment of spliced chromatin associated RNA.
  k4=parms[["k4"]] # Detachment of unspliced chromatin associated RNA.
  k5=parms[["k5"]] # Post-transcriptional splicing.
  k6=parms[["k6"]] # Export.
  k7=parms[["k7"]] # Translation.
  k8=parms[["k8"]] # Cytoplasmic decay
  k9=parms[["k9"]] # Nuclear decay.
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

dataGenerationFullND <- function(exampleRates # Set of rate.
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
  k9=exampleRates[["k9"]] # Nuclear decay.
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
                    ,func=fullModelSimulationND # Model.
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

### Nuclear decay on nucleoplasmic premature
fullModelSimulationND_prem <- function(t,y,parms)
{
  k1=parms[["k1"]] # Synthesis.
  k2=parms[["k2"]] # Co-transcriptional splicing.
  k3=parms[["k3"]] # Detachment of spliced chromatin associated RNA.
  k4=parms[["k4"]] # Detachment of unspliced chromatin associated RNA.
  k5=parms[["k5"]] # Post-transcriptional splicing.
  k6=parms[["k6"]] # Export.
  k7=parms[["k7"]] # Translation.
  k8=parms[["k8"]] # Cytoplasmic decay
  k11=parms[["k11"]] # Nuclear decay.
  k10=parms[["k10"]] # Polysomal decay

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

dataGenerationFullND_prem <- function(exampleRates # Set of rate.
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
  k11=exampleRates[["k11"]] # Nuclear decay.
  k10=exampleRates[["k10"]] # Polysomal decay

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
                    ,func=fullModelSimulationND_prem # Model.
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

### Data simulations
##  - Function to mimic replicates,
##  - Data simulation function.
dataAverageFunction <- function(v
                               ,nRep
                               ,CV)
{
  generated_gaussian <- sapply(v, function(x){rnorm(nRep,m=x,sd=abs(CV*x))})
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
    deviazioni <- generated_gaussian[1,]
  }
  list(medie,deviazioni)
}

# Function to generate simulated data
simulateData <- function(exampleRates # Mean rates to sample.
                        ,TauFractions # Time points with cellular fractionation.
                        ,TauPoly # Time points with polysomal profiling.
                        ,TauTotal # Time points with total nascent and pre-existing RNA profiling.
                        ,noise # TRUE to add noise to data.
                        ,CV # Variation Coefficient for noise.
                        ,Reps # Number of replicates.
                        ,nGenes # Number of genes to be simulated.
                        ,seed # Seed for reproducibility.
                        ,ZeroThresh # Minimum expression value.
                        ,MultFact) # Number of nGenes to be modeled.
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
  { #Full Model with also nuclear decay
    modelSimulation <- fullModelSimulationND
    dataGeneration <- dataGenerationFullND
    print("Full Model with nuclear decay")
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
  elements <- c("Noise","Reps","CV","TauFractions","TauPoly","TauTotal",names(exampleRates))
  values <- c(noise,Reps,CV,paste0(TauFractions,collapse="_"),paste0(TauPoly,collapse = "_"),paste0(TauTotal,collapse = "_"),unname(exampleRates))
  nameTmp <- character(0)
  for (i in 1:length(elements)){
    name <- paste0(c(elements[i],values[i]),collapse = "_")
    nameTmp <- paste0(c(nameTmp,name),collapse ="_")
  }

  # Sampling of the rates.
  set.seed(seed)
  multipleExampleRates <- sapply(exampleRates,function(x)rnorm(nGenes*MultFact,mean=x,sd=5*x))

  #Selection of positive rate-sets.
  multipleExampleRates <- multipleExampleRates[apply(multipleExampleRates,1,function(row)all(row>=0)),]
  
  # Inclusion of the original rates set.
  multipleExampleRates <- rbind(exampleRates,multipleExampleRates)
  rownames(multipleExampleRates) <- 1:nrow(multipleExampleRates)

  # Data generation and addition of noise at the expression level.
  set.seed(seed)
  Obj <- t(apply(multipleExampleRates,1,function(par)
  {
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
                                        ,CV=CV) # Variation Coefficient.
      return(ListDataTmp)
    }else{
      return(exampleDataTmp)
    }
  }))
  
  # Extraction of the quantities of interest.
  exampleData <- t(sapply(Obj,"[[",1))          #means 
  rownames(exampleData) <- 1:nrow(exampleData)
  DevDataTmp <- t(sapply(Obj,"[[",2))           #deviations
  rownames(DevDataTmp) <- 1:nrow(DevDataTmp)
  
  # Conditions with low coverage are set to the minimum coverage value.
  exampleData[abs(exampleData)<ZeroThresh] <- ZeroThresh

  # Selection of suitable examples without negative expression levels (not biologically reasonable).
  exampleData[which(exampleData<0)] <- NaN
  exampleData <- exampleData[apply(exampleData,1,function(r)all(is.finite(r))),]

  # Selection of at last nGenes examples.
  if(nrow(exampleData)<nGenes){print(paste0("Warning: less than ",nGenes," genes provided."))}
  exampleData <- exampleData[1:min(nrow(exampleData),nGenes),]
  multipleExampleRates <- multipleExampleRates[rownames(exampleData),]
  DevDataTmp <- DevDataTmp[rownames(exampleData),]

  rownames(exampleData) <- 1:nrow(exampleData)
  rownames(DevDataTmp) <- 1:nrow(DevDataTmp)
  rownames(multipleExampleRates) <- 1:nrow(multipleExampleRates)

  # Output.
  return(list("exampleData"=exampleData,"DevDataTmp" = DevDataTmp
              ,"exampleRates"=multipleExampleRates,"name"=nameTmp
              ,"TauGlobals"=TauGlobals,"TauFractions"=TauFractions
              ,"TauPoly"=TauPoly,"TauTotal"=TauTotal
              ,"modelSimulation"=modelSimulation,"dataGeneration"=dataGeneration))
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
  AllRates <- colnames(inferedRates)

  #names for the titles of the plots
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

  RateNames <- RateNames[c("k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11")%in%colnames(inferedRates)] #select actually present rates
  myPalette <- myPalette[c("k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11")%in%colnames(inferedRates)] #set a color for every rate
  plot_list <- list()

  names(AllRates) <- AllRates
  AllRates <- AllRates[ratesOrder[ratesOrder%in%AllRates]]
  
  inferedRates <- inferedRates[,AllRates] #select only useful columns

  #generate plot objects
  for (i in 1:length(AllRates)){
    if(length(grep(AllRates[i], MissingRates))==0){
      inferedRates_curr <- log10(data.frame(rate = inferedRates[, AllRates[i]]))
      rownames(inferedRates_curr) <- c()
      
      plot_list[[i]] <- ggplot(inferedRates_curr, aes(x = rate)) +
        geom_histogram(aes(y = ..density..), fill=myPalette[i]) +
        geom_density(alpha=.4, fill="grey") +
        theme_classic() +
        labs(x= "", y="density") +
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

  #actually present rates
  ActualRates <- names(which(table(unlist(sapply(inferedRates,colnames)))==length(objects)))

  #names for the titles of the plots
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

  if(!allRates) # if some rate misses, then exclude it
  {
    inferedRates <- lapply(inferedRates,function(i)i[,ActualRates])
    RateNames <- RateNames[names(RateNames)%in%ActualRates]
  }
  #creating palette
  myPalette <- RColorBrewer::brewer.pal(10,"Set3")
  names(myPalette) <- c("k1","k2","k3","k4","k5","k6","k7","k8","k9","k10")
  myPalette <- myPalette[names(RateNames)]

  sortedRates <- names(RateNames)
  names(sortedRates) <- sortedRates

  #order rates
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
  ratesTmp1 <- object1$inferedRates #ref object
  ratesTmp2 <- object2$inferedRates #object to compare with ref

  #pick the right columns order
  if(ncol(ratesTmp1)==9){
    ratesTmp1 <- ratesTmp1[,c(1:6,8,7,9)]
    colnames(ratesTmp1) <- paste0("k",1:9)
  
    ratesTmp2 <- ratesTmp2[,c(1:6,8,7,9)]
    colnames(ratesTmp2) <- paste0("k",1:9)
  }else{
    colnames(ratesTmp1) <- paste0("k",1:7)
    colnames(ratesTmp2) <- paste0("k",1:7)
  }

  #build dataframes
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
                        ,ratesOrder=c("k1","k2","k3","k4","k5","k11","k9","k6","k8","k7","k10"))
{
  ExampleData <- object$expressionData
  polyFlag <- ("pyt"%in%colnames(ExampleData))
  
  if(polyFlag) #there is polysomal RNA
  {
    ExampleData <- ExampleData[,c("chpt","chmt","npt","nmt","cyt","pyt")]
  }else{
    #no polysomal RNA
    ExampleData <- ExampleData[,c("chpt","chmt","npt","nmt","cyt")]
  }

  #set colnames to ExampleData
  colnames(ExampleData) <- c("Chp","Chm","Np","Nm","C","P")[1:ncol(ExampleData)]

  inferedRates <- object$inferedRates

  sortedRates <- colnames(inferedRates)
  names(sortedRates) <- sortedRates

  #pick the right order
  sortedRates <- sortedRates[ratesOrder[ratesOrder%in%sortedRates]]

  inferedRates <- inferedRates[,sortedRates]

  RateNames <- c("k1","k2","k3","k4","k5","k6","k8","k7","k10","k9")
  names(RateNames) <- paste0("k",1:10)

  colnames(inferedRates) <- RateNames[colnames(inferedRates)]

  #Heatmap of rates
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

  #Heatmap of species and rates
  pheatmap(cbind(ExampleDataBinned,inferedRatesBinned)[ord,]
  ,cluster_rows = FALSE,cluster_cols = FALSE,annotation_row=rowAnnotation
  ,show_rownames = FALSE,color = colorRampPalette(c("white","pink","red"))(101)
  ,filename = paste0("expressionAndRatesHeatmap_",obj_name,".pdf"),width=5,height=4)
  
  saveRDS(ord2,file=paste0("Clustering_",obj_name,".rds"))
  return(ord2)
}

DifferentialHeatmapsPlot <- function(refObject
                                    ,object
                                    ,name
                                    ,n_clust
                                    ,ylim=c(-5,5)
                                    ,ratesWeight=1
                                    ,ratesOrder=c("k1","k2","k3","k4","k5","k11","k9","k6","k8","k7","k10")
                                    ,show_rownames=FALSE)
{
  ExampleData1 <- refObject$expressionData #ref object
  polyFlag1 <- ("pyt"%in%colnames(ExampleData1))
  if(polyFlag1) #if there is polysomal RNA
  {
    ExampleData1 <- ExampleData1[,c("chpt","chmt","npt","nmt","cyt","pyt")]
  }else{
    #if polysomal RNA is missing
    ExampleData1 <- ExampleData1[,c("chpt","chmt","npt","nmt","cyt")]
  }

  inferedRates1 <- refObject$inferedRates

  ExampleData2 <- object$expressionData #object to be compared with ref
  polyFlag2 <- ("pyt"%in%colnames(ExampleData2))
  if(polyFlag2)# if there is polysomal RNA
  {
    ExampleData2 <- ExampleData2[,c("chpt","chmt","npt","nmt","cyt","pyt")]
  }else{
    #if polysomal RNA is missing
    ExampleData2 <- ExampleData2[,c("chpt","chmt","npt","nmt","cyt")]
  }

  inferedRates2 <- object$inferedRates

  #Select common genes, rates, species
  commonGenes <- intersect(rownames(ExampleData1),rownames(ExampleData2))
  commonRates <- intersect(colnames(inferedRates1),colnames(inferedRates2))
  commonSpecies <- intersect(colnames(ExampleData1),colnames(ExampleData2))

  names(commonRates) <- commonRates

  commonRates <- commonRates[ratesOrder[ratesOrder%in%commonRates]]

  print(paste0(length(commonGenes)," common genes."))

  #compute correlations
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

  RateNames <- c("k1","k2","k3","k4","k5","k6","k8","k7","k10","k9")
  names(RateNames) <- paste0("k",1:10)

  colnames(corB) <- rownames(corB) <- RateNames[colnames(corB)]

  #Heatmap of differential correlations
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

  #Extract the clustering
  foe <- pheatmap(cbind(ExampleData,do.call("cbind",lapply(1:ratesWeight,function(i)inferedRates)))
  ,cluster_rows = TRUE,cluster_cols = FALSE,
  ,show_rownames = FALSE,color = colorRampPalette(c("blue","white","red"))(101)
  ,silent=TRUE)

  ord <- foe$tree_row$order
  ord2 <- cutree(foe$tree_row,n_clust)[ord]

  ord2 <- split(ord2,ord2)[unique(ord2)]

  ord2 <- unlist(sapply(seq_along(ord2),function(i){rep(LETTERS[i],length(ord2[[i]]))}))
  names(ord2) <- rownames(ExampleData)[ord]

  rowAnnotation <- as.data.frame(list("Cluster"=ord2))

  #Heatmap of data and rates
  pheatmap(cbind(ExampleData,inferedRates)[ord,]
  ,cluster_rows = FALSE,cluster_cols = FALSE,annotation_row=rowAnnotation
  ,show_rownames = show_rownames,color = colorRampPalette(c("blue","white","red"))(101)
  ,filename = paste0("differentialExpressionAndRatesHeatmap_",name,".pdf"),width=5,height=4)
  
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

    #Genes characterization of genes with fastest and slowest dynamics

    ratesNames <- c(rep(unlist(lapply(colnames(sortedRates),function(i){rep(i,top)})),2),
                        unlist(lapply(colnames(sortedRates),function(i){rep(i,nrow(sortedRates))})))
    df <- data.frame(cbind("Rank"=c(rep("TOP",top*ncol(sortedRates)),rep("BOT",bot*ncol(sortedRates)),rep("ALL",nrow(sortedRates)*9)),
                           "Rate"=ratesNames,
                           "Genes"=c(matrix(sortedRates[1:top,],ncol=1),
                                     matrix(sortedRates[(nrow(sortedRates)-bot+1):nrow(sortedRates),],ncol=1),
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

      MatCol[names(Wtest[[i]]),i] <- c("blue","red")[(Wtest[[i]]<1e-5)+1]
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

  # list of bams
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

        CGTmp <- sapply(fastaTmp,function(i)sum(strsplit(i,"")[[1]]=="C"|strsplit(i,"")[[1]]=="G")/length(strsplit(i,"")[[1]]))
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
GoodnessOfFit <- function(inferedRates,width=7,height=12,name="",lowSat=0,upSat=1)
{
  expressionDataTmp <- inferedRates$expressionData
  inferedDataTmp <- inferedRates$inferedData

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

    if(any(!is.finite(x)))
    {
      print("Not finite data-points.")
      y <- y[is.finite(x)]
      x <- x[is.finite(x)]
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
      smoothScatter(x,y
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

PlotCorrelations <- function(object,width,height,obj_name)
{
        tmp <- object$relativeErrors
        nc <- ncol(tmp)/3
        AllRates <- c("k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11")
        rates <- colnames(object$inferedRates)
        MissingRates <- setdiff(AllRates,rates)
        Present <- setdiff(c(1:11),sapply(seq_along(MissingRates),function(i){grep(MissingRates[i],AllRates)}))
  

        pdf(paste0(obj_name,".pdf"),width=width,height=height)
        par(mfrow=c(1,nc),cex=1.1)
        mains = c("Synthesis\n[fg/(M. cells*h)]"
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
        mains = mains[Present]
        for(j in 1:nc)
        {
            tryCatch(smoothScatter(log10(tmp[,j]),log10(tmp[,(nc+j)]),xlab="Inferred",ylab="Expected",main=mains[[j]],mgp=c(2,1,0),cex.lab = 1.2),error=function(e)plot(1:10,1:10,xlab="ERROR",ylab="ERROR",main=mains[[j]]))
            abline(0,1,col=2,lwd=2)
            tryCatch(points(lowess(log10(tmp[,j]),log10(tmp[,(nc+j)])),col=3,lwd=2,type="l"),error=function(e){return(NULL)})
            legend("bottomright",legend=paste0("Cor.=",round(cor(tmp[,j],tmp[,(nc+j)],method="s",use="c"),2)),bty = "n",text.font=2)
        }
  dev.off()
}

replicatesDataCorrelations <- function(object1,object2,width=7,height=12,name="",lowSat=0,upSat=1)
{
  expressionLevels1 <- object1$expressionData
  expressionLevels2 <- object2$expressionData

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

  pdf(paste0("rawDataCorrelations_",name,".pdf"),width=7,height=12)
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

    outList <- append(outList,corTmp)

    if(is.finite(corTmp))
    {
      smoothScatter(x,y
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

  pdf(paste0("ratesCorrelations_",name,".pdf"),width=7,height=12)
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

    xSat <- xSat[xSat<xSatBounds[[1]]|xSat>xSatBounds[[2]]|ySat<ySatBounds[[1]]|ySat>ySatBounds[[2]]]
    ySat <- ySat[names(xSat)]

    xSat[xSat<xSatBounds[[1]]] <- xSatBounds[[1]]
    xSat[xSat>xSatBounds[[2]]] <- xSatBounds[[2]]

    ySat[ySat<ySatBounds[[1]]] <- ySatBounds[[1]]
    ySat[ySat>ySatBounds[[2]]] <- ySatBounds[[2]]

    commonGenes <- intersect(rownames(inferredRates1),rownames(inferredRates2))

    outList <- append(outList,corTmp)

    if(is.finite(corTmp))
    {
      smoothScatter(x,y
               ,xlab="Replicate 1",ylab="Replicate 2",main=RateNames[[i]],xlim=xSatBounds,ylim=ySatBounds)
      abline(0,1,col=2,lwd=2)
      points(xSat,ySat,pch=16,col=2)
      tryCatch(points(lowess(x,y),col=3,lwd=2,type="l"),error=function(e){return(NULL)})
      legend("bottomright",legend=paste0("Cor.=",round(corTmp,2)),bty = "n",text.font=2)
    }
  }
  dev.off()
  names(outList) <- colnames(inferredRates1)[!grepl("0$",colnames(inferredRates1))]
  outList
}

simulatedRatesCorrelations <- function(object1,simObject,width=7,height=12,name="",lowSat=0,upSat=1,ratesOrder=c("k1","k2","k3","k4","k5","k11","k9","k6","k8","k7","k10"))
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

  pdf(paste0("simRatesCorrelations_",name,".pdf"),width=7,height=12)
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

    xSat <- xSat[xSat<xSatBounds[[1]]|xSat>xSatBounds[[2]]|ySat<ySatBounds[[1]]|ySat>ySatBounds[[2]]]
    ySat <- ySat[names(xSat)]

    xSat[xSat<xSatBounds[[1]]] <- xSatBounds[[1]]
    xSat[xSat>xSatBounds[[2]]] <- xSatBounds[[2]]

    ySat[ySat<ySatBounds[[1]]] <- ySatBounds[[1]]
    ySat[ySat>ySatBounds[[2]]] <- ySatBounds[[2]]

    commonGenes <- intersect(rownames(inferredRates1),rownames(inferredRates2))

    outList <- append(outList,corTmp)

    if(is.finite(corTmp))
    {
      smoothScatter(x,y
               ,xlab="Inferred",ylab="Expected",main=RateNames[[i]],xlim=xSatBounds,ylim=ySatBounds)
      abline(0,1,col=2,lwd=2)
      points(xSat,ySat,pch=16,col=2)
      tryCatch(points(lowess(x,y),col=3,lwd=2,type="l"),error=function(e){return(NULL)})
      legend("bottomright",legend=paste0("Cor.=",round(corTmp,2)),bty = "n",text.font=2)
    }
  }
  dev.off()
  names(outList) <- colnames(inferredRates1)[!grepl("0$",colnames(inferredRates1))]
  outList
}