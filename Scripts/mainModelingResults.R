### Libraries
source("/path/to/allInternalFunctions.R")

### All BAMs
## Untreated
bamPaths_Untreated <- c(Chr_Untreated1="/path/to/bam.RData"
                       ,Nuc_Untreated1="/path/to/bam.RData"
                       ,Cyt_Untreated1="/path/to/bam.RData"
                       ,Poly_Untreated1="path/to/bam.RData"
                       ,Chr_Untreated2="/path/to/bam.RData"
                       ,Nuc_Untreated2="/path/to/bam.RData"
                       ,Cyt_Untreated2="/path/to/bam.RData"
                       ,Poly_Untreated2="/path/to/bam.RData"
                       )

## PlaB
bamPaths_PlaB <- c(Chr_PlaB1="/path/to/bam.RData"
                  ,Nuc_PlaB1="/path/to/bam.RData"
                  ,Cyt_PlaB1="/path/to/bam.RData"
                  ,Poly_PlaB1="/path/to/bam.RData"
                  ,Chr_PlaB2="/path/to/bam.RData"
                  ,Nuc_PlaB2="/path/to/bam.RData"
                  ,Cyt_PlaB2="/path/to/bam.RData"
                  ,Poly_PlaB2="/path/to/bam.RData"
                  )

## HAR
bamPaths_HAR <- c(Chr_HAR1="/path/to/bam.RData"
                 ,Nuc_HAR1="/path/to/bam.RData"
                 ,Cyt_HAR1="/path/to/bam.RData"
                 ,Chr_HAR2="/path/to/bam.RData"
                 ,Nuc_HAR2="/path/to/bam.RData"
                 ,Cyt_HAR2="/path/to/bam.RData"
                 )

## LEP
bamPaths_LEP <- c(Chr_LEP1="/path/to/bam.RData"
                 ,Nuc_LEP1="/path/to/bam.RData"
                 ,Cyt_LEP1="/path/to/bam.RData"
                 ,Poly_LEP1="/path/to/bam.RData"
                 ,Chr_LEP2="/path/to/bam.RData"
                 ,Nuc_LEP2="/path/to/bam.RData"
                 ,Cyt_LEP2="/path/to/bam.RData"
                 ,Poly_LEP2="/path/to/bam.RData"
                )

### All Tail Files
## Untreated
tailsPaths_Untreated <- c(Chr_Untreated1="/path/to/polya_results.tsv"
                         ,Nuc_Untreated1="/path/to/polya_results.tsv"
                         ,Cyt_Untreated1="/path/to/polya_results.tsv"
                         ,Poly_Untreated1="/path/to/polya_results.tsv"
                         ,Poly_Untreated1="/path/to/polya_results.tsv"
                         ,Chr_Untreated2="/path/to/polya_results.tsv"
                         ,Nuc_Untreated2="/path/to/polya_results.tsv"
                         ,Cyt_Untreated2="/path/to/polya_results.tsv"
                         ,Poly_Untreated2="/path/to/polya_results.tsv"
                         )

## PlaB
tailsPaths_PlaB <- c(Chr_PlaB1="/path/to/polya_results.tsv"
                    ,Nuc_PlaB1="/path/to/polya_results.tsv"
                    ,Cyt_PlaB1="/path/to/polya_results.tsv"
                    ,Poly_PlaB1="/path/to/polya_results.tsv"
                    ,Chr_PlaB2="/path/to/polya_results.tsv"
                    ,Chr_PlaB2="/path/to/polya_results.tsv"
                    ,Nuc_PlaB2="/path/to/polya_results.tsv"
                    ,Cyt_PlaB2="/path/to/polya_results.tsv"
                    ,Poly_PlaB2="/path/to/polya_results.tsv"
                    )

## HAR
tailsPaths_HAR <- c(Chr_HAR1="/path/to/polya_results.tsv"
                   ,Nuc_HAR1="/path/to/polya_results.tsv"
                   ,Cyt_HAR1="/path/to/polya_results.tsv"
                   ,Chr_HAR2="/path/to/polya_results.tsv"
                   ,Nuc_HAR2="/path/to/polya_results.tsv"
                   ,Cyt_HAR2="/path/to/polya_results.tsv"
                   )

## LEP
tailsPaths_LEP <- c(Chr_LEP1="/path/to/polya_results.tsv"
                   ,Nuc_LEP1="/path/to/polya_results.tsv"
                   ,Cyt_LEP1="/path/to/polya_results.tsv"
                   ,Poly_LEP1="/path/to/polya_results.tsv"
                   ,Chr_LEP2="/path/to/polya_results.tsv"
                   ,Chr_LEP2="/path/to/polya_results.tsv"
                   ,Nuc_LEP2="/path/to/polya_results.tsv"
                   ,Cyt_LEP2="/path/to/polya_results.tsv"
                   ,Poly_LEP2="/path/to/polya_results.tsv"
                   )
### All expression levels
## Untreated
expressionLevelsUntreated <- readRDS("/path/to/Results/Untreated/expressionLevels.rds")
expressionLevelsUntreated_NoPoly <- readRDS("/path/to/Results/Untreated/expressionLevels_noPolysomal.rds")

## PladienolideB
expressionLevelsPlaB <- readRDS("/path/to/Results/PladienolideB/expressionLevels.rds")

## LeptomycinB
expressionLevelsLepB <- readRDS("/path/to/Results/LeptomycinB/expressionLevels.rds")

## Harringtonine
expressionLevelsHAR <- readRDS("/path/to/Results/Harringtonine/expressionLevels.rds")

### All models
## Untreated
inferedRatesUntreated1_yesChpNpP_multi <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreated1_YesChpNpP_multi.rds")
inferedRatesUntreated2_yesChpNpP_multi <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreated2_YesChpNpP_multi.rds")
inferedRatesUntreatedMerged_yesChpNpP_single <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreatedMerged_YesChpNpP_single.rds")
inferedRatesUntreatedMerged_yesChpNpP_multi <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreatedMerged_YesChpNpP_multi.rds")

inferedRatesUntreatedMerged_noPoly <- readRDS("/path/to/Results/Untreated/fourthRun_NoPolysomal/inferedRatesUntreatedMerged_YesChpNpNoP_multi.rds")

inferedRatesUntreatedMerged_nucleoplasmicPrematureExport <- readRDS("/path/to/Results/Untreated/fifthRun_NucleoplasmicPrematureExport/inferedRatesUntreatedMerged_YesChpNpP_multi.rds")

## PlaB
inferedRatesPlaB1_yesChpNpP_multi <- readRDS("/path/to/Results/PladienolideB/FirstRun_FullModel/inferedRatesPlaB1_YesChpNpP_multi.rds")
inferedRatesPlaB2_yesChpNpP_multi <- readRDS("/path/to/Results/PladienolideB/FirstRun_FullModel/inferedRatesPlaB2_YesChpNpP_multi.rds")
inferedRatesPlaBMerged_yesChpNpP_multi <- readRDS("/path/to/Results/PladienolideB/FirstRun_FullModel/inferedRatesPlaBMerged_YesChpNpP_multi.rds")

inferedRatesPlaBMerged_nucleoplasmicPrematureExport <- readRDS("/path/to/Results/PladienolideB/SecondRun_NucleoplasmicPrematureExport/inferedRatesPlaBMerged_YesChpNpP_multi.rds")

## HAR
inferedRatesHAR1_yesChpNpNoP_multi <- readRDS("/path/to/Results/Harringtonine/FirstRun_FullModel/inferedRatesHAR1_YesChpNpNoP_multi.rds")
inferedRatesHAR2_yesChpNpNoP_multi <- readRDS("/path/to/Results/Harringtonine/FirstRun_FullModel/inferedRatesHAR2_YesChpNpNoP_multi.rds")
inferedRatesHARMerged_yesChpNpNoP_multi <- readRDS("/path/to/Results/Harringtonine/FirstRun_FullModel/inferedRatesHARMerged_YesChpNpNoP_multi.rds")

## LepB
inferedRatesLEP1_yesChpNpP_multi <- readRDS("/path/to/Results/LeptomycinB/FirstRun_FullModel/inferedRatesLepB1_YesChpNpP_multi.rds")
inferedRatesLEP2_yesChpNpP_multi <- readRDS("/path/to/Results/LeptomycinB/FirstRun_FullModel/inferedRatesLepB2_YesChpNpP_multi.rds")
inferedRatesLEPMerged_yesChpNpP_multi <- readRDS("/path/to/Results/LeptomycinB/FirstRun_FullModel/inferedRatesLepBMerged_YesChpNpP_multi.rds")

### Replicates
## Expression levels - Figures S7, S18, S23, S28
replicatesDataCorrelationsUntreated1 <- unlist(replicatesDataCorrelations(object1=inferedRatesUntreated1_yesChpNpP_multi,object2=inferedRatesUntreated2_yesChpNpP_multi,expressionFlag=1,width=7,height=12,name="Untreated",lowSat=0.05,upSat=0.95))
round(c(min(replicatesDataCorrelationsUntreated1),max(replicatesDataCorrelationsUntreated1),median(replicatesDataCorrelationsUntreated1)),2)
replicatesDataCorrelationsPlaB1 <- unlist(replicatesDataCorrelations(object1=inferedRatesPlaB1_yesChpNpP_multi,object2=inferedRatesPlaB2_yesChpNpP_multi,expressionFlag=1,width=7,height=12,name="PlaB",lowSat=0.05,upSat=0.95))
round(c(min(replicatesDataCorrelationsPlaB1),max(replicatesDataCorrelationsPlaB1),median(replicatesDataCorrelationsPlaB1)),2)
replicatesDataCorrelationsHAR1 <- unlist(replicatesDataCorrelations(object1=inferedRatesHAR1_yesChpNpNoP_multi,object2=inferedRatesHAR2_yesChpNpNoP_multi,expressionFlag=1,width=7,height=12,name="HAR",lowSat=0.05,upSat=0.95))
round(c(min(replicatesDataCorrelationsHAR1),max(replicatesDataCorrelationsHAR1),median(replicatesDataCorrelationsHAR1)),2)
replicatesDataCorrelationsLEP1 <- unlist(replicatesDataCorrelations(object1=inferedRatesLEP1_yesChpNpP_multi,object2=inferedRatesLEP2_yesChpNpP_multi,expressionFlag=1,width=7,height=12,name="LEP",lowSat=0.05,upSat=0.95))
round(c(min(replicatesDataCorrelationsLEP1),max(replicatesDataCorrelationsLEP1),median(replicatesDataCorrelationsLEP1)),2)

## Kinetic rates -  Figures S8, S19, S24, S29
replicatesRatesCorrelationsUntreated1 <- unlist(replicatesRatesCorrelations(object1=inferedRatesUntreated1_yesChpNpP_multi,object2=inferedRatesUntreated2_yesChpNpP_multi,width=7,height=12,name="Untreated",lowSat=0,upSat=1))
round(c(min(replicatesRatesCorrelationsUntreated1),max(replicatesRatesCorrelationsUntreated1),median(replicatesRatesCorrelationsUntreated1)),2)
replicatesRatesCorrelationsPlaB1 <- unlist(replicatesRatesCorrelations(object1=inferedRatesPlaB1_yesChpNpP_multi,object2=inferedRatesPlaB2_yesChpNpP_multi,width=7,height=12,name="PlaB",lowSat=0,upSat=1))
round(c(min(replicatesRatesCorrelationsPlaB1),max(replicatesRatesCorrelationsPlaB1),median(replicatesRatesCorrelationsPlaB1)),2)
replicatesRatesCorrelationsHAR1 <- unlist(replicatesRatesCorrelations(object1=inferedRatesHAR1_yesChpNpNoP_multi,object2=inferedRatesHAR2_yesChpNpNoP_multi,width=7,height=12,name="HAR",lowSat=0,upSat=1))
round(c(min(replicatesRatesCorrelationsHAR1),max(replicatesRatesCorrelationsHAR1),median(replicatesRatesCorrelationsHAR1)),2)
replicatesRatesCorrelationsLEP1 <- unlist(replicatesRatesCorrelations(object1=inferedRatesLEP1_yesChpNpP_multi,object2=inferedRatesLEP2_yesChpNpP_multi,width=7,height=12,name="LEP",lowSat=0,upSat=1))
round(c(min(replicatesRatesCorrelationsLEP1),max(replicatesRatesCorrelationsLEP1),median(replicatesRatesCorrelationsLEP1)),2)

### Goodness of fit -  Figures S9, S20
goodnessOfFitUntreatedMergedNoPoly <- unlist(GoodnessOfFit(inferedRates=inferedRatesUntreatedMerged_noPoly,expressionFlag="MEAN",width=7,height=12,name="Untreated_NoPoly",lowSat=0.025,upSat=0.975))
round(c(min(goodnessOfFitUntreatedMergedNoPoly),max(goodnessOfFitUntreatedMergedNoPoly),median(goodnessOfFitUntreatedMergedNoPoly)),2)
goodnessOfFitUntreatedMerged <- unlist(GoodnessOfFit(inferedRates=inferedRatesUntreatedMerged_yesChpNpP_multi,expressionFlag="MEAN",width=7,height=12,name="Untreated",lowSat=0.025,upSat=0.975))
round(c(min(goodnessOfFitUntreatedMerged),max(goodnessOfFitUntreatedMerged),median(goodnessOfFitUntreatedMerged)),2)
goodnessOfFitPlaBMerged <- unlist(GoodnessOfFit(inferedRates=inferedRatesPlaBMerged_yesChpNpP_multi,expressionFlag="MEAN",width=7,height=12,name="PlaB",lowSat=0.025,upSat=0.975))
round(c(min(goodnessOfFitPlaBMerged),max(goodnessOfFitPlaBMerged),median(goodnessOfFitPlaBMerged)),2)

### Distributions - Figure 2C and S12B
RatesDistributions(object=inferedRatesUntreatedMerged_yesChpNpP_multi,width=24,height=4,obj_name="Untreated",xlimTmp=c(-2,6)) 
RatesDistributions(object=inferedRatesUntreatedMerged_yesChpNpP_single,width=24,height=4,obj_name="Untreated_Single",xlimTmp=c(-2,6))

### Bimodality conservation
## Untreated
  intervals <- list(c(0.25,1),c(0,1),c(-0.5,1),c(0.25,1.25))
  names(intervals) <- paste0("k",2:5)
  commonGenes <- intersect(rownames(inferedRatesUntreated1_yesChpNpP_multi$inferedRates),rownames(inferedRatesUntreated2_yesChpNpP_multi$inferedRates))

  for(i in paste0("k",2:5))
  {
    y <- hist(log10(inferedRatesUntreated1_yesChpNpP_multi$inferedRates[commonGenes,i]),100)
    thUntreated1 <- optimize(approxfun(y$mids,y$density),interval=intervals[[i]])$minimum
    
    y <- hist(log10(inferedRatesUntreated2_yesChpNpP_multi$inferedRates[commonGenes,i]),100)
    thUntreated2 <- optimize(approxfun(y$mids,y$density),interval=intervals[[i]])$minimum
    
    print(i)
    print(round(thUntreated1,2))
    print(round(thUntreated2,2))
  
    print(round((sum(log10(inferedRatesUntreated1_yesChpNpP_multi$inferedRates[commonGenes,i])>thUntreated1&log10(inferedRatesUntreated2_yesChpNpP_multi$inferedRates[commonGenes,i])>thUntreated2)+
                 sum(log10(inferedRatesUntreated1_yesChpNpP_multi$inferedRates[commonGenes,i])<thUntreated1&log10(inferedRatesUntreated2_yesChpNpP_multi$inferedRates[commonGenes,i])<thUntreated2))/length(commonGenes),2))

    print(fisher.test(table(log10(inferedRatesUntreated1_yesChpNpP_multi$inferedRates[commonGenes,i])>thUntreated1
                           ,log10(inferedRatesUntreated2_yesChpNpP_multi$inferedRates[commonGenes,i])>thUntreated2))$p.value)
  }

## PlaB
  intervals <- list(c(0.25,1),c(0,1),c(-0.5,1),c(0.25,1.25))
  names(intervals) <- paste0("k",2:5)
  commonGenes <- intersect(rownames(inferedRatesPlaB1_yesChpNpP_multi$inferedRates),rownames(inferedRatesPlaB2_yesChpNpP_multi$inferedRates))

  for(i in paste0("k",2:5))
  {
    y <- hist(log10(inferedRatesPlaB1_yesChpNpP_multi$inferedRates[commonGenes,i]),100)
    thPlaB1 <- optimize(approxfun(y$mids,y$density),interval=intervals[[i]])$minimum
    
    y <- hist(log10(inferedRatesPlaB2_yesChpNpP_multi$inferedRates[commonGenes,i]),100)
    thPlaB2 <- optimize(approxfun(y$mids,y$density),interval=intervals[[i]])$minimum
    
    print(i)
    print(round(thPlaB1,2))
    print(round(thPlaB2,2))
  
    print(round((sum(log10(inferedRatesPlaB1_yesChpNpP_multi$inferedRates[commonGenes,i])>thPlaB1&log10(inferedRatesPlaB2_yesChpNpP_multi$inferedRates[commonGenes,i])>thPlaB2)+
                 sum(log10(inferedRatesPlaB1_yesChpNpP_multi$inferedRates[commonGenes,i])<thPlaB1&log10(inferedRatesPlaB2_yesChpNpP_multi$inferedRates[commonGenes,i])<thPlaB2))/length(commonGenes),2))
 
    print(fisher.test(table(log10(inferedRatesPlaB1_yesChpNpP_multi$inferedRates[commonGenes,i])>thPlaB1
                           ,log10(inferedRatesPlaB2_yesChpNpP_multi$inferedRates[commonGenes,i])>thPlaB2))$p.value)
  }

## LepB
  intervals <- list(c(-0.25,1),c(0.25,1),c(-1,1),c(0.25,1.25),c(1,3))
  names(intervals) <- paste0("k",c(2:5,7))
  commonGenes <- intersect(rownames(inferedRatesLEP1_yesChpNpP_multi$inferedRates),rownames(inferedRatesLEP2_yesChpNpP_multi$inferedRates))

  for(i in paste0("k",c(2:5,7)))
  {
    y <- hist(log10(inferedRatesLEP1_yesChpNpP_multi$inferedRates[commonGenes,i]),100)
    thLEP1 <- optimize(approxfun(y$mids,y$density),interval=intervals[[i]])$minimum
    
    y <- hist(log10(inferedRatesLEP2_yesChpNpP_multi$inferedRates[commonGenes,i]),100)
    thLEP2 <- optimize(approxfun(y$mids,y$density),interval=intervals[[i]])$minimum
    
    print(i)
    print(round(thLEP1,2))
    print(round(thLEP2,2))
  
    print(round((sum(log10(inferedRatesLEP1_yesChpNpP_multi$inferedRates[commonGenes,i])>thLEP1&log10(inferedRatesLEP2_yesChpNpP_multi$inferedRates[commonGenes,i])>thLEP2)+
                 sum(log10(inferedRatesLEP1_yesChpNpP_multi$inferedRates[commonGenes,i])<thLEP1&log10(inferedRatesLEP2_yesChpNpP_multi$inferedRates[commonGenes,i])<thLEP2))/length(commonGenes),2))
 
    print(fisher.test(table(log10(inferedRatesLEP1_yesChpNpP_multi$inferedRates[commonGenes,i])>thLEP1
                           ,log10(inferedRatesLEP2_yesChpNpP_multi$inferedRates[commonGenes,i])>thLEP2))$p.value)
  }

## Harr
  intervals <- list(c(-0.25,1),c(0,1),c(-1,1),c(0.25,1.25))
  names(intervals) <- paste0("k",2:5)
  commonGenes <- intersect(rownames(inferedRatesHAR1_yesChpNpNoP_multi$inferedRates),rownames(inferedRatesHAR2_yesChpNpNoP_multi$inferedRates))

  for(i in paste0("k",2:5))
  {
    y <- hist(log10(inferedRatesHAR1_yesChpNpNoP_multi$inferedRates[commonGenes,i]),100)
    thHAR1 <- optimize(approxfun(y$mids,y$density),interval=intervals[[i]])$minimum
    
    y <- hist(log10(inferedRatesHAR2_yesChpNpNoP_multi$inferedRates[commonGenes,i]),100)
    thHAR2 <- optimize(approxfun(y$mids,y$density),interval=intervals[[i]])$minimum
    
    print(i)
    print(round(thHAR1,2))
    print(round(thHAR2,2))
  
    print(round((sum(log10(inferedRatesHAR1_yesChpNpNoP_multi$inferedRates[commonGenes,i])>thHAR1&log10(inferedRatesHAR2_yesChpNpNoP_multi$inferedRates[commonGenes,i])>thHAR2)+
                 sum(log10(inferedRatesHAR1_yesChpNpNoP_multi$inferedRates[commonGenes,i])<thHAR1&log10(inferedRatesHAR2_yesChpNpNoP_multi$inferedRates[commonGenes,i])<thHAR2))/length(commonGenes),2))
 
    print(fisher.test(table(log10(inferedRatesHAR1_yesChpNpNoP_multi$inferedRates[commonGenes,i])>thHAR1
                           ,log10(inferedRatesHAR2_yesChpNpNoP_multi$inferedRates[commonGenes,i])>thHAR2))$p.value)
  }

### Statistical comparisons
## Synthesis Untreated vs PlaB
ks.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k1"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k1"])$p.value
wilcox.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k1"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k1"],alternative="g")$p.value

## Co-transcriptional splicing Untreated vs PlaB
ks.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k2"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k2"])$p.value
wilcox.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k2"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k2"],alternative="g")$p.value

## Mature chromatin RNA detachment Untreated vs PlaB
ks.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k3"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k3"])$p.value
wilcox.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k3"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k3"],alternative="g")$p.value

## Premature chromatin RNA detachment Untreated vs PlaB
ks.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k4"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k4"])$p.value
wilcox.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k4"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k4"],alternative="l")$p.value

## Post-transcriptional splicing Untreated vs PlaB
ks.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k5"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k5"])$p.value
wilcox.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k5"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k5"],alternative="l")$p.value

### Violins - Figures 3E 4E and 5E
differentialDistributionPlot(object1=inferedRatesUntreatedMerged_yesChpNpP_multi,object2=inferedRatesPlaBMerged_yesChpNpP_multi,name="PlaBvsUntreated",width=6,height=4,breaks=100,col1="palevioletred",col2="lightblue",name1="Control",name2="PlaB",pchMed="-")
differentialDistributionPlot(object1=inferedRatesUntreatedMerged_noPoly,object2=inferedRatesHARMerged_yesChpNpNoP_multi,name="HarrvsUntreated",width=6,height=4,breaks=100,col1="palevioletred",col2="lightblue",name1="Control",name2="Harr",pchMed="-")
differentialDistributionPlot(object1=inferedRatesUntreatedMerged_yesChpNpP_multi,object2=inferedRatesLEPMerged_yesChpNpP_multi,name="LepBvsUntreated",width=6,height=4,breaks=100,col1="palevioletred",col2="lightblue",name1="Control",name2="LepB",pchMed="-")

### Heatmaps
## Untreated - Figure 2D
HeatmapsPlotUntreated <- HeatmapsPlot(object=inferedRatesUntreatedMerged_yesChpNpP_multi
                                     ,obj_name="Untreated",n_clust=6,ratesWeight=3,expressionFlag="MEAN")

## PlaB - Figures 3F and 3G
HeatmapsPlotPlaB <- HeatmapsPlot(object=inferedRatesPlaBMerged_yesChpNpP_multi
                                ,obj_name="PlaB",n_clust=7,ratesWeight=3,expressionFlag="MEAN")

DifferentialHeatmapsPlotPlaB <- DifferentialHeatmapsPlot(refObject=inferedRatesUntreatedMerged_yesChpNpP_multi
                                                        ,object=inferedRatesPlaBMerged_yesChpNpP_multi
                                                        ,name="PlaBvsUntreated",n_clust=6,ylim=c(-2,2),ratesWeight=3,expressionFlag="MEAN",width=5,height=4)

# Genes switching from co- to post- transcriptional regulations
commonGenesTmp <- intersect(rownames(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates),rownames(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates))

sum(log2(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k2"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k2"])<0
   &log2(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k3"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k3"])<0
   &log2(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k4"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k4"])>0
   &log2(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k5"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k5"])>0)

sum(log2(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k2"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k2"])>0
   &log2(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k3"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k3"])>0
   &log2(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k4"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k4"])<0
   &log2(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k5"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k5"])<0)

## HAR - Figures 5C and 5D
HeatmapsPlotHAR <- HeatmapsPlot(object=inferedRatesHARMerged_yesChpNpNoP_multi
                               ,obj_name="HAR",n_clust=7,ratesWeight=3,expressionFlag="MEAN")

DifferentialHeatmapsPlotHAR <- DifferentialHeatmapsPlot(refObject=inferedRatesUntreatedMerged_noPoly
                                                       ,object=inferedRatesHARMerged_yesChpNpNoP_multi
                                                       ,name="HARvsUntreated",n_clust=6,ylim=c(-2,2),ratesWeight=3,expressionFlag="MEAN",width=5,height=4)

commonGenesTmp <- intersect(rownames(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates),rownames(inferedRatesUntreatedMerged_noPoly$inferedRates))
sum(log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k1"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k1"])>0)/length(commonGenesTmp)

# Genes switching from co- to post- transcriptional regulations
sum(log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k2"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k2"])>0
   &log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k3"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k3"])>0
   &log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k4"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k4"])<0
   &log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k5"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k5"])<0)/length(commonGenesTmp)

sum(log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k2"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k2"])<0
   &log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k3"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k3"])<0
   &log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k4"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k4"])>0
   &log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k5"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k5"])>0)/length(commonGenesTmp)

## LEP - Figures 4B and 4D
HeatmapsPlotLEP <- HeatmapsPlot(object=inferedRatesLEPMerged_yesChpNpP_multi
                               ,obj_name="LEP",n_clust=7,ratesWeight=3,expressionFlag="MEAN")

DifferentialHeatmapsPlotLEP <- DifferentialHeatmapsPlot(refObject=inferedRatesUntreatedMerged_yesChpNpP_multi
													   ,object=inferedRatesLEPMerged_yesChpNpP_multi
													   ,name="LEPvsUntreated",n_clust=9,ylim=c(-2,2),ratesWeight=3,expressionFlag="MEAN",width=5,height=4)
table(DifferentialHeatmapsPlotLEP)

# Genes switching from co- to post- transcriptional regulations
commonGenesTmp <- intersect(rownames(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates),rownames(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates))

sum(log2(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k2"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k2"])>0
   &log2(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k3"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k3"])>0
   &log2(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k4"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k4"])<0
   &log2(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k5"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k5"])<0)/length(commonGenesTmp)

sum(log2(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k2"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k2"])<0
   &log2(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k3"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k3"])<0
   &log2(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k4"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k4"])>0
   &log2(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k5"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k5"])>0)/length(commonGenesTmp)

# Genes involved in other regulations
sum(log2(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k7"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k7"])<0)/length(commonGenesTmp)
sum(log2(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k8"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k8"])<0)/length(commonGenesTmp)
sum(log2(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k10"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k10"])>0)/length(commonGenesTmp)

# Export modulations
## Down in export
# Genes from Engel et al
EngelGenes <- read.table("/path/to/TableS4_LepBTargetGenes.txt",header=TRUE,sep="\t")
EngelGenesDown <- EngelGenes[EngelGenes[,6]<0.1&EngelGenes[,5]>0,1]
EngelGenesUp <- EngelGenes[EngelGenes[,6]<0.1&EngelGenes[,5]<0,1]

EngelGenesDown <- getBM(attributes = c("ensembl_gene_id","entrezgene_id"),filter="ensembl_gene_id",values=EngelGenesDown,mart = useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl")))
EngelGenesDown <- setdiff(unique(EngelGenesDown[,2]),NA)

EngelGenesUp <- getBM(attributes = c("ensembl_gene_id","entrezgene_id"),filter="ensembl_gene_id",values=EngelGenesUp,mart = useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl")))
EngelGenesUp <- setdiff(unique(EngelGenesUp[,2]),NA)

length(intersect(c(EngelGenesDown,EngelGenesUp),rownames(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates)))
length(intersect(c(EngelGenesDown,EngelGenesUp),rownames(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates)))
length(intersect(c(EngelGenesDown,EngelGenesUp),intersect(rownames(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates),rownames(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates))))

# Our genes
refObject=inferedRatesUntreatedMerged_yesChpNpP_multi
object=inferedRatesLEPMerged_yesChpNpP_multi

ExampleData1 <- lapply(refObject$expressionData,function(i){
  i[i==1e-10] <- NA
  i
})
ExampleData1 <- sapply(colnames(ExampleData1[[1]]),function(i)apply(cbind(ExampleData1[[1]][,i],ExampleData1[[2]][,i]),1,mean,na.rm=TRUE))
polyFlag1 <- ("pyt"%in%colnames(ExampleData1))
if(polyFlag1)
{
ExampleData1 <- ExampleData1[,c("chpt","chmt","npt","nmt","cyt","pyt")]
}else{
ExampleData1 <- ExampleData1[,c("chpt","chmt","npt","nmt","cyt")]
}

inferedRates1 <- refObject$inferedRates

ExampleData2 <- lapply(object$expressionData,function(i){
      i[i==1e-10] <- NA
      i
    })
    ExampleData2 <- sapply(colnames(ExampleData2[[1]]),function(i)apply(cbind(ExampleData2[[1]][,i],ExampleData2[[2]][,i]),1,mean,na.rm=TRUE))

polyFlag2 <- ("pyt"%in%colnames(ExampleData2))
if(polyFlag2)
{
ExampleData2 <- ExampleData2[,c("chpt","chmt","npt","nmt","cyt","pyt")]
}else{
ExampleData2 <- ExampleData2[,c("chpt","chmt","npt","nmt","cyt")]
}

inferedRates2 <- object$inferedRates

commonGenes <- intersect(rownames(ExampleData1),rownames(ExampleData2))
commonRates <- intersect(colnames(inferedRates1),colnames(inferedRates2))
commonSpecies <- intersect(colnames(ExampleData1),colnames(ExampleData2))

ExampleData1 <- ExampleData1[commonGenes,commonSpecies]
inferedRates1 <- inferedRates1[commonGenes,commonRates]
ExampleData2 <- ExampleData2[commonGenes,commonSpecies]
inferedRates2 <- inferedRates2[commonGenes,commonRates]

ExampleData <- log2(ExampleData2/ExampleData1)  
inferedRates <- log2(inferedRates2/inferedRates1)  

ExampleData[ExampleData<(-2)] <- (-2)
inferedRates[inferedRates<(-2)] <- (-2)
ExampleData[ExampleData>(2)] <- (2)
inferedRates[inferedRates>(2)] <- (2)

OurGenesDown <- names(which(inferedRates[,"k6"]<0))
OurGenesUp <- names(which(inferedRates[,"k6"]>0))

OurGenesHighDown <- names(which(inferedRates[,"k6"]<(-log2(2.5))))
OurGenesHighUp <- names(which(inferedRates[,"k6"]>(log2(2.5))))

OurGenesHighDownSymbol <- getBM(attributes = c("entrezgene_id","ensembl_gene_id","external_gene_name"),filter="entrezgene_id",values=OurGenesHighDown,mart = useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl")))
rownames(OurGenesHighDownSymbol) <- OurGenesHighDownSymbol[,1]
OurGenesHighDownSymbol <- dplyr::select(as.data.frame(OurGenesHighDownSymbol),columns='external_gene_name')

setdiff(OurGenesHighDown,rownames(OurGenesHighDownSymbol))
# "entrezid""

OurGenesHighDownSymbol <- rbind(OurGenesHighDownSymbol,data.frame(columns="GeneName"))
rownames(OurGenesHighDownSymbol)[nrow(OurGenesHighDownSymbol)] <- "entrezid"

OurGenesHighUpSymbol <- getBM(attributes = c("entrezgene_id","ensembl_gene_id","external_gene_name"),filter="entrezgene_id",values=OurGenesHighUp,mart = useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl")))
rownames(OurGenesHighUpSymbol) <- OurGenesHighUpSymbol[,1]
OurGenesHighUpSymbol <- dplyr::select(as.data.frame(OurGenesHighUpSymbol),columns='external_gene_name')

setdiff(OurGenesHighUp,rownames(OurGenesHighUpSymbol))

OurGenesHighSymbol <- rbind(OurGenesHighDownSymbol,OurGenesHighUpSymbol)

functionTmp <- function(OurGenesHighSymbol,name,n_clust,show_rownames)
{
  inferedRatesTmp <- inferedRatesLEPMerged_yesChpNpP_multi
  inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp <- inferedRatesUntreatedMerged_yesChpNpP_multi

  inferedRatesTmp$inferedRates <- inferedRatesTmp$inferedRates[rownames(OurGenesHighSymbol),]
  inferedRatesTmp$inferedData <- inferedRatesTmp$inferedData[rownames(OurGenesHighSymbol),]
  inferedRatesTmp$expressionData <- lapply(inferedRatesTmp$expressionData,function(i){
        qwe <- i[rownames(OurGenesHighSymbol),]
        rownames(qwe) <- OurGenesHighSymbol[rownames(qwe),1]
        qwe
    })
 
  rownames(inferedRatesTmp$inferedRates) <- OurGenesHighSymbol[rownames(inferedRatesTmp$inferedRates),1]
  rownames(inferedRatesTmp$inferedData) <- OurGenesHighSymbol[rownames(inferedRatesTmp$inferedData),1]
  
  inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedRates <- inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedRates[rownames(OurGenesHighSymbol),]
  inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedData <- inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedData[rownames(OurGenesHighSymbol),]
  inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$expressionData <- lapply(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$expressionData,function(i){i[rownames(OurGenesHighSymbol),]})
  
  rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedRates) <- OurGenesHighSymbol[rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedRates),1]
  rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedData) <- OurGenesHighSymbol[rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedData),1]
  rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$expressionData[[1]]) <- OurGenesHighSymbol[rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$expressionData[[1]]),1]
  rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$expressionData[[2]]) <- OurGenesHighSymbol[rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$expressionData[[2]]),1]
  
  DifferentialHeatmapsPlotLEPExtreme <- DifferentialHeatmapsPlot(refObject=inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp
                                      ,object=inferedRatesTmp
                                      ,expressionFlag="MEAN"
                                      ,name=name,n_clust=n_clust,ylim=c(-2,2),ratesWeight=3,show_rownames=show_rownames,width=10,height=10)
}

# Heatmaps - Figure 4C
functionTmp(OurGenesHighDownSymbol,name="LEPvsUntreated_Extreme_Down",n_clust=2,show_rownames=TRUE)
functionTmp(OurGenesHighUpSymbol,name="LEPvsUntreated_Extreme_Up",n_clust=2,show_rownames=TRUE)

### Untreated rates in Pladienolide B derived cluster - Figure 3H
df <- data.frame(cbind(DifferentialHeatmapsPlotPlaB,inferedRatesUntreatedMerged_yesChpNpP_multi[names(DifferentialHeatmapsPlotPlaB),][,1:5]))
par(mfrow=c(1,5))
sapply(colnames(df)[-1],function(i){
  boxplot(as.numeric(df[[i]])~df$ord,las=2,main=i,notch=TRUE,ylab="Log10(rates)",outline=FALSE)
})

### k1 vs k9 Log2FCs between Pladienolide B and Untreated cells - Figure S30
commonGenes <- intersect(rownames(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates),rownames(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates))

ratesUntreated <- inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenes,c("k1","k10")]
ratesTreated <- inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[commonGenes,c("k1","k10")]
ratesLog2FC <- log2(ratesTreated/ratesUntreated)

smoothScatter(ratesLog2FC[,"k10"],ratesLog2FC[,"k1"],ylab="Synthesis [Log2 FC]",xlab="Polysomal Degradation [Log2 FC]")
abline(0,1,lwd=2,col=2);abline(h=0);abline(v=0)
legend("bottomright",box.lwd=-1,legend=paste0("Cor.=",round(cor(ratesLog2FC[,"k1"],ratesLog2FC[,"k10"],method="s"),2)))

### GO Enrichments
## Untreated
# Clusters
GO_ALL_Untreated <- lapply(sapply(split(HeatmapsPlotUntreated,HeatmapsPlotUntreated),names),function(i)
{
        enrichGO(gene = i,
				 OrgDb = org.Hs.eg.db,
				 keyType = 'ENTREZID',
				 ont = "ALL",
				 pAdjustMethod = "BH",
				 universe = names(HeatmapsPlotUntreated),
				 qvalueCutoff = 0.1,
				 readable = TRUE,
				 maxGSSize = 3000)
})

# Top/Bottom genes
sortedRatesUntreated <- apply(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates,2,function(i)
{
  names(sort(i,decreasing=TRUE))
})

GO_ALL_Untreated_TopBottom <- apply(sortedRatesUntreated,2,function(i)
{

  list(enrichGO(gene = i[1:100],
                OrgDb = org.Hs.eg.db,
                keyType = 'ENTREZID',
                ont = "ALL",
                pAdjustMethod = "BH",
                universe = i[101:(length(i)-100)],
                qvalueCutoff = 0.1,
                readable = TRUE,
                maxGSSize = 3000),
       enrichGO(gene = i[(length(i)-99):length(i)],
                OrgDb = org.Hs.eg.db,
                keyType = 'ENTREZID',
                ont = "ALL",
                pAdjustMethod = "BH",
                universe = i[101:(length(i)-100)],
                qvalueCutoff = 0.1,
                readable = TRUE,
                maxGSSize = 3000))
})

## PlaB
GO_ALL_PlaB <- lapply(sapply(split(HeatmapsPlotPlaB,HeatmapsPlotPlaB),names),function(i)
{
        enrichGO(gene = i,
				 OrgDb = org.Hs.eg.db,
				 keyType = 'ENTREZID',
				 ont = "ALL",
				 pAdjustMethod = "BH",
				 universe = names(HeatmapsPlotPlaB),
				 qvalueCutoff = 0.1,
				 readable = TRUE,
				 maxGSSize = 3000)
})

GO_ALL_PlaB_DIFF <- lapply(sapply(split(DifferentialHeatmapsPlotPlaB,DifferentialHeatmapsPlotPlaB),names),function(i)
{
        enrichGO(gene = i,
				 OrgDb = org.Hs.eg.db,
				 keyType = 'ENTREZID',
				 ont = "ALL",
				 pAdjustMethod = "BH",
				 universe = names(DifferentialHeatmapsPlotPlaB),
				 qvalueCutoff = 0.1,
				 readable = TRUE,
         maxGSSize = 3000)
})

## HAR
GO_ALL_HAR <- lapply(sapply(split(HeatmapsPlotHAR,HeatmapsPlotHAR),names),function(i)
{
        enrichGO(gene = i,
				 OrgDb = org.Hs.eg.db,
				 keyType = 'ENTREZID',
				 ont = "ALL",
				 pAdjustMethod = "BH",
				 universe = names(HeatmapsPlotHAR),
				 qvalueCutoff = 0.1,
				 readable = TRUE,
         maxGSSize = 3000)
})

GO_ALL_HAR_DIFF <- lapply(sapply(split(DifferentialHeatmapsPlotHAR,DifferentialHeatmapsPlotHAR),names),function(i)
{
        enrichGO(gene = i,
				 OrgDb = org.Hs.eg.db,
				 keyType = 'ENTREZID',
				 ont = "ALL",
				 pAdjustMethod = "BH",
				 universe = names(DifferentialHeatmapsPlotHAR),
				 qvalueCutoff = 0.1,
				 readable = TRUE,
         maxGSSize = 3000)
})

## LEP
GO_ALL_LEP <- lapply(sapply(split(HeatmapsPlotLEP,HeatmapsPlotLEP),names),function(i)
{
        enrichGO(gene = i,
				 OrgDb = org.Hs.eg.db,
				 keyType = 'ENTREZID',
				 ont = "ALL",
				 pAdjustMethod = "BH",
				 universe = names(HeatmapsPlotLEP),
				 qvalueCutoff = 0.1,
				 readable = TRUE,
         maxGSSize = 3000)
})

GO_ALL_LEP_DIFF <- lapply(sapply(split(DifferentialHeatmapsPlotLEP,DifferentialHeatmapsPlotLEP),names),function(i)
{
        enrichGO(gene = i,
				 OrgDb = org.Hs.eg.db,
				 keyType = 'ENTREZID',
				 ont = "ALL",
				 pAdjustMethod = "BH",
				 universe = names(DifferentialHeatmapsPlotLEP),
				 qvalueCutoff = 0.1,
				 readable = TRUE,
         maxGSSize = 3000)
})

### Genes structural characterization.
## Longest UTRs sequences
cgUtr3 <- cgContentFunction("/path/to/Files/utr3.fa")
entropyUtr3 <- entropyFunction("/path/to/Files/utr3.fa")

cgUtr5 <- cgContentFunction("/path/to/Files/utr5.fa")
entropyUtr5 <- entropyFunction("/path/to/Files/utr5.fa")

## Figure S14
GenesCharacterization(txdb,HeatmapsPlotUntreated,obj_name="Untreated",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12)
GenesCharacterization(txdb,HeatmapsPlotPlaB,obj_name="PlaB",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12)
GenesCharacterization(txdb,DifferentialHeatmapsPlotPlaB,obj_name="PlaB_DIFF",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12)
GenesCharacterization(txdb,HeatmapsPlotHAR,obj_name="HAR",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12)
GenesCharacterization(txdb,DifferentialHeatmapsPlotHAR,obj_name="HAR_DIFF",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12)
GenesCharacterization(txdb,HeatmapsPlotLEP,obj_name="LEP",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12)
GenesCharacterization(txdb,DifferentialHeatmapsPlotLEP,obj_name="LEP_DIFF",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12)

## Untreated - Figure S11
GenesCharacterization(txdb=txdb,sortedRates=sortedRatesUntreated,obj_name="Untreated",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),top=100,bottom=100,width=12,height=6)

### PolyA tails characterization
## Untreated - Figure 2E, Figure S13
tailsUntreated <- TailsCharacterization(txdb=txdb,clustering=HeatmapsPlotUntreated,obj_name="Untreated"
                                       ,bam_untreated=NULL,bam_treated=bamPaths_Untreated,tailsPaths_untreatetd=NULL,tailsPaths_treated=tailsPaths_Untreated
                                       ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)

## Pladienolide B - Figure 3I, Figure S22 
tailsPlaB <- TailsCharacterization(txdb=txdb,clustering=HeatmapsPlotPlaB,obj_name="PlaB"
                                  ,bam_untreated=NULL,bam_treated=bamPaths_PlaB,tailsPaths_untreatetd=NULL,tailsPaths_treated=tailsPaths_PlaB
                                  ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)
tailsPlaBvsUntreated <- TailsCharacterization(txdb=txdb,clustering=DifferentialHeatmapsPlotPlaB,obj_name="PlaBvsUntreated"
                                      ,bam_untreated=bamPaths_Untreated,bam_treated=bamPaths_PlaB
                                      ,tailsPaths_untreatetd=tailsPaths_Untreated,tailsPaths_treated=tailsPaths_PlaB
                                      ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)

## Harringtonine - Figure 5F
tailsHAR <- TailsCharacterization(txdb=txdb,clustering=HeatmapsPlotHAR,obj_name="HAR"
                                ,bam_untreated=NULL,bam_treated=bamPaths_HAR,tailsPaths_untreatetd=NULL,tailsPaths_treated=tailsPaths_HAR
                                ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)

## Leptomycin B - Figure 4G
tailsLEP <- TailsCharacterization(txdb=txdb,clustering=HeatmapsPlotLEP,obj_name="LEP"
                                ,bam_untreated=NULL,bam_treated=bamPaths_LEP,tailsPaths_untreatetd=NULL,tailsPaths_treated=tailsPaths_LEP
                                ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)

tailsLEPvsUntreated <- TailsCharacterization(txdb=txdb,clustering=DifferentialHeatmapsPlotLEP,obj_name="LepBvsUntreated"
                                      ,bam_untreated=bamPaths_Untreated,bam_treated=bamPaths_LEP
                                      ,tailsPaths_untreatetd=tailsPaths_Untreated,tailsPaths_treated=tailsPaths_LEP
                                      ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)

## Correlations between tails and rates modulations
commonGenesPlaB <- intersect(rownames(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates),rownames(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates))
differentialRatesPlaB <- log2(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[commonGenesPlaB,]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesPlaB,])
tailsFCPlaB <- log2(sapply(tailsPlaBvsUntreated$Treated[[1]][[3]][commonGenesPlaB],median)/sapply(tailsPlaBvsUntreated$Untreated[[1]][[3]][commonGenesPlaB],median))

commonGenesLEP <- intersect(rownames(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates),rownames(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates))
differentialRatesLEP <- log2(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[commonGenesLEP,]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesLEP,])
tailsFCLEP <- log2(sapply(tailsLEPvsUntreated$Treated[[1]][[3]][commonGenesLEP],median)/sapply(tailsLEPvsUntreated$Untreated[[1]][[3]][commonGenesLEP],median))

commonGenesHAR <- intersect(rownames(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates),rownames(inferedRatesUntreatedMerged_noPoly$inferedRates))
differentialRatesHAR <- log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesHAR,]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesHAR,])
tailsFCHAR <- log2(sapply(tailsHARvsUntreated$Treated[[1]][[3]][commonGenesHAR],median)/sapply(tailsHARvsUntreated$Untreated[[1]][[3]][commonGenesHAR],median))

corMatrix <- significanceMatrix <- matrix(NaN,nrow=3,ncol=9)
colnames(corMatrix) <- colnames(significanceMatrix) <- colnames(differentialRatesPlaB)
rownames(corMatrix) <- rownames(significanceMatrix) <- c("Pladienolide B","Leptomycin B","Harringtonine")

corMatrix[1,] <- apply(differentialRatesPlaB,2,function(i)cor(i,tailsFCPlaB,method="s"))
significanceMatrix[1,] <- apply(differentialRatesPlaB,2,function(i)cor.test(i,tailsFCPlaB,method="s")$p.value)

corMatrix[2,] <- apply(differentialRatesLEP,2,function(i)cor(i,tailsFCLEP,method="s"))
significanceMatrix[2,] <- apply(differentialRatesLEP,2,function(i)cor.test(i,tailsFCLEP,method="s")$p.value)

corMatrix[3,colnames(differentialRatesHAR)] <- apply(differentialRatesHAR,2,function(i)cor(i,tailsFCHAR,method="s"))
significanceMatrix[3,colnames(differentialRatesHAR)] <- apply(differentialRatesHAR,2,function(i)cor.test(i,tailsFCHAR,method="s")$p.value)

colnames(corMatrix)[7:9] <- colnames(significanceMatrix)[7:9] <- c("k8","k7","k9")

corMatrix <- corMatrix[,sort(colnames(corMatrix))]
significanceMatrix <- significanceMatrix[,sort(colnames(significanceMatrix))]

corMatrixFinal <- corMatrix
corMatrixFinal[significanceMatrix>0.01] <- NaN

### k6 modulated genes sequence features - Figure S26
regulated <- c(OurGenesHighDown,OurGenesHighUp)
OurGenesNeutral <- inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[rownames(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates)%in%setdiff(rownames(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates),regulated),]
 
neutral <- cbind("Not-regulated",OurGenesNeutral[,"k6"])
up <- cbind("Up-regulated",inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[rownames(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates)%in%OurGenesHighUp,"k6"])
down <- cbind("Down-regulated",inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[rownames(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates)%in%OurGenesHighDown,"k6"])
df <- data.frame(rbind(down,neutral,up))
colnames(df) <- c("Class","Values")
 
clustering <- rownames(df)
names(clustering) <- df$Class
 
GenesCharacterization(txdb,clustering=clustering,obj_name="LEP_REG",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12)

### Polysomal RNA yield in Harringtonine - Figure 5B
df <- data.frame(cbind(c("UNT","UNT","HARR","HARR"),c(1,2,1,2),c(98/11,116.25/13,18/9,25/9)))
colnames(df) <- c("Condition","Replicate","Total")
df$Condition <- factor(df$Condition,levels = c("UNT","HARR"))
df$Total <- as.numeric(df$Total)

pl <- ggplot(df,aes(x = Condition,y = Total,fill = Replicate)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=myPal[4:7]) + ylab(label = "Total RNA extracted (ug)/ M.cells") +
  theme(legend.position="top") + theme( axis.title = element_text(size = 12), axis.text = element_text(size = 8)) 
ggsave("Polysomal.pdf",pl,width=4, height = 3)

### Premature Nucleoplasmic Export - Figure S39
DifferentialHeatmapsPlotLEP <- DifferentialHeatmapsPlot(refObject=inferedRatesUntreatedMerged_nucleoplasmicPrematureExport
                             ,object=inferedRatesPlaBMerged_nucleoplasmicPrematureExport
                             ,name="PlaBvsUntreated_NPE",n_clust=9,ylim=c(-2,2),ratesWeight=3)

### 5'TOP Factors - Figure S21
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

lookupUntreated <- getBM(mart = mart,
                         filters = "entrezgene_id",
                         attributes = c('entrezgene_id','ensembl_gene_id','gene_biotype','hgnc_symbol'),
                         values = rownames(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates))

lookupPlaB <- getBM(mart = mart,
                    filters = "entrezgene_id",
                    attributes = c('entrezgene_id','ensembl_gene_id','gene_biotype','hgnc_symbol'),
                    values = rownames(inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates))

lookupLepB <- getBM(mart = mart,
                    filters = "entrezgene_id",
                    attributes = c('entrezgene_id','ensembl_gene_id','gene_biotype','hgnc_symbol'),
                    values = rownames(inferedRatesLEPMerged_yesChpNpP_multi$inferedRates))

lookupUntreated$hgnc_symbol <- as.character(lookupUntreated$hgnc_symbol)
lookupPlaB$hgnc_symbol <- as.character(lookupPlaB$hgnc_symbol)
lookupLepB$hgnc_symbol <- as.character(lookupLepB$hgnc_symbol)

factors <- c("RPSA","RPS2","RPS3","RPS3A","RPS4X","RPS4Y","RPS5","RPS6","RPS7","RPS8","RPS9","RPS10","RPS11","RPS12","RPS13",
             "RPS14","RPS15","RPS15A","RPS16","RPS17","RPS18","RPS19","RPS20","RPS21","RPS23","RPS24","RPS25","RPS26","RPS27",
             "RPS27A","RPS28","RPS29","RPS30","RPP0","RPP1","RPP2","RPL3","RPL4","RPL5","RPL6","RPL7","RPL7A","RPL8","RPL9",
             "RPL10","RPL11","RPL10A","RPL12","RPL13","RPL13A","RPL14","RPL15","RPL17","RPL18","RPL18A","RPL19","RPL21","RPL22",
             "RPL23","RPL23A","RPL24","RPL26","RPL27","RPL27A","RPL30","RPL31","RPL32","RPL34","RPL35","RPL36","RPL36A","RPL37","RPL39","RPL40","RPL41")

results <- lapply(factors,function(i)
{
  untreated <- inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[which(lookupUntreated$hgnc_symbol==i),"k7"]
  plab <- inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[which(lookupPlaB$hgnc_symbol==i),"k7"]
  lepb <- inferedRatesLEPMerged_yesChpNpP_multi$inferedRates[which(lookupLepB$hgnc_symbol==i),"k7"]

  if(length(which(lookupUntreated$hgnc_symbol==i)>0)&&length(which(lookupPlaB$hgnc_symbol==i)>0)){
    FCwp <- log2(plab/untreated)
  }else{
    FCwp <- NA
  }
  if(length(which(lookupUntreated$hgnc_symbol==i)>0)&&length(which(lookupLepB$hgnc_symbol==i)>0)){
    FCwl <- log2(lepb/untreated)
  }else{
    FCwl <- NA
  }
  return(list("Untreated-PlaB"=FCwp,"Untreated-LepB"=FCwl))
})
names(results) <- factors

PlFC <- unlist(lapply(results,function(i){i[["Untreated-PlaB"]]}))
PlFC <- PlFC[!is.na(PlFC)]
LepFC <- unlist(lapply(results,function(i){i[["Untreated-LepB"]]}))
LepFC <- LepFC[!is.na(LepFC)]

par(mfrow=c(1,2))
barplot(PlFC,horiz=TRUE,ylab="Gene",xlab="Log2FC",col=c("blue","white","red")[sign(PlFC)+2],las=2)
barplot(LepFC,horiz=TRUE,ylab="Gene",xlab="Log2FC",col=c("blue","white","red")[sign(LepFC)+2],las=2)

### Gene example - Figure 2B
x <- inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates["2950",]
x <- x[c(1:6,8,7,9)]
names(x) <- paste0("k",1:9)

y <- apply(sapply(inferedRatesUntreatedMerged_yesChpNpP_multi$expressionData,function(i)i["2950",]),1,mean)
y <- y[grepl("t$",names(y))]
y <- y[c(2,1,5,4,3,6)]
names(y) <- c("Chp","Chm","Np","Nm","C","P")

y2 <- sapply(inferedRatesUntreatedMerged_yesChpNpP_multi$expressionData,function(i)i["2950",])
y2 <- y2[grepl("t$",rownames(y2)),]
y2 <- y2[c(2,1,5,4,3,6),]

z <- c(y,NaN,NaN,x[1],NaN,NaN,x[2:9])

barplotTmp <- barplot(z,ylim=c(1e-2,1e6),log="y",las=2,yaxt='n')
axis(2, at = 10^(-2:6),las=1)

points(barplotTmp[1:nrow(y2)],y2[,1],pch=16)
points(barplotTmp[1:nrow(y2)],y2[,2],pch=16)

### Couplings
## Compute the FC matrices: 0 weak correlation/anticorrelation, -1 anticorrelation, 0.75 correlation down-down, 1.25 correlation up-up
## Figure 6C and S31
plab_thr <- Mapping(mat=inferedRatesPlaBMerged_yesChpNpP_multi,refMat=inferedRatesUntreatedMerged_yesChpNpP_multi,thr=0.5,plot=FALSE,nclust=12,name="PlaB")
lepb_thr <- Mapping(mat=inferedRatesLEPMerged_yesChpNpP_multi,refMat=inferedRatesUntreatedMerged_yesChpNpP_multi,thr=0.5,plot=FALSE,nclust=10,name="LepB")
harr_thr <- Mapping(mat=inferedRatesHARMerged_yesChpNpNoP_multi,refMat=inferedRatesUntreatedMerged_noPoly,thr=0.5,plot=FALSE,nclust=7,name="Harr")

## Clusters
plab_clust <- plab_thr[[4]] 
lepb_clust <- lepb_thr[[4]]
harr_clust <- harr_thr[[4]] 

## Networks - Figures 6A and 6B
plab_corr <- PlotNetwork(matTmp=inferedRatesPlaBMerged_yesChpNpP_multi,refMat=inferedRatesUntreatedMerged_yesChpNpP_multi,name="PlaB",plot=TRUE)
lepb_corr <- PlotNetwork(matTmp=inferedRatesLEPMerged_yesChpNpP_multi,refMat=inferedRatesUntreatedMerged_yesChpNpP_multi,name="LepB",plot=TRUE)
harr_corr <- PlotNetwork(matTmp=inferedRatesHARMerged_yesChpNpNoP_multi,refMat=inferedRatesUntreatedMerged_noPoly,name="Harr",plot=TRUE)
 
intersection_corr <- PlotNetwork(matTmp=list(inferedRatesPlaBMerged_yesChpNpP_multi,inferedRatesLEPMerged_yesChpNpP_multi),refMat=inferedRatesUntreatedMerged_yesChpNpP_multi,name="Mean",plot=TRUE)

## Stacked barplot with colors
bar_p <- apply(plab_thr[[1]],2,function(i){c(zero=sum(i==0),up_pos=sum(i==1.25),up_neg=sum(i==0.75),anticorr=sum(i==-1))})
bar_l <- apply(lepb_thr[[1]],2,function(i){c(zero=sum(i==0),up_pos=sum(i==1.25),up_neg=sum(i==0.75),anticorr=sum(i==-1))}) 
bar_h <- apply(harr_thr[[1]],2,function(i){c(zero=sum(i==0),up_pos=sum(i==1.25),up_neg=sum(i==0.75),anticorr=sum(i==-1))})

PlotBar(bar_p,"PlaB_stacked")
PlotBar(bar_l,"LepB_stacked")
PlotBar(bar_h,"Harr_stacked")

## Tails characterization - Figures S32
tailsPlaB <- TailsCharacterization(txdb=txdb,clustering=plab_clust,obj_name="PlaBvsWT"
                                      ,bam_untreated=bamPaths_WT,bam_treated=bamPaths_PlaB
                                      ,tailsPaths_untreatetd=tailsPaths_WT,tailsPaths_treated=tailsPaths_PlaB
                                      ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)
tailsLepB <- TailsCharacterization(txdb=txdb,clustering=lepb_clust,obj_name="LepBvsWT"
                                      ,bam_untreated=bamPaths_WT,bam_treated=bamPaths_LepB
                                      ,tailsPaths_untreatetd=tailsPaths_WT,tailsPaths_treated=tailsPaths_LepB
                                      ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)
tailsHarr <- TailsCharacterization(txdb=txdb,clustering=harr_clust,obj_name="HarrvsWT"
                                      ,bam_untreated=bamPaths_WT,bam_treated=bamPaths_Harr
                                      ,tailsPaths_untreatetd=tailsPaths_WT,tailsPaths_treated=tailsPaths_Harr
                                      ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)

## Structural features - Figures S33
GenesCharacterization(txdb=txdb,clustering=plab_clust,obj_name="PlaB",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12,height=6,ids=LETTERS)
GenesCharacterization(txdb=txdb,clustering=lepb_clust,obj_name="LepB",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12,height=6,ids=LETTERS)
GenesCharacterization(txdb=txdb,clustering=harr_clust,obj_name="Harr",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12,height=6,ids=LETTERS)

## GO on clusters
plab_clust <- split(plab_clust,plab_clust)[unique(plab_clust)]
lepb_clust <- split(lepb_clust,lepb_clust)[unique(lepb_clust)] 
harr_clust <- split(harr_clust,harr_clust)[unique(harr_clust)]

GO_plab <- EnrichFunc(cluster=plab_clust,binaryMat=plab_thr[[1]],"PlaB")
GO_lepb <- EnrichFunc(cluster=lepb_clust,binaryMat=lepb_thr[[1]],"LepB")
GO_harr <- EnrichFunc(cluster=harr_clust,binaryMat=harr_thr[[1]],"Harr")

## Most similar and most different genes between Pladienolide B and Leptomycin B - Figures S34, S35
all_commons <- Reduce(intersect,list(rownames(plab_thr[[1]]),rownames(lepb_thr[[1]])))
PlaBc <- inferedRatesPlaBMerged_yesChpNpP_multi
PlaBc$inferedRates <- PlaBc$inferedRates[all_commons,]
LepBc <- inferedRatesLEPMerged_yesChpNpP_multi
LepBc$inferedRates <- LepBc$inferedRates[all_commons,]

plab_c <- Mapping(mat=PlaBc,refMat=inferedRatesUntreatedMerged_yesChpNpP_multi,thr=0.5,plot=FALSE,nclust=10,name="PlaB_commons")
lepb_c <- Mapping(mat=LepBc,refMat=inferedRatesUntreatedMerged_yesChpNpP_multi,thr=0.5,plot=FALSE,nclust=10,name="LepB_commons")
plabC <- plab_c[[1]][,intersect(colnames(plab_c[[1]]),colnames(lepb_c[[1]]))]
lepbC <- lepb_c[[1]][,intersect(colnames(plab_c[[1]]),colnames(lepb_c[[1]]))]

HammDist <- sapply(1:nrow(plabC),function(i){hamming.distance(plabC[i,],lepbC[i,])})
names(HammDist) <- rownames(plabC)
hist(HammDist,ylab="Hamming distance",xlab=NULL,breaks=30,main=NULL)

NullHamm1 <- NullModelHD(mat1=plabC,mat2=lepbC,nShuffle=1000,method=1)
NullHamm2 <- NullModelHD(mat1=plabC,mat2=lepbC,nShuffle=1000,method=2)
hist(NullHamm2,ylab="Hamming distance",xlab=NULL,breaks=30,main="Null model")

LowHamDist <- quantile(NullHamm2,probs = 0.05)
ConservedGenes <- names(HammDist[HammDist<=LowHamDist])
SymbConservedGenes <- getBM(attributes = c("hgnc_symbol","entrezgene_id"),filter="entrezgene_id",values=ConservedGenes,mart = useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl")))

PlabCons <- plab_thr[[1]][ConservedGenes,intersect(colnames(plab_thr[[1]]),colnames(lepb_thr[[1]]))]
colnames(PlabCons) <- sapply(colnames(PlabCons),function(i){paste0(i,"_plab")})
rownames(PlabCons) <- SymbConservedGenes[,1]

LepbCons <- lepb_thr[[1]][ConservedGenes,intersect(colnames(plab_thr[[1]]),colnames(lepb_thr[[1]]))]
colnames(LepbCons) <- sapply(colnames(LepbCons),function(i){paste0(i,"_lepb")})
rownames(LepbCons) <- SymbConservedGenes[,1]

ConsMat <- cbind(PlabCons,LepbCons)
pheatmap(ConsMat,cluster_cols=FALSE,color=c("navyblue","beige","darkolivegreen3","gold"),breaks=c(-1.5,-0.5,0.5,1,1.5),show_rownames=TRUE)

# 100 less conserved genes
HighHamDist <- sort(HammDist)[length(HammDist):(length(HammDist)-99)]
DiffGenes <- names(HighHamDist)
SymbDiffGenes <- getBM(attributes = c("hgnc_symbol","entrezgene_id"),filter="entrezgene_id",values=DiffGenes,mart = useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl")))

PlabDiff <- plab_thr[[1]][DiffGenes,intersect(colnames(plab_thr[[1]]),colnames(lepb_thr[[1]]))]
colnames(PlabDiff) <- sapply(colnames(PlabDiff),function(i){paste0(i,"_plab")})
rownames(PlabDiff) <- SymbDiffGenes[,1]
LepbDiff <- lepb_thr[[1]][DiffGenes,intersect(colnames(plab_thr[[1]]),colnames(lepb_thr[[1]]))]
colnames(LepbDiff) <- sapply(colnames(LepbDiff),function(i){paste0(i,"_lepb")})
rownames(LepbDiff) <- SymbDiffGenes[,1]

DiffMat <- cbind(PlabDiff,LepbDiff)
foe <- pheatmap(DiffMat,cluster_cols=FALSE,color=c("navyblue","beige","darkolivegreen3","gold"),breaks=c(-1.5,-0.5,0.5,1,1.5),show_rownames=TRUE)

# Enrichments of RBPs/TFs - GSEA analyses
WTrates <- inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates
WTrates_nopoly <- inferedRatesUntreatedMerged_noPoly$inferedRates
PlaBrates <- inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates; PlaBrates <- PlaBrates[intersect(rownames(WTrates),rownames(PlaBrates)),]
LepBrates <- inferedRatesLEPMerged_yesChpNpP_multi$inferedRates; LepBrates <- LepBrates[intersect(rownames(WTrates),rownames(LepBrates)),]
Harrates <- inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates; Harrates <- Harrates[intersect(rownames(WTrates_nopoly),rownames(Harrates)),]

# Matrices of log2FC
PlabFC <- log2(PlaBrates/WTrates[intersect(rownames(WTrates),rownames(PlaBrates)),])
LepbFC <- log2(LepBrates/WTrates[intersect(rownames(WTrates),rownames(LepBrates)),])
HarrFC <- log2(Harrates/WTrates_nopoly[intersect(rownames(WTrates_nopoly),rownames(Harrates)),])

# Reordering columns
colnames(PlabFC) <- c("k1","k2","k3","k4","k5","k6","k8","k7","k9")
colnames(LepbFC) <- c("k1","k2","k3","k4","k5","k6","k8","k7","k9")
colnames(HarrFC) <- c("k1","k2","k3","k4","k5","k6","k7")
PlabFC <- PlabFC[,sort(colnames(PlabFC))]
LepbFC <- LepbFC[,sort(colnames(LepbFC))]
HarrFC <- HarrFC[,sort(colnames(HarrFC))]

PlabEdges <- FoldChanges(PlabFC)
LepbEdges <- FoldChanges(LepbFC)
HarrEdges <- FoldChanges(HarrFC)

SortedPlabEdges <- lapply(colnames(PlabEdges),function(i)sort(PlabEdges[,i],decreasing=TRUE));names(SortedPlabEdges) <- colnames(PlabEdges)
SortedLepbEdges <- lapply(colnames(LepbEdges),function(i)sort(LepbEdges[,i],decreasing=TRUE));names(SortedLepbEdges) <- colnames(LepbEdges)
SortedHarrEdges <- lapply(colnames(HarrEdges),function(i)sort(HarrEdges[,i],decreasing=TRUE));names(SortedHarrEdges) <- colnames(HarrEdges)

# RNA Binding Proteins 
proteins <- list.files(path="/path/to/Encode/BED/Folder/RBPs",pattern = ".bed") # See encodeFilesRBPs.txt
metadata <- read.table("path/to/metadata.tsv",sep="\t",header=TRUE)
cell_type <- metadata$Biosample.term.name
prot_names <- metadata$Experiment.target

# Exons
exonsDB <- reduce(exonsBy(txdb ,'gene'))
exonsDB <- exonsDB[elementNROWS(range(exonsDB))==1]

# Build matrix
mat <- matrix(0,nrow=length(exonsDB),ncol=length(proteins))
rownames(mat) <- names(exonsDB)
colnames(mat) <- prot_names

# Fill matrix with the number of peaks per gene
for(i in 1:length(proteins)){
  bedTmp <- read.table(proteins[i],header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  names(bedTmp) <- c("seqnames","start","end","V4","V5","strand","V7","V8","V9")
  bedTmp <- makeGRangesFromDataFrame(bedTmp,
                           keep.extra.columns=FALSE,
                           ignore.strand=FALSE,
                           seqinfo=NULL,
                           seqnames.field=c("seqnames", "seqname",
                                            "chromosome", "chrom",
                                            "chr", "chromosome_name",
                                            "seqid"),
                           start.field="start",
                           end.field=c("end", "stop"),
                           strand.field="strand",
                           starts.in.df.are.0based=FALSE)
  overlap <- findOverlaps(bedTmp,unlist(range(exonsDB)),minoverlap=10) 
  foe <- sapply(split(overlap@from,overlap@to),length)
  names(foe) <- names(exonsDB)[as.numeric(names(foe))]
  mat[names(foe),i] <- foe
}

# Put together duplicated proteins
mat[mat<25]=0
mat[mat!=0]=1
qwe <- sapply(unique(colnames(mat[,duplicated(colnames(mat))])),function(i){
  x <- mat[,grep(i,colnames(mat))]
  x[rowSums(x)!=ncol(x),]=0
  x[rowSums(x)==ncol(x),]=1
  x[,1]
})
rbpMat <- cbind(mat[,!colnames(mat)%in%colnames(mat[,duplicated(colnames(mat))])],qwe)

# Transcription factors 
TF <- list.files(path="/path/to/Encode/BED/Folder/TFs",pattern = ".bed") # See encodeFilesTFs.txt
metadata <- read.table("metadata.tsv",sep="\t",header=TRUE)
cell_type <- metadata$Biosample.term.name
prot_names <- metadata$Experiment.target

# Exons
exonsDB <- reduce(exonsBy(txdb ,'gene'))
exonsDB <- exonsDB[elementNROWS(range(exonsDB))==1]
exonsDB <- unlist(range(exonsDB))
start(exonsDB)[strand(exonsDB)@values=="+"] <- start(exonsDB)[strand(exonsDB)@values=="+"] -2000
end(exonsDB)[strand(exonsDB)@values=="+"] <- end(exonsDB)[strand(exonsDB)@values=="+"] +1000
start(exonsDB)[strand(exonsDB)@values=="-"] <- start(exonsDB)[strand(exonsDB)@values=="-"] -1000
end(exonsDB)[strand(exonsDB)@values=="-"] <- end(exonsDB)[strand(exonsDB)@values=="-"] +2000

# Build matrix
mat <- matrix(0,nrow=length(exonsDB),ncol=length(TF))
rownames(mat) <- names(exonsDB)
colnames(mat) <- prot_names

# Fill matrix with the number of peaks per gene
for(i in 1:length(TF)){
  bedTmp <- read.table(TF[i],header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  names(bedTmp) <- c("seqnames","start","end","V4","V5","strand","V7","V8","V9")
  bedTmp <- makeGRangesFromDataFrame(bedTmp,
                           keep.extra.columns=FALSE,
                           ignore.strand=FALSE,
                           seqinfo=NULL,
                           seqnames.field=c("seqnames", "seqname",
                                            "chromosome", "chrom",
                                            "chr", "chromosome_name",
                                            "seqid"),
                           start.field="start",
                           end.field=c("end", "stop"),
                           strand.field="strand",
                           starts.in.df.are.0based=FALSE)
  overlap <- findOverlaps(bedTmp,exonsDB,minoverlap=10)
  foe <- sapply(split(overlap@from,overlap@to),length)
  names(foe) <- names(exonsDB)[as.numeric(names(foe))]
  mat[names(foe),i] <- foe
}

mat[mat<1] = 0
qwe <- sapply(unique(colnames(mat[,duplicated(colnames(mat))])),function(i){
  x <- mat[,grep(i,colnames(mat))]
  x[rowSums(x)!=ncol(x),]=0
  x[rowSums(x)==ncol(x),]=1
  x[,1]
})
tfMat<- cbind(mat[,!colnames(mat)%in%colnames(mat[,duplicated(colnames(mat))])],qwe)

## Pladienolide B enrichments
#RBPs
OurList <- sapply(colnames(rbpMat),function(i){
  geni <- names(which(rbpMat[rownames(PlabEdges),i]!=0))
  cbind(rep(i,length(geni)),geni)
})
OurList <- data.frame(do.call("rbind",OurList))

set.seed(1)
GSEA_PlaB <- lapply(SortedPlabEdges,function(i){
     GSEA(i,
          exponent = 1,
          minGSSize = 10,
          maxGSSize = 500,
          eps = 1e-10,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          gson = NULL,
          TERM2GENE = OurList,
          TERM2NAME = NA,
          verbose = TRUE,
          seed = TRUE,
          by = "fgsea"
          )
      })
names(GSEA_PlaB) <- names(SortedPlabEdges)
GSEA_PlaB <- GSEA_PlaB[unlist(lapply(GSEA_PlaB,function(i){length(i$ID)!=0}))]

# TFs
OurList <- sapply(colnames(tfMat),function(i){
  geni <- names(which(tfMat[rownames(PlabEdges),i]!=0))
  cbind(rep(i,length(geni)),geni)
})
OurList <- data.frame(do.call("rbind",OurList))

set.seed(1)
GSEA_PlaBTF <- lapply(SortedPlabEdges,function(i){
     GSEA(i,
          exponent = 1,
          minGSSize = 10,
          maxGSSize = 500,
          eps = 1e-10,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          gson = NULL,
          TERM2GENE = OurList,
          TERM2NAME = NA,
          verbose = TRUE,
          seed = TRUE,
          by = "fgsea"
          )
      })
names(GSEA_PlaBTF) <- names(SortedPlabEdges)
GSEA_PlaBTF <- GSEA_PlaBTF[unlist(lapply(GSEA_PlaBTF,function(i){length(i$ID)!=0}))]

## Leptomycin B enrichments
OurList <- sapply(colnames(rbpMat),function(i){
  geni <- names(which(rbpMat[rownames(LepbEdges),i]!=0))
  cbind(rep(i,length(geni)),geni)
})
OurList <- data.frame(do.call("rbind",OurList))

set.seed(1)
GSEA_LepB <- lapply(SortedLepbEdges,function(i){
     GSEA(i,
          exponent = 1,
          minGSSize = 10,
          maxGSSize = 500,
          eps = 1e-10,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          gson = NULL,
          TERM2GENE = OurList,
          TERM2NAME = NA,
          verbose = TRUE,
          seed = TRUE,
          by = "fgsea"
          )
      })
names(GSEA_LepB) <- names(SortedLepbEdges)
GSEA_LepB <- GSEA_LepB[unlist(lapply(GSEA_LepB,function(i){length(i$ID)!=0}))]

OurList <- sapply(colnames(tfMat),function(i){
  geni <- names(which(tfMat[rownames(LepbEdges),i]!=0))
  cbind(rep(i,length(geni)),geni)
})
OurList <- data.frame(do.call("rbind",OurList))

set.seed(1)
GSEA_LepBTF <- lapply(SortedLepbEdges,function(i){
     GSEA(i,
          exponent = 1,
          minGSSize = 10,
          maxGSSize = 500,
          eps = 1e-10,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          gson = NULL,
          TERM2GENE = OurList,
          TERM2NAME = NA,
          verbose = TRUE,
          seed = TRUE,
          by = "fgsea"
          )
      })
names(GSEA_LepBTF) <- names(SortedLepbEdges)
GSEA_LepBTF <- GSEA_LepBTF[unlist(lapply(GSEA_LepBTF,function(i){length(i$ID)!=0}))]

## Enrichments barplots - Figure 6D and S36
significant_PlaB = correlationSignificance(matTmp=inferedRatesPlaBMerged_yesChpNpP_multi,refMat=inferedRatesUntreatedMerged_yesChpNpP_multi,method="s")[[3]] 
significant_LepB = correlationSignificance(matTmp=inferedRatesLEPMerged_yesChpNpP_multi,refMat=inferedRatesUntreatedMerged_yesChpNpP_multi,method="s")[[3]] 

PlotBars(rbpMat = GSEA_PlaB)
         ,tfMat = GSEA_PlaBTF)
         ,name="PlaB"
         ,significant=significant_PlaB)
PlotBars(rbpMat = GSEA_LepB)
         ,tfMat = GSEA_LepBTF)
         ,name="LepB"
         ,significant=significant_LepB)
### Overfit check
OverfitFunctionTmp <- function(inferedData,expressionData,width=7,height=12,name="",lowSat=0,upSat=1,NpQ=NULL)
{
  commonGenesTmp <- intersect(rownames(inferedData),rownames(expressionData))
  
  inferedDataTmp <- inferedData[commonGenesTmp,]
  expressionDataTmp <- expressionData[commonGenesTmp,]

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

  if(!is.null(NpQ))
  {
    expressionDataTmp <- expressionDataTmp[which(rank(apply(cbind(rank(expressionDataTmp[,"Np_Nas_20min"]),rank(expressionDataTmp[,"Np_PreEx_20min"])),1,mean))>=NpQ*nrow(expressionDataTmp)),]
    #expressionDataTmp <- expressionDataTmp[expressionDataTmp[,"Np_Nas_20min"]>quantile(expressionDataTmp[,"Np_Nas_20min"],NpQ),]
    inferedDataTmp <- inferedDataTmp[rownames(expressionDataTmp),]    
  }

  print(paste0(nrow(inferedDataTmp)," genes left."))

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

Real1 <- inferedRatesUntreatedMerged_yesChpNpP_multi$expressionData[[1]]
Real2 <- inferedRatesUntreatedMerged_yesChpNpP_multi$expressionData[[2]]
Real <- rbind(Real1,Real2)
InfTmp <- inferedRatesUntreatedMerged_yesChpNpP_multi$inferedData[rownames(Real),]

Real[Real==1e-10] <- NaN

rownames(Real) <- rownames(InfTmp) <- 1:nrow(Real)

corReal1vsReal2 <- OverfitFunctionTmp(inferedData=inferedRatesUntreated1_yesChpNpP_multi$expressionData[[1]],expressionData=inferedRatesUntreated2_yesChpNpP_multi$expressionData[[1]],width=7,height=12,name="foe",lowSat=0,upSat=1)
corInfMvsReal <- OverfitFunctionTmp(inferedData=InfTmp,expressionData=Real,width=7,height=12,name="foe",lowSat=0,upSat=1)
corInf1vsReal2 <- OverfitFunctionTmp(inferedData=inferedRatesUntreated1_yesChpNpP_multi$inferedData,expressionData=inferedRatesUntreated2_yesChpNpP_multi$expressionData[[1]],width=7,height=12,name="foe",lowSat=0,upSat=1)
corInf2vsReal1 <- OverfitFunctionTmp(inferedData=inferedRatesUntreated2_yesChpNpP_multi$inferedData,expressionData=inferedRatesUntreated1_yesChpNpP_multi$expressionData[[1]],width=7,height=12,name="foe",lowSat=0,upSat=1)

corMatrix <- rbind(as.numeric(corInfMvsReal)
          ,as.numeric(corReal1vsReal2)
          ,as.numeric(corInf1vsReal2)
          ,as.numeric(corInf2vsReal1))
colnames(corMatrix) <- names(corInfMvsReal)

## Highly expressed
corReal1vsReal2_HE <- OverfitFunctionTmp(inferedData=inferedRatesUntreated1_yesChpNpP_multi$expressionData[[1]],expressionData=inferedRatesUntreated2_yesChpNpP_multi$expressionData[[1]],width=7,height=12,name="foe",lowSat=0,upSat=1,NpQ=0.9)
corInfMvsReal_HE <- OverfitFunctionTmp(inferedData=InfTmp,expressionData=Real,width=7,height=12,name="foe",lowSat=0,upSat=1,NpQ=0.9)
corInf1vsReal2_HE <- OverfitFunctionTmp(inferedData=inferedRatesUntreated1_yesChpNpP_multi$inferedData,expressionData=inferedRatesUntreated2_yesChpNpP_multi$expressionData[[1]],width=7,height=12,name="foe",lowSat=0,upSat=1,NpQ=0.9)
corInf2vsReal1_HE <- OverfitFunctionTmp(inferedData=inferedRatesUntreated2_yesChpNpP_multi$inferedData,expressionData=inferedRatesUntreated1_yesChpNpP_multi$expressionData[[1]],width=7,height=12,name="foe",lowSat=0,upSat=1,NpQ=0.9)

corMatrix_HE <- rbind(as.numeric(corInfMvsReal_HE)
          ,as.numeric(corReal1vsReal2_HE)
          ,as.numeric(corInf1vsReal2_HE)
          ,as.numeric(corInf2vsReal1_HE))
colnames(corMatrix_HE) <- names(corInfMvsReal_HE)

## Figure S40
par(mfrow=c(2,1))
par(mar = c(9, 4, 4, 2))

barplot(corMatrix,beside=TRUE,las=2,ylim=c(0,1.5),col=c("black","red","blue","darkgreen"),ylab="Spearman Correlation",main="New models",border=NULL)
legend("top",pch=15,col=c("black","red","blue","darkgreen"),legend=c("Mod_vs_Exp","Exp1_vs_Exp2","Mod1_vs_Exp2","Mod2_vs_Exp1"),ncol=2,box.lwd="-1")
abline(h=1,col=1,lwd=1,lty=3)

barplot(corMatrix_HE,beside=TRUE,las=2,ylim=c(0,1.5),col=c("black","red","blue","darkgreen"),ylab="Spearman Correlation",main="New models\nHigh Np Nascent",border=NULL)
legend("top",pch=15,col=c("black","red","blue","darkgreen"),legend=c("Mod_vs_Exp","Exp1_vs_Exp2","Mod1_vs_Exp2","Mod2_vs_Exp1"),ncol=2,box.lwd="-1")
abline(h=1,col=1,lwd=1,lty=3)

### Differential species - Figure S41
## Untreated 
UntreatedRep1vsRep2<-differentialSpecies(CIControl=expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT1
                                        ,CITreatment=expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT2
                                        ,nameTmp="Untreated Rep1 vs Rep2")

UntreatedRep1vsRep2_NoPoly<-differentialSpecies(CIControl=expressionLevelsUntreated_NoPoly$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$WT1
            ,CITreatment=expressionLevelsUntreated_NoPoly$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$WT2
            ,nameTmp="Untreated Rep1 vs Rep2")

## Pladienolide B
PladienolideBRep1vsRep2<-differentialSpecies(CIControl=expressionLevelsPlaB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$PlaB1
                                            ,CITreatment=expressionLevelsPlaB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$PlaB2
                                            ,nameTmp="PladienolideB Rep1 vs Rep2")

UntreatedRep1vsPladienolideBRep1<-differentialSpecies(CIControl=expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT1
            ,CITreatment=expressionLevelsPlaB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$PlaB1
            ,nameTmp="Untreated Rep1 vs PladienolideB Rep1")

UntreatedRep2vsPladienolideBRep2<-differentialSpecies(CIControl=expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT2
            ,CITreatment=expressionLevelsPlaB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$PlaB2
            ,nameTmp="Untreated Rep2 vs PladienolideB Rep2")

UntreatedRep2vsPladienolideBRep1<-differentialSpecies(CIControl=expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT2
            ,CITreatment=expressionLevelsPlaB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$PlaB1
            ,nameTmp="Untreated Rep2 vs PladienolideB Rep1")

UntreatedRep1vsPladienolideBRep2<-differentialSpecies(CIControl=expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT1
            ,CITreatment=expressionLevelsPlaB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$PlaB2
            ,nameTmp="Untreated Rep1 vs PladienolideB Rep2")

UntreatedRep2vsRep1<-differentialSpecies(CIControl=expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT2
                ,CITreatment=expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT1
                ,nameTmp="Untreated Rep2 vs Rep1")

## Leptomycin B
LeptomycinBRep2vsRep1<-differentialSpecies(CIControl=expressionLevelsLepB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$LepB2
            ,CITreatment=expressionLevelsLepB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$LepB1
            ,nameTmp="LeptomycinB Rep2 vs Rep1")

LeptomycinBRep1vsRep2<-differentialSpecies(CIControl=expressionLevelsLepB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$LepB1
            ,CITreatment=expressionLevelsLepB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$LepB2
            ,nameTmp="LeptomycinB Rep1 vs Rep2")

UntreatedRep1vsLeptomycinBRep1<-differentialSpecies(CIControl=expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT1
            ,CITreatment=expressionLevelsLepB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$LepB1
            ,nameTmp="Untreated Rep1 vs LeptomycinB Rep1")

UntreatedRep2vsLeptomycinBRep2<-differentialSpecies(CIControl=expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT2
            ,CITreatment=expressionLevelsLepB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$LepB2
            ,nameTmp="Untreated Rep2 vs LeptomycinB Rep2")

UntreatedRep2vsLeptomycinBRep1<-differentialSpecies(CIControl=expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT2
            ,CITreatment=expressionLevelsLepB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$LepB1
            ,nameTmp="Untreated Rep2 vs LeptomycinB Rep1")

UntreatedRep1vsLeptomycinBRep2<-differentialSpecies(CIControl=expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT1
            ,CITreatment=expressionLevelsLepB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$LepB2
            ,nameTmp="Untreated Rep1 vs LeptomycinB Rep2")

## Harringtonine
HarringtonineRep2vsRep1<-differentialSpecies(CIControl=expressionLevelsHAR$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$HAR2
            ,CITreatment=expressionLevelsHAR$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$HAR1
            ,nameTmp="Harringtonine Rep2 vs Rep1")

HarringtonineRep1vsRep2<-differentialSpecies(CIControl=expressionLevelsHAR$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$HAR1
            ,CITreatment=expressionLevelsHAR$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$HAR2
            ,nameTmp="Harringtonine Rep1 vs Rep2")

UntreatedRep1vsHarringtonineRep1<-differentialSpecies(CIControl=expressionLevelsUntreated_NoPoly$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$WT1
            ,CITreatment=expressionLevelsHAR$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$HAR1
            ,nameTmp="Untreated Rep1 vs Harringtonine Rep1")

UntreatedRep2vsHarringtonineRep2<-differentialSpecies(CIControl=expressionLevelsUntreated_NoPoly$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$WT2
            ,CITreatment=expressionLevelsHAR$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$HAR2
            ,nameTmp="Untreated Rep2 vs Harringtonine Rep2")

UntreatedRep2vsHarringtonineRep1<-differentialSpecies(CIControl=expressionLevelsUntreated_NoPoly$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$WT2
            ,CITreatment=expressionLevelsHAR$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$HAR1
            ,nameTmp="Untreated Rep2 vs Harringtonine Rep1")

UntreatedRep1vsHarringtonineRep2<-differentialSpecies(CIControl=expressionLevelsUntreated_NoPoly$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$WT1
            ,CITreatment=expressionLevelsHAR$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpNoP$HAR2
            ,nameTmp="Untreated Rep1 vs Harringtonine Rep2")
dev.off()

pdf("DifferentialSpeciesSingleBarplot.pdf",width=10,height=10)
par(mfrow=c(2,2))
barplot(rbind(Untreated_Rep1_vs_Rep2=apply(abs(UntreatedRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_Rep2=apply(abs(UntreatedRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_Rep2=apply(abs(UntreatedRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_Rep2=apply(abs(UntreatedRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,PladienolideB_Rep1_vs_Rep2=apply(abs(PladienolideBRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,PladienolideB_Rep1_vs_Rep2=apply(abs(PladienolideBRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,PladienolideB_Rep1_vs_Rep2=apply(abs(PladienolideBRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,PladienolideB_Rep1_vs_Rep2=apply(abs(PladienolideBRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_PladienolideB_Rep1=apply(abs(UntreatedRep1vsPladienolideBRep1),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep2_vs_PladienolideB_Rep2=apply(abs(UntreatedRep2vsPladienolideBRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep2_vs_PladienolideB_Rep1=apply(abs(UntreatedRep2vsPladienolideBRep1),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_PladienolideB_Rep2=apply(abs(UntreatedRep1vsPladienolideBRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i))))
,las=2,ylab="Regulated genes [%]",ylim=c(0,1),beside=TRUE,border=NaN,col=c(rep("red",4),rep("blue",4),rep("darkgreen",4)),main="PladienolideB")
legend("topleft",col=c("red","blue","darkgreen"),pch=15,box.lwd=-1,legend=c("Untreated replicates","PladienolideB replicates","Untreated vs PladienolideB"))

barplot(rbind(Untreated_Rep1_vs_Rep2=apply(abs(UntreatedRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_Rep2=apply(abs(UntreatedRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_Rep2=apply(abs(UntreatedRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_Rep2=apply(abs(UntreatedRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,LeptomycinB_Rep1_vs_Rep2=apply(abs(LeptomycinBRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,LeptomycinB_Rep1_vs_Rep2=apply(abs(LeptomycinBRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,LeptomycinB_Rep1_vs_Rep2=apply(abs(LeptomycinBRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,LeptomycinB_Rep1_vs_Rep2=apply(abs(LeptomycinBRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_LeptomycinB_Rep1=apply(abs(UntreatedRep1vsLeptomycinBRep1),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep2_vs_LeptomycinB_Rep2=apply(abs(UntreatedRep2vsLeptomycinBRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep2_vs_LeptomycinB_Rep1=apply(abs(UntreatedRep2vsLeptomycinBRep1),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_LeptomycinB_Rep2=apply(abs(UntreatedRep1vsLeptomycinBRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i))))
,las=2,ylab="Regulated genes [%]",ylim=c(0,1),beside=TRUE,border=NaN,col=c(rep("red",4),rep("blue",4),rep("darkgreen",4)),main="LeptomycinB")
legend("topleft",col=c("red","blue","darkgreen"),pch=15,box.lwd=-1,legend=c("Untreated replicates","LeptomycinB replicates","Untreated vs LeptomycinB"))

barplot(rbind(Untreated_Rep1_vs_Rep2=apply(abs(UntreatedRep1vsRep2_NoPoly),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_Rep2=apply(abs(UntreatedRep1vsRep2_NoPoly),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_Rep2=apply(abs(UntreatedRep1vsRep2_NoPoly),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_Rep2=apply(abs(UntreatedRep1vsRep2_NoPoly),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Harringtonine_Rep1_vs_Rep2=apply(abs(HarringtonineRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Harringtonine_Rep1_vs_Rep2=apply(abs(HarringtonineRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Harringtonine_Rep1_vs_Rep2=apply(abs(HarringtonineRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Harringtonine_Rep1_vs_Rep2=apply(abs(HarringtonineRep1vsRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_Harringtonine_Rep1=apply(abs(UntreatedRep1vsHarringtonineRep1),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep2_vs_Harringtonine_Rep2=apply(abs(UntreatedRep2vsHarringtonineRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep2_vs_Harringtonine_Rep1=apply(abs(UntreatedRep2vsHarringtonineRep1),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i)))
             ,Untreated_Rep1_vs_Harringtonine_Rep2=apply(abs(UntreatedRep1vsHarringtonineRep2),2,function(i)sum(i,na.rm=TRUE)/sum(is.finite(i))))
,las=2,ylab="Regulated genes [%]",ylim=c(0,1),beside=TRUE,border=NaN,col=c(rep("red",4),rep("blue",4),rep("darkgreen",4)),main="Harringtonine")
legend("topleft",col=c("red","blue","darkgreen"),pch=15,box.lwd=-1,legend=c("Untreated replicates","Harringtonine replicates","Untreated vs Harringtonine"))
dev.off()

### Genes modulated in k6 vs Differential expressions
j=2.5 # Fold Change

upRegulatedGenesTmp <- names(which(log2(inferedRatesLepBMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k6"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k6"])>=log2(j)))
downRegulatedGenesTmp <- names(which(log2(inferedRatesLepBMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k6"]/inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[commonGenesTmp,"k6"])<=(-log2(j))))
notRegulatedGenesTmp <- setdiff(commonGenesTmp,c(upRegulatedGenesTmp,downRegulatedGenesTmp))

i=2 # Number of modulated conditions
differentialGenesLepB <- sapply(colnames(differentialGenesListLepB[[1]]),function(i)apply(sapply(differentialGenesListLepB,function(j)j[,i]),1,function(z)sum(z,na.rm=TRUE)))
differentialGenesLepB[abs(differentialGenesLepB)<i] <- 0
differentialGenesLepB <- differentialGenesLepB[,grepl("^nm",colnames(differentialGenesLepB))|grepl("^cy",colnames(differentialGenesLepB))]
differentialGenesLepB <- differentialGenesLepB[,!grepl("0$",colnames(differentialGenesLepB))]
  
differentialGenesLepB <- sign(differentialGenesLepB[,"nmt"])

fisher.test(table(commonGenesTmp%in%downRegulatedGenesTmp,differentialGenesLepB==1))$p.value
fisher.test(table(commonGenesTmp%in%upRegulatedGenesTmp,differentialGenesLepB==(-1)))$p.value
fisher.test(table(commonGenesTmp%in%upRegulatedGenesTmp,differentialGenesLepB==1))$p.value
fisher.test(table(commonGenesTmp%in%downRegulatedGenesTmp,differentialGenesLepB==(-1)))$p.value,3
