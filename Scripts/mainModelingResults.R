### Libraries
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
library("DescTools")
library("vioplot")
require('biomaRt')

source("/path/to/allInternalFunctions.R")

### All BAMs
## Untreated
bamPaths_Untreated <- c(Chr_Untreated1="/path/to/bam.RData"
                       ,Nuc_Untreated1="/path/to/bam.RData"
                       ,Cyt_Untreated1="/path/to/bam.RData"
                       ,Poly_Untreated1="/path/to/bam.RData" # BAM merged from multiple samples.
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
                  ,Chr_PlaB2="/path/to/bam.RData" # BAM merged from multiple samples.
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
                 ,Chr_LEP2="/path/to/bam.RData" # BAM merged from multiple samples.
                 ,Nuc_LEP2="/path/to/bam.RData"
                 ,Cyt_LEP2="/path/to/bam.RData"
                 ,Poly_LEP2="/path/to/bam.RData"
                )

### All Tail Files
## Untreated
tailsPaths_Untreated <- c(Chr_Untreated1="path/to/polyAtails.tsv"
                         ,Nuc_Untreated1="path/to/polyAtails.tsv"
                         ,Cyt_Untreated1="path/to/polyAtails.tsv"
                         ,Poly_Untreated1="path/to/polyAtails.tsv"
                         ,Poly_Untreated1="path/to/polyAtails.tsv"
                         ,Chr_Untreated2="path/to/polyAtails.tsv"
                         ,Nuc_Untreated2="path/to/polyAtails.tsv"
                         ,Cyt_Untreated2="path/to/polyAtails.tsv"
                         ,Poly_Untreated2="path/to/polyAtails.tsv"
                         )

## PlaB
tailsPaths_PlaB <- c(Chr_PlaB1="path/to/polyAtails.tsv"
                    ,Nuc_PlaB1="path/to/polyAtails.tsv"
                    ,Cyt_PlaB1="path/to/polyAtails.tsv"
                    ,Poly_PlaB1="path/to/polyAtails.tsv"
                    ,Chr_PlaB2="path/to/polyAtails.tsv"
                    ,Chr_PlaB2="path/to/polyAtails.tsv"
                    ,Nuc_PlaB2="path/to/polyAtails.tsv"
                    ,Cyt_PlaB2="path/to/polyAtails.tsv"
                    ,Poly_PlaB2="path/to/polyAtails.tsv"
                    )

## HAR
tailsPaths_HAR <- c(Chr_HAR1="path/to/polyAtails.tsv"
                   ,Nuc_HAR1="path/to/polyAtails.tsv"
                   ,Cyt_HAR1="path/to/polyAtails.tsv"
                   ,Chr_HAR2="path/to/polyAtails.tsv"
                   ,Nuc_HAR2="path/to/polyAtails.tsv"
                   ,Cyt_HAR2="path/to/polyAtails.tsv"
                   )

## LEP
tailsPaths_LEP <- c(Chr_LEP1="path/to/polyAtails.tsv"
                   ,Nuc_LEP1="path/to/polyAtails.tsv"
                   ,Cyt_LEP1="path/to/polyAtails.tsv"
                   ,Poly_LEP1="path/to/polyAtails.tsv"
                   ,Chr_LEP2="path/to/polyAtails.tsv"
                   ,Chr_LEP2="path/to/polyAtails.tsv"
                   ,Nuc_LEP2="path/to/polyAtails.tsv"
                   ,Cyt_LEP2="path/to/polyAtails.tsv"
                   ,Poly_LEP2="path/to/polyAtails.tsv"
                   )

### All models
## Untreated
inferedRatesUntreated1_yesChpNpP_multi <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreated1_yesChpNpP_multi.rds")
inferedRatesUntreated1_yesChpNpP_single <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreated1_yesChpNpP_single.rds")
inferedRatesUntreated2_yesChpNpP_multi <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreated2_yesChpNpP_multi.rds")
inferedRatesUntreated2_yesChpNpP_single <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreated2_yesChpNpP_single.rds")
inferedRatesUntreatedMerged_yesChpNpP_multi <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreatedMerged_yesChpNpP_multi.rds")
inferedRatesUntreatedMerged_yesChpNpP_single <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreatedMerged_yesChpNpP_single.rds")

inferedRatesUntreatedMerged_noPoly <- readRDS("/path/to/Results/Untreated/thirdRun_sameGenesAllModels/inferedRatesUntreatedMerged_yesChpNpNoP.rds")

inferedRatesUntreatedMerged_noChpP <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_noChpP.rds")
inferedRatesUntreatedMerged_yesChpNoNpP <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_yesChpNoNpP.rds")
inferedRatesUntreatedMerged_yesChpNpNoP <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_yesChpNpNoP.rds")
inferedRatesUntreatedMerged_yesChpPNoNp <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_yesChpPNoNp.rds")
inferedRatesUntreatedMerged_yesPNoChp <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_yesPNoChp.rds")

inferedRatesUntreatedMerged_nucleoplasmicPrematureExport <- readRDS("/path/to/Results/Untreated/fiftRun_NucleoplasmicPrematureExport/inferedRatesUntreatedMerged_yesChpNpP.rds")

## PlaB
inferedRatesPlaB1_yesChpNpP_multi <- readRDS("/path/to/Results/PladienolideB/firstRun_FullModel/inferedRatesPlaB1_yesChpNpP.rds")
inferedRatesPlaB2_yesChpNpP_multi <- readRDS("/path/to/Results/PladienolideB/firstRun_FullModel/inferedRatesPlaB2_yesChpNpP.rds")
inferedRatesPlaBMerged_yesChpNpP_multi <- readRDS("/path/to/Results/PladienolideB/firstRun_FullModel/inferedRatesPlaBMerged_yesChpNpP.rds")

inferedRatesPlaBMerged_nucleoplasmicPrematureExport <- readRDS("/path/to/Results/PladienolideB/secondRun_NucleoplasmicPrematureExport/inferedRatesPlaBMerged_yesChpNpP.rds")

## HAR
inferedRatesHAR1_yesChpNpNoP_multi <- readRDS("/path/to/Results/Harringtonin/firstRun_FullModel/inferedRatesHAR1_yesChpNpNoP_multi.rds")
inferedRatesHAR2_yesChpNpNoP_multi <- readRDS("/path/to/Results/Harringtonin/firstRun_FullModel/inferedRatesHAR2_yesChpNpNoP_multi.rds")
inferedRatesHARMerged_yesChpNpNoP_multi <- readRDS("/path/to/Results/Harringtonin/firstRun_FullModel/inferedRatesHARMerged_yesChpNpNoP_multi.rds")

## LepB
inferedRatesLEP1_yesChpNpP_multi <- readRDS("/path/to/Results/LeptomycinB/firstRun_FullModel/inferedRatesLEP1_yesChpNpP_multi.rds")
inferedRatesLEP2_yesChpNpP_multi <- readRDS("/path/to/Results/LeptomycinB/firstRun_FullModel/inferedRatesLEP2_yesChpNpP_multi.rds")
inferedRatesLEPMerged_yesChpNpP_multi <- readRDS("/path/to/Results/LeptomycinB/firstRun_FullModel/inferedRatesLEPMerged_yesChpNpP_multi.rds")

### Replicates
## Expression levels - Figure 2A, Figures S7, S17, S22, S25
replicatesDataCorrelationsUntreated1 <- unlist(replicatesDataCorrelations(object1=inferedRatesUntreated1_yesChpNpP_multi,object2=inferedRatesUntreated2_yesChpNpP_multi,width=7,height=12,name="Untreated",lowSat=0.05,upSat=0.95))
round(c(min(replicatesDataCorrelationsUntreated1),max(replicatesDataCorrelationsUntreated1),median(replicatesDataCorrelationsUntreated1)),2)

replicatesDataCorrelationsPlaB1 <- unlist(replicatesDataCorrelations(object1=inferedRatesPlaB1_yesChpNpP_multi,object2=inferedRatesPlaB2_yesChpNpP_multi,width=7,height=12,name="PlaB",lowSat=0.05,upSat=0.95))
round(c(min(replicatesDataCorrelationsPlaB1),max(replicatesDataCorrelationsPlaB1),median(replicatesDataCorrelationsPlaB1)),2)

replicatesDataCorrelationsHAR1 <- unlist(replicatesDataCorrelations(object1=inferedRatesHAR1_yesChpNpNoP_multi,object2=inferedRatesHAR2_yesChpNpNoP_multi,width=7,height=12,name="HAR",lowSat=0.05,upSat=0.95))
round(c(min(replicatesDataCorrelationsHAR1),max(replicatesDataCorrelationsHAR1),median(replicatesDataCorrelationsHAR1)),2)

replicatesDataCorrelationsLEP1 <- unlist(replicatesDataCorrelations(object1=inferedRatesLEP1_yesChpNpP_multi,object2=inferedRatesLEP2_yesChpNpP_multi,width=7,height=12,name="LEP",lowSat=0.05,upSat=0.95))
round(c(min(replicatesDataCorrelationsLEP1),max(replicatesDataCorrelationsLEP1),median(replicatesDataCorrelationsLEP1)),2)

# Barplot
file1 <- inferedRatesUntreated1_yesChpNpP_multi$expressionData
file2 <- inferedRatesUntreated2_yesChpNpP_multi$expressionData

file1 <- file1[,-grep("0$",colnames(i))]
file2 <- file2[,-grep("0$",colnames(i))]

common_genes <- intersect(rownames(files1),rownames(files2))
cor <- diag(round(cor(files1[common_genes,],files2[common_genes,],method="s",use="c"),2))
times <- c(rep("Nascent",6),rep("Pre-existing",6),rep("Total",6))
Species <- rep(c("Chp","Chm","Np","Nm","C","P"),3)
ord <- unlist(lapply(c("chp","chm","np","nm","cy","py"),function(i){lapply(c("n0.33","p0.33","t"),function(j){paste0(i,j)})}))
cor <- cor[match(names(cor),ord)]
df <- data.frame(cbind("times" <- as.factor(times),"Species" =Species,"cor"=unname(cor)))
df $Species <- factor(df$Species,levels=c("Chp","Chm","Np","Nm","C","P"))
df$cor <- as.numeric(df$cor)

pl <- ggplot(df,aes(x =times , y = cor,fill =Species )) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=myPal[3:8]) +
  theme(legend.position="top") + theme( axis.title = element_text(size = 12), axis.text.x = element_text(size = 10),axis.title.x=element_blank()) + 
  ylab("Spearman Correlation") + 
  labs(text= element_text("Labeling pulses")) 

## Kinetic rates - Figure 2B, Figures S8, S18, S23, S26
replicatesRatesCorrelationsUntreated1 <- unlist(replicatesRatesCorrelations(object1=inferedRatesUntreated1_yesChpNpP_multi,object2=inferedRatesUntreated2_yesChpNpP_multi,width=7,height=12,name="Untreated",lowSat=0,upSat=1))
round(c(min(replicatesRatesCorrelationsUntreated1),max(replicatesRatesCorrelationsUntreated1),median(replicatesRatesCorrelationsUntreated1)),2)

replicatesRatesCorrelationsPlaB1 <- unlist(replicatesRatesCorrelations(object1=inferedRatesPlaB1_yesChpNpP_multi,object2=inferedRatesPlaB2_yesChpNpP_multi,width=7,height=12,name="PlaB",lowSat=0,upSat=1))
round(c(min(replicatesRatesCorrelationsPlaB1),max(replicatesRatesCorrelationsPlaB1),median(replicatesRatesCorrelationsPlaB1)),2)

replicatesRatesCorrelationsHAR1 <- unlist(replicatesRatesCorrelations(object1=inferedRatesHAR1_yesChpNpNoP_multi,object2=inferedRatesHAR2_yesChpNpNoP_multi,width=7,height=12,name="HAR",lowSat=0,upSat=1))
round(c(min(replicatesRatesCorrelationsHAR1),max(replicatesRatesCorrelationsHAR1),median(replicatesRatesCorrelationsHAR1)),2)

replicatesRatesCorrelationsLEP1 <- unlist(replicatesRatesCorrelations(object1=inferedRatesLEP1_yesChpNpP_multi,object2=inferedRatesLEP2_yesChpNpP_multi,width=7,height=12,name="LEP",lowSat=0,upSat=1))
round(c(min(replicatesRatesCorrelationsLEP1),max(replicatesRatesCorrelationsLEP1),median(replicatesRatesCorrelationsLEP1)),2)

# Barplot
common_genes <- intersect(rownames(inferedRatesUntreated1_yesChpNpP_multi$inferedRates),rownames(inferedRatesUntreated2_yesChpNpP_multi$inferedRates))
cor <- data.frame(unname(diag(round(cor(inferedRatesUntreated1_yesChpNpP_multi$inferedRates[common_genes,],inferedRatesUntreated2_yesChpNpP_multi$inferedRates[common_genes,],method="s",use="c"),2))))
cor <- cbind(c("k1","k2","k3","k4","k5","k6","k8","k7","k9"),cor)
colnames(cor) <- c("Rates","Cor")
myPalette <- brewer.pal(9,"Set3")
cor$Rates <- factor(cor$Rates,levels=c("k1","k2","k3","k4","k5","k6","k7","k8","k9"))

pl <- ggplot(cor,aes(x=Rates,y=Cor, fill=Rates)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=myPalette[1:9]) +
  theme( axis.title = element_text(size = 12), axis.text.x = element_text(size = 10), legend.position = 'none') + 
  ylab("Spearman Correlation") 

### Goodness of fit - Figure 2C, Figures S4, S9, S19
goodnessOfFitUntreatedMergedNoPoly <- unlist(GoodnessOfFit(inferedRates=inferedRatesUntreatedMerged_noPoly,width=7,height=12,name="Untreated_NoPoly",lowSat=0.025,upSat=0.975))
round(c(min(goodnessOfFitUntreatedMergedNoPoly),max(goodnessOfFitUntreatedMergedNoPoly),median(goodnessOfFitUntreatedMergedNoPoly)),2)

goodnessOfFitUntreatedMerged <- unlist(GoodnessOfFit(inferedRates=inferedRatesUntreatedMerged_yesChpNpP_multi,width=7,height=12,name="Untreated",lowSat=0.025,upSat=0.975))
round(c(min(goodnessOfFitUntreatedMerged),max(goodnessOfFitUntreatedMerged),median(goodnessOfFitUntreatedMerged)),2)

goodnessOfFitPlaBMerged <- unlist(GoodnessOfFit(inferedRates=inferedRatesPlaBMerged_yesChpNpP_multi,width=7,height=12,name="PlaB",lowSat=0.025,upSat=0.975))
round(c(min(goodnessOfFitPlaBMerged),max(goodnessOfFitPlaBMerged),median(goodnessOfFitPlaBMerged)),2)

# Barplot
expected <- inferedRatesUntreatedMerged_yesChpNpP_multi$expressionData[,-grep("0$",colnames(inferedRatesUntreatedMerged_yesChpNpP_multi$expressionData))]
inferred <- inferedRatesUntreatedMerged_yesChpNpP_multi$inferedData[,-grep("0$",colnames(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedData))]
cor <- diag(round(cor(expected,inferred,method="s",use="c"),2))
times <- c(rep("Nascent",6),rep("Pre-existing",6),rep("Total",6))
Species <- rep(c("Chp","Chm","Np","Nm","C","P"),3)
ord <- unlist(lapply(c("chp","chm","np","nm","cy","py"),function(i){lapply(c("n0.33","p0.33","t"),function(j){paste0(i,j)})}))
cor <- cor[match(names(cor),ord)]
df <- data.frame(cbind("times" <- as.factor(times),"Species" =Species,"cor"=unname(cor)))
df$Species <- factor(df$Species,levels=c("Chp","Chm","Np","Nm","C","P"))
df$cor <- as.numeric(df$cor)

pl <- ggplot(df,aes(x =times , y = cor,fill =Species )) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=myPal[3:8]) +
  theme(legend.position="top") + theme( axis.title = element_text(size = 12), axis.text.x = element_text(size = 10),axis.title.x=element_blank()) + 
  ylab("Spearman Correlation") + 
  labs(text= element_text("Labeling pulses")) 

### Distributions - Figure 2D, Figure S12 B
RatesDistributions(object=inferedRatesUntreatedMerged_yesChpNpP_multi,width=24,height=4,obj_name="Untreated",xlimTmp=c(-2,6)) #FIG. 2D
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
ks.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k1"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k1"])$p.value # 4.209688e-11
wilcox.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k1"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k1"],alternative="g")$p.value # 1.900747e-16

## Co-transcriptional splicing Untreated vs PlaB
ks.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k2"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k2"])$p.value # 0
wilcox.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k2"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k2"],alternative="g")$p.value # 3.895045e-140

## Mature chromatin RNA detachment Untreated vs PlaB
ks.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k3"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k3"])$p.value # 0
wilcox.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k3"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k3"],alternative="g")$p.value # 5.181007e-78

## Premature chromatin RNA detachment Untreated vs PlaB
ks.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k4"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k4"])$p.value # 0
wilcox.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k4"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k4"],alternative="l")$p.value # 1.676424e-05

## Post-transcriptional splicing Untreated vs PlaB
ks.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k5"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k5"])$p.value # 0
wilcox.test(inferedRatesUntreatedMerged_yesChpNpP_multi$inferedRates[,"k5"],inferedRatesPlaBMerged_yesChpNpP_multi$inferedRates[,"k5"],alternative="l")$p.value # 4.514181e-20

### Violins - Figures 3E, 4D, 5D
differentialDistributionPlot(object1=inferedRatesUntreatedMerged_yesChpNpP_multi,object2=inferedRatesPlaBMerged_yesChpNpP_multi,name="PlaBvsUntreated",width=6,height=4,breaks=100,col1="palevioletred",col2="lightblue",name1="Control",name2="PlaB",pchMed="-")
differentialDistributionPlot(object1=inferedRatesUntreatedMerged_noPoly,object2=inferedRatesHARMerged_yesChpNpNoP_multi,name="HarrvsUntreated",width=6,height=4,breaks=100,col1="palevioletred",col2="lightblue",name1="Control",name2="Harr",pchMed="-")
differentialDistributionPlot(object1=inferedRatesUntreatedMerged_yesChpNpP_multi,object2=inferedRatesLEPMerged_yesChpNpP_multi,name="LepBvsUntreated",width=6,height=4,breaks=100,col1="palevioletred",col2="lightblue",name1="Control",name2="LepB",pchMed="-")

### Heatmaps
## Untreated - Figure 2E, Figure S27A
HeatmapsPlotUntreated <- HeatmapsPlot(object=inferedRatesUntreatedMerged_yesChpNpP_multi
                                     ,obj_name="Untreated",n_clust=6,ratesWeight=3)

## PlaB - Figures 3D, 3F
HeatmapsPlotPlaB <- HeatmapsPlot(object=inferedRatesPlaBMerged_yesChpNpP_multi
                                ,obj_name="PlaB",n_clust=7,ratesWeight=3)

DifferentialHeatmapsPlotPlaB <- DifferentialHeatmapsPlot(refObject=inferedRatesUntreatedMerged_yesChpNpP_multi
                                                        ,object=inferedRatesPlaBMerged_yesChpNpP_multi
                                                        ,name="PlaBvsUntreated",n_clust=6,ylim=c(-2,2),ratesWeight=3)

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

## HAR - Figures 5C, 5E
HeatmapsPlotHAR <- HeatmapsPlot(object=inferedRatesHARMerged_yesChpNpNoP_multi
                               ,obj_name="HAR",n_clust=7,ratesWeight=3)

DifferentialHeatmapsPlotHAR <- DifferentialHeatmapsPlot(refObject=inferedRatesUntreatedMerged_noPoly
                                                       ,object=inferedRatesHARMerged_yesChpNpNoP_multi
                                                       ,name="HARvsUntreated",n_clust=6,ylim=c(-2,2),ratesWeight=3)

commonGenesTmp <- intersect(rownames(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates),rownames(inferedRatesUntreatedMerged_noPoly$inferedRates))
sum(log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k1"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k1"])>0)/length(commonGenesTmp)
# 0.4071429

# Genes switching from co- to post- transcriptional regulations
sum(log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k2"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k2"])>0
   &log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k3"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k3"])>0
   &log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k4"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k4"])<0
   &log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k5"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k5"])<0)/length(commonGenesTmp)

sum(log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k2"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k2"])<0
   &log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k3"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k3"])<0
   &log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k4"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k4"])>0
   &log2(inferedRatesHARMerged_yesChpNpNoP_multi$inferedRates[commonGenesTmp,"k5"]/inferedRatesUntreatedMerged_noPoly$inferedRates[commonGenesTmp,"k5"])>0)/length(commonGenesTmp)

## LEP - Figures 4A, 4E
HeatmapsPlotLEP <- HeatmapsPlot(object=inferedRatesLEPMerged_yesChpNpP_multi
                               ,obj_name="LEP",n_clust=8,ratesWeight=3)

DifferentialHeatmapsPlotLEP <- DifferentialHeatmapsPlot(refObject=inferedRatesUntreatedMerged_yesChpNpP_multi
													   ,object=inferedRatesLEPMerged_yesChpNpP_multi
													   ,name="LEPvsUntreated",n_clust=9,ylim=c(-2,2),ratesWeight=3)
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
refObject=inferedRatesUntreatedMerged_yesChpNpP_multi
object=inferedRatesLEPMerged_yesChpNpP_multi

ExampleData1 <- refObject$expressionData
polyFlag1 <- ("pyt"%in%colnames(ExampleData1))
if(polyFlag1)
{
ExampleData1 <- ExampleData1[,c("chpt","chmt","npt","nmt","cyt","pyt")]
}else{
ExampleData1 <- ExampleData1[,c("chpt","chmt","npt","nmt","cyt")]
}

inferedRates1 <- refObject$inferedRates

ExampleData2 <- object$expressionData
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

OurGenesHighDown <- names(which(inferedRates[,"k6"]<(-1.5)))
OurGenesHighUp <- names(which(inferedRates[,"k6"]>(1.5)))

OurGenesHighDownSymbol <- getBM(attributes = c("entrezgene_id","ensembl_gene_id","external_gene_name"),filter="entrezgene_id",values=OurGenesHighDown,mart = useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl")))
rownames(OurGenesHighDownSymbol) <- OurGenesHighDownSymbol[,1]
OurGenesHighDownSymbol <- dplyr::select(as.data.frame(OurGenesHighDownSymbol),columns='external_gene_name')

setdiff(OurGenesHighDown,rownames(OurGenesHighDownSymbol))
# "132241"

OurGenesHighDownSymbol <- rbind(OurGenesHighDownSymbol,data.frame(columns="RPL32P3")) # From genecards.org
rownames(OurGenesHighDownSymbol)[nrow(OurGenesHighDownSymbol)] <- "132241"

OurGenesHighUpSymbol <- getBM(attributes = c("entrezgene_id","ensembl_gene_id","external_gene_name"),filter="entrezgene_id",values=OurGenesHighUp,mart = useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl")))
rownames(OurGenesHighUpSymbol) <- OurGenesHighUpSymbol[,1]
OurGenesHighUpSymbol <- dplyr::select(as.data.frame(OurGenesHighUpSymbol),columns='external_gene_name')

setdiff(OurGenesHighUp,rownames(OurGenesHighUpSymbol))
# character(0)

OurGenesHighSymbol <- rbind(OurGenesHighDownSymbol,OurGenesHighUpSymbol)

functionTmp <- function(OurGenesHighSymbol,name,n_clust,show_rownames)
{
  inferedRatesTmp <- inferedRatesLEPMerged_yesChpNpP_multi
  inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp <- inferedRatesUntreatedMerged_yesChpNpP_multi

  inferedRatesTmp$inferedRates <- inferedRatesTmp$inferedRates[rownames(OurGenesHighSymbol),]
  inferedRatesTmp$inferedData <- inferedRatesTmp$inferedData[rownames(OurGenesHighSymbol),]
  inferedRatesTmp$expressionData <- inferedRatesTmp$expressionData[rownames(OurGenesHighSymbol),]

  rownames(inferedRatesTmp$inferedRates) <- OurGenesHighSymbol[rownames(inferedRatesTmp$inferedRates),1]
  rownames(inferedRatesTmp$inferedData) <- OurGenesHighSymbol[rownames(inferedRatesTmp$inferedData),1]
  rownames(inferedRatesTmp$expressionData) <- OurGenesHighSymbol[rownames(inferedRatesTmp$expressionData),1]
  
  inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedRates <- inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedRates[rownames(OurGenesHighSymbol),]
  inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedData <- inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedData[rownames(OurGenesHighSymbol),]
  inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$expressionData <- inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$expressionData[rownames(OurGenesHighSymbol),]
  
  rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedRates) <- OurGenesHighSymbol[rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedRates),1]
  rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedData) <- OurGenesHighSymbol[rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$inferedData),1]
  rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$expressionData) <- OurGenesHighSymbol[rownames(inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp$expressionData),1]
  
  DifferentialHeatmapsPlotLEPExtreme <- DifferentialHeatmapsPlot(refObject=inferedRatesUntreatedMerged_yesChpNpP_multi_Tmp
                                      ,object=inferedRatesTmp
                                      ,name=name,n_clust=n_clust,ylim=c(-2,2),ratesWeight=3,show_rownames=show_rownames)
}

# Heatmaps - Figure 4B
functionTmp(OurGenesHighDownSymbol,name="LEPvsUntreated_Extreme_Down",n_clust=2,show_rownames=TRUE)
functionTmp(OurGenesHighUpSymbol,name="LEPvsUntreated_Extreme_Up",n_clust=2,show_rownames=TRUE)

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
cgUtr3 <- cgContentFunction("/path/to/utr3.fa")
entropyUtr3 <- entropyFunction("/path/to/utr3.fa")

cgUtr5 <- cgContentFunction("/path/to/utr5.fa")
entropyUtr5 <- entropyFunction("/path/to/utr5.fa")

## Untreated - Figures S11, S14
GenesCharacterization(txdb=txdb,sortedRates=sortedRatesUntreated,obj_name="Untreated",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12,height=6)
GenesCharacterization(txdb,HeatmapsPlotUntreated,"Untreated",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12)

## Leptomycin B - Figure S24
GenesCharacterization(txdb,DifferentialHeatmapsPlotLEP,"LEP_DIFF",cgFeatures=list(cgUtr3,cgUtr5),entropyFeatures=list(entropyUtr3,entropyUtr5),width=12)

### PolyA tails characterization
## Untreated - Figure 2F, Figure S13
tailsUntreated <- TailsCharacterization(txdb=txdb,clustering=HeatmapsPlotUntreated,obj_name="Untreated"
                                       ,bam_untreated=NULL,bam_treated=bamPaths_Untreated,tailsPaths_untreatetd=NULL,tailsPaths_treated=tailsPaths_Untreated
                                       ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)
## Pladienolide B - Figure 3G, Figure S21 
tailsPlaB <- TailsCharacterization(txdb=txdb,clustering=HeatmapsPlotPlaB,obj_name="PlaB"
                                  ,bam_untreated=NULL,bam_treated=bamPaths_PlaB,tailsPaths_untreatetd=NULL,tailsPaths_treated=tailsPaths_PlaB
                                  ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)
tailsPlaBvsUntreated <- TailsCharacterization(txdb=txdb,clustering=DifferentialHeatmapsPlotPlaB,obj_name="PlaBvsUntreated"
                                      ,bam_untreated=bamPaths_Untreated,bam_treated=bamPaths_PlaB
                                      ,tailsPaths_untreatetd=tailsPaths_Untreated,tailsPaths_treated=tailsPaths_PlaB
                                      ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)

## Harringtonin - Figure 4F
tailsHAR <- TailsCharacterization(txdb=txdb,clustering=HeatmapsPlotHAR,obj_name="HAR"
                                ,bam_untreated=NULL,bam_treated=bamPaths_HAR,tailsPaths_untreatetd=NULL,tailsPaths_treated=tailsPaths_HAR
                                ,minoverlap_I=10,minoverlap_E=10,mergeReplicates=TRUE,width=12,height=6)

## Leptomycin B - Figure 5F 
tailsLEP <- TailsCharacterization(txdb=txdb,clustering=HeatmapsPlotLEP,obj_name="LEP"
                                ,bam_untreated=NULL,bam_treated=bamPaths_LEP,tailsPaths_untreatetd=NULL,tailsPaths_treated=tailsPaths_LEP
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
rownames(corMatrix) <- rownames(significanceMatrix) <- c("Pladienolide B","Leptomycin B","Harringtonin")

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

corMatrixFinal

### Polysomal RNA yield in Harringtonin - Figure 5A
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

### Premature Nucleoplasmic Export - Figure S32
DifferentialHeatmapsPlotLEP <- DifferentialHeatmapsPlot(refObject=inferedRatesUntreatedMerged_nucleoplasmicPrematureExport
                             ,object=inferedRatesPlaBMerged_nucleoplasmicPrematureExport
                             ,name="PlaBvsUntreated_NPE",n_clust=9,ylim=c(-2,2),ratesWeight=3)

### 5'TOP Factors - Figure S20
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
  untreated <- Untreated$inferedRates[which(lookupUntreated$hgnc_symbol==i),"k7"]
  plab <- PlaB$inferedRates[which(lookupPlaB$hgnc_symbol==i),"k7"]
  lepb <- LepB$inferedRates[which(lookupLepB$hgnc_symbol==i),"k7"]

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

pdf("RPS_RPL_Factors.pdf",width=8,height=8)
par(mfrow=c(2,2))
boxplot(PlFC,outline=FALSE,main="Log2FC Untreated-PlaB",xlab="PlaB")
boxplot(LepFC,outline=FALSE,main="Log2FC Untreated-LepB",xlab="LepB")
hist(PlFC,xlab="PlaB",main=NULL)
hist(LepFC,xlab="LepB",main=NULL)
dev.off()