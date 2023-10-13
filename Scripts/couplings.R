### Required libraries
library(combinat)
library(pheatmap)
library(igraph)

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

source("/path/to/allInternalFunctions.R")

### All inferred rates
## Simulated - Full Model
Simulated <- readRDS("/path/to/Results/dataSimulation/FullModel_Replicates/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")$exampleRates
## Simulated - Model without polysomal
Simulated_NoPoly <- readRDS("/path/to/Results/dataSimulation/noPolysomal/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k8_1.3.rds")$exampleRates
## Untreated - Full Model
Untreated <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreatedMerged_yesChpNpP_multi.rds")$inferedRates
## Untreated - Model without polysomal
Untreated_NoPoly <- readRDS("/path/to/Results/Untreated/fourthRun_NoPolysomal/inferedRatesUntreatedMerged_yesChpNpNoP_multi.rds")$inferedRates
## Untreated - Pladienolide B
PlaB <- readRDS("/path/to/Results/PladienolideB/firstRun_FullModel/inferedRatesPlaBMerged_yesChpNpP.rds")$inferedRates
## Untreated - Harringtonin
Harr <- readRDS("/path/to/Results/Harringtonin/firstRun_FullModel/inferedRatesHARMerged_yesChpNpNoP_multi.rds")$inferedRates
## Untreated - Leptomycin B
LepB <- readRDS("/path/to/Results/LeptomycinB/firstRun_FullModel/inferedRatesLEPMerged_yesChpNpP_multi.rds")$inferedRates

### First part of the analysis - correlative network at steady state
# Compute rates correlations and the corresponding significances
correlationSignificance <- function(matTmp,method="s")
{
    combi <- apply(combn(colnames(matTmp),2),2,paste,collapse='_')

    CorTmp <- cor(matTmp,method=method)
    RatesCorTmp <- CorTmp[t(upper.tri(CorTmp, diag = FALSE))]; names(RatesCorTmp) <- combi

    RatesCorSigTmp <- sapply(names(RatesCorTmp),function(i)
    {
        x <- strsplit(i,"_")[[1]]
        cor.test(matTmp[,x[[1]]],matTmp[,x[[2]]],method=method)$p.value
    })

    list(cor=RatesCorTmp,sig=RatesCorSigTmp)
}

# Full model-simulated
RatesCorSimulated <- correlationSignificance(Simulated)
RatesCorSimulatedSig <- RatesCorSimulated[[2]]
RatesCorSimulated <- RatesCorSimulated[[1]]
RatesCorSimulated[RatesCorSimulatedSig>=1e-5] <- 0  # Not-significant correlations are set to 0

# Full model-simulated without polysomal
RatesCorSimulated_NoPoly <- correlationSignificance(Simulated_NoPoly)
RatesCorSimulatedSig_NoPoly <- RatesCorSimulated_NoPoly[[2]]
RatesCorSimulated_NoPoly <- RatesCorSimulated_NoPoly[[1]]
RatesCorSimulated_NoPoly[RatesCorSimulatedSig_NoPoly>=1e-5] <- 0  # Not-significant correlations are set to 0

# Untreated
RatesCorUntreated <- correlationSignificance(Untreated)
RatesCorUntreatedSig <- RatesCorUntreated[[2]]
RatesCorUntreated <- RatesCorUntreated[[1]]
RatesCorUntreated[RatesCorUntreatedSig>=1e-5] <- 0  # Not-significant correlations are set to 0

# Untreated without polysomal
RatesCorUntreated_NoPoly <- correlationSignificance(Untreated_NoPoly)
RatesCorUntreatedSig_NoPoly <- RatesCorUntreated_NoPoly[[2]]
RatesCorUntreated_NoPoly <- RatesCorUntreated_NoPoly[[1]]
RatesCorUntreated_NoPoly[RatesCorUntreatedSig_NoPoly>=1e-5] <- 0  # Not-significant correlations are set to 0

# PlaB
RatesCorPlaB <- correlationSignificance(PlaB)
RatesCorPlaBSig <- RatesCorPlaB[[2]]
RatesCorPlaB <- RatesCorPlaB[[1]]
RatesCorPlaB[RatesCorPlaBSig>=1e-5] <- 0  # Not-significant correlations are set to 0

# LepB
RatesCorLepB <- correlationSignificance(LepB)
RatesCorLepBSig <- RatesCorLepB[[2]]
RatesCorLepB <- RatesCorLepB[[1]]
RatesCorLepB[RatesCorLepBSig>=1e-5] <- 0  # Not-significant correlations are set to 0

# Harringtonin
RatesCorHarr <- correlationSignificance(Harr)
RatesCorHarrSig <- RatesCorHarr[[2]]
RatesCorHarr <- RatesCorHarr[[1]]
RatesCorHarr[RatesCorHarrSig>=1e-5] <- 0  # Not-significant correlations are set to 0

# Rate-specific shuffling to identify random correlations
shufflingFunction <- function(matTmp,iterations=1000,seedTmp=0)
{
    outTmp2 <- t(sapply(1:iterations,function(j)
    {
        set.seed(j+seedTmp)
        mat <- sapply(1:ncol(matTmp),function(i){
            sampleTmp <- sample(1:nrow(matTmp),nrow(matTmp))
            matTmp[,i] <- matTmp[sampleTmp,i]
            })
        rownames(mat) <- c(1:nrow(matTmp))
        colnames(mat) <- colnames(matTmp)
        mat

        outTmp <- cor(mat,method="s")
        outTmp[upper.tri(outTmp,diag=FALSE)]
    }))

    colnames(outTmp2) <- apply(combn(colnames(matTmp),2),2,paste,collapse='_')
    outTmp2
}

NullSimulated <- shufflingFunction(Simulated)
NullUntreated <- shufflingFunction(Untreated)
NullSimulated_NoPoly <- shufflingFunction(Simulated_NoPoly)
NullUntreated_NoPoly <- shufflingFunction(Untreated_NoPoly)
NullPlaB <- shufflingFunction(PlaB)
NullHarr <- shufflingFunction(Harr)
NullLepB <- shufflingFunction(LepB)

### Plot - Figure 6A and Figure S27B
pdf("allGenesNetworks.pdf",width=10,height=6)
par(mfrow=c(2,4))
    # Plot of networks
    functionTmp <- function(RatesCor,Null,sfTh,mainTmp="",relevantRates=NULL,randomCorrelations=NULL)
    {
        # Rates re-name
        names(RatesCor) <- tryCatch(gsub("k10","k9",names(RatesCor)),error=function(e)names(RatesCor))
        names(RatesCor) <- tryCatch(gsub("k7","keightnew",names(RatesCor)),error=function(e)names(RatesCor))
        names(RatesCor) <- tryCatch(gsub("k8","ksevennew",names(RatesCor)),error=function(e)names(RatesCor))
        names(RatesCor) <- tryCatch(gsub("keightnew","k8",names(RatesCor)),error=function(e)names(RatesCor))
        names(RatesCor) <- tryCatch(gsub("ksevennew","k7",names(RatesCor)),error=function(e)names(RatesCor))

        colnames(Null) <- tryCatch(gsub("k10","k9",colnames(Null)),error=function(e)colnames(Null))
        colnames(Null) <- tryCatch(gsub("k7","keightnew",colnames(Null)),error=function(e)colnames(Null))
        colnames(Null) <- tryCatch(gsub("k8","ksevennew",colnames(Null)),error=function(e)colnames(Null))
        colnames(Null) <- tryCatch(gsub("keightnew","k8",colnames(Null)),error=function(e)colnames(Null))
        colnames(Null) <- tryCatch(gsub("ksevennew","k7",colnames(Null)),error=function(e)colnames(Null))
        
        # Random correlations re-name
        if(!is.null(randomCorrelations))
        {
            names(randomCorrelations) <- tryCatch(gsub("k10","k9",names(randomCorrelations)),error=function(e)names(randomCorrelations))
            names(randomCorrelations) <- tryCatch(gsub("k7","keightnew",names(randomCorrelations)),error=function(e)names(randomCorrelations))
            names(randomCorrelations) <- tryCatch(gsub("k8","ksevennew",names(randomCorrelations)),error=function(e)names(randomCorrelations))
            names(randomCorrelations) <- tryCatch(gsub("keightnew","k8",names(randomCorrelations)),error=function(e)names(randomCorrelations))
            names(randomCorrelations) <- tryCatch(gsub("ksevennew","k7",names(randomCorrelations)),error=function(e)names(randomCorrelations))
        }

        # Significant edges
        foe <- sapply(names(RatesCor),function(i)
        {
            x <- tryCatch(sum(Null[,i]<=RatesCor[[i]])/nrow(Null),error=function(e)0.5);x<(sfTh/2) | x>(1-sfTh/2)
        })
        foe <- RatesCor[foe]
        
        if(!is.null(randomCorrelations))
        {
            namesTmp <- intersect(names(randomCorrelations),names(foe))
            namesTmp <- namesTmp[sign(randomCorrelations[namesTmp])==sign(foe[namesTmp])]
            foe <- foe[setdiff(names(foe),namesTmp)]
        }

        if(!is.null(relevantRates))
        {
            foe <- foe[names(foe)%in%relevantRates]
        }

        # Plot
        qwe <- graph.data.frame(data.frame("from"=sapply(strsplit(names(foe),"_"),"[[",1)
                                          ,"to"=sapply(strsplit(names(foe),"_"),"[[",2))
                                ,directed=FALSE,vertices=paste0("k",1:9))
        igraph.options(plot.layout=layout.circle)
        colTmp <- c("blue","white","red")[sign(foe)+2]
        plot.igraph(qwe,edge.width=abs(foe)*5,edge.color=colTmp,vertex.size=35,main=mainTmp)
        foe
    }

    # Full model-simulated
    significantCorSimulated <- functionTmp(RatesCorSimulated,NullSimulated,sfTh=0.001,mainTmp="Simulated")
    # Untreated
    significantCorUntreated <- functionTmp(RatesCor=RatesCorUntreated,Null=NullUntreated,sfTh=0.001,mainTmp="Untreated")
    # PlaB
    significantCorPlaB <- functionTmp(RatesCor=RatesCorPlaB,Null=NullPlaB,sfTh=0.001,mainTmp="Pladienolide B")
    # LepB
    significantCorLepB <- functionTmp(RatesCor=RatesCorLepB,Null=NullLepB,sfTh=0.001,mainTmp="Leptomycin B")

    # Simulated without polysomal
    significantCorSimulated_NoPoly <- functionTmp(RatesCor=RatesCorSimulated_NoPoly,Null=NullSimulated_NoPoly,sfTh=0.001,mainTmp="Simulated")
    # Full model-simulated without polysomal
    significantCorUntreated_NoPoly <- functionTmp(RatesCor=RatesCorUntreated_NoPoly,Null=NullUntreated_NoPoly,sfTh=0.001,mainTmp="Untreated")
    # Harringtonin
    significantCorHarr <- functionTmp(RatesCor=RatesCorHarr,Null=NullHarr,sfTh=0.001,mainTmp="Harringtonine")

    par(mfrow=c(2,3))
    significantCleanedCorUntreated <- functionTmp(RatesCor=RatesCorUntreated,Null=NullUntreated,sfTh=0.001,mainTmp="Untreated",relevantRates=NULL,randomCorrelations=significantCorSimulated)
    significantCleanedCorPlaB <- functionTmp(RatesCor=RatesCorPlaB,Null=NullPlaB,sfTh=0.001,mainTmp="Pladienolide B",relevantRates=NULL,randomCorrelations=significantCorSimulated)
    significantCleanedCorLepB <- functionTmp(RatesCor=RatesCorLepB,Null=NullLepB,sfTh=0.001,mainTmp="Leptomycin B",relevantRates=NULL,randomCorrelations=significantCorSimulated)

    significantCleanedCorUntreated_NoPoly <- functionTmp(RatesCor=RatesCorUntreated_NoPoly,Null=NullUntreated_NoPoly,sfTh=0.001,mainTmp="Untreated",relevantRates=NULL,randomCorrelations=significantCorSimulated_NoPoly)
    significantCleanedCorHarr <- functionTmp(RatesCor=RatesCorHarr,Null=NullHarr,sfTh=0.001,mainTmp="Harringtonine",relevantRates=NULL,randomCorrelations=significantCorSimulated_NoPoly)
dev.off()

### Second part of the analysis - network refinement with perturbation data
### Figure 6B and Figures S28 A-D
GenesCorr <- function(mat1,mat2,ratioTh,expectedModel,freqFilt=0,iterations=1000,seedTmp=0,sfTh=0.01,redifineExpectedModel=FALSE)
{
    # Expected rates from the correlative analysis
    expectedModel <- sign(expectedModel)
    common_genes <- intersect(rownames(mat1),rownames(mat2))
    mat1 <- mat1[common_genes,]
    mat2 <- mat2[common_genes,]
    n_genes <- nrow(mat1)

    mat1 <- mat1[,intersect(colnames(mat2),colnames(mat1))]

    # Rates re-name
    colnames(mat1) <- tryCatch(gsub("k10","k9",colnames(mat1)),error=function(e)colnames(mat1))
    colnames(mat2) <- tryCatch(gsub("k10","k9",colnames(mat2)),error=function(e)colnames(mat2))
    colnames(mat1) <- tryCatch(gsub("k7","keightnew",colnames(mat1)),error=function(e)colnames(mat1))
    colnames(mat2) <- tryCatch(gsub("k7","keightnew",colnames(mat2)),error=function(e)colnames(mat2))
    colnames(mat1) <- tryCatch(gsub("k8","ksevennew",colnames(mat1)),error=function(e)colnames(mat1))
    colnames(mat2) <- tryCatch(gsub("k8","ksevennew",colnames(mat2)),error=function(e)colnames(mat2))
    colnames(mat1) <- tryCatch(gsub("keightnew","k8",colnames(mat1)),error=function(e)colnames(mat1))
    colnames(mat2) <- tryCatch(gsub("keightnew","k8",colnames(mat2)),error=function(e)colnames(mat2))
    colnames(mat1) <- tryCatch(gsub("ksevennew","k7",colnames(mat1)),error=function(e)colnames(mat1))
    colnames(mat2) <- tryCatch(gsub("ksevennew","k7",colnames(mat2)),error=function(e)colnames(mat2))

    # Number of genes supporting a specific modulation
    functionTmp <- function(mat1,mat2,ratioTh)
    {
        ratios <- mat2/mat1
        ratios[ratios<(1+ratioTh)&ratios>(1-ratioTh)]=0 
    
        ratios[ratios>(1+ratioTh)&ratios!=0]=1
        ratios[ratios<(1-ratioTh)&ratios!=0]=(-1)
    
        combos <- combn(colnames(mat2),2)
        Freq <- sapply(1:ncol(combos),function(i)
        {
            rbind(length(which(ratios[,combos[1,i]]*ratios[,combos[2,i]]>0)),
                  length(which(ratios[,combos[1,i]]*ratios[,combos[2,i]]<0)),
                  length(which(ratios[,combos[1,i]]*ratios[,combos[2,i]]==0)) - length(which(ratios[,combos[1,i]]==0&ratios[,combos[2,i]]==0)),
                  length(which(ratios[,combos[1,i]]==0&ratios[,combos[2,i]]==0)))
            })
        rownames(Freq) <- c("Positive","Negative","NotSupportive","Zeros")
        colnames(Freq) <- apply(combn(colnames(mat1),2),2,paste,collapse="_")
        Freq
    }

    # Real data
    Freq <- functionTmp(mat1,mat2,ratioTh)
    
    # Random data
    RandFreq <- lapply(1:iterations,function(j)
    {
        set.seed(j+seedTmp)

        mat1Tmp <- sapply(1:ncol(mat1),function(i){
            sampleTmp <- sample(1:nrow(mat1),nrow(mat1))
            mat1[sampleTmp,i]
            })
        rownames(mat1Tmp) <- c(1:nrow(mat1))
        colnames(mat1Tmp) <- colnames(mat1)

        mat2Tmp <- sapply(1:ncol(mat2),function(i){
            sampleTmp <- sample(1:nrow(mat2),nrow(mat2))
            mat2[sampleTmp,i]
            })
        rownames(mat2Tmp) <- c(1:nrow(mat2))
        colnames(mat2Tmp) <- colnames(mat2)

        functionTmp(mat1Tmp,mat2Tmp,ratioTh)
    })

    if(redifineExpectedModel)
    {
        ### If we want to define the model according to the specific response
        ### we re-define the expectedModel object (keeping the edges but changing the sing)
        expectedModelOld <- expectedModel
        expectedModel[names(expectedModel)] <- c(-1,1)[(Freq[,names(expectedModel)]["Positive",]>Freq[,names(expectedModel)]["Negative",])+1]

        switchedEdges <- names(which(expectedModel!=expectedModelOld))

    }

    estimatedFreq <- sapply(names(expectedModel),function(i)Freq[min(((-1*sign(expectedModel[i]))+2),2),i]/sum(Freq[1:2,i]))
    estimatedFreqRand <- t(sapply(RandFreq,function(j)sapply(names(expectedModel),function(i)j[min(((-1*sign(expectedModel[i]))+2),2),i]/sum(j[1:2,i]))))

    notSupportiveFreq <- sapply(names(expectedModel),function(i)Freq[3,i]/sum(Freq[,i]))
    notSupportiveFreqRand <- t(sapply(RandFreq,function(j)sapply(names(expectedModel),function(i)j[3,i]/sum(j[,i]))))

    boolFreqRand <- sapply(names(estimatedFreq),function(i)(sum(estimatedFreqRand[,i]>estimatedFreq[[i]])/iterations)<sfTh)
    boolNotSuppFreqRand <- sapply(names(notSupportiveFreq),function(i)(sum(notSupportiveFreqRand[,i]>notSupportiveFreq[[i]])/iterations)<sfTh)

    # Median values for the random configurations
    medianFreqRand <- sapply(colnames(RandFreq[[1]]),function(i)apply(sapply(RandFreq,function(j)j[,i]),1,median))

    # Random bars
    barplot(medianFreqRand[,names(expectedModel)],las=2,col=c("red","blue","green","white"),main="Median random")
    
    # Real data bars
    barplot(Freq[,names(expectedModel)],las=2,col=c("red","blue","green","white"),main="")

    # Genes supporting the external model/the most frequent model
    colTmp <- c("blue",NA,"red")[sign(expectedModel)+2]
    names(colTmp) <- names(expectedModel)
    colTmp[names(colTmp)%in%switchedEdges] <- "gray"
    barplot(estimatedFreq,las=2,ylim=c(0,1),ylab="Supportive Genes [%]",col=colTmp,main="All Candidate Edges")
    
    # Final edges
    colTmp <- c("blue",NA,"red")[sign(expectedModel[boolFreqRand])+2]
    names(colTmp) <- names(expectedModel[boolFreqRand])
    colTmp[names(colTmp)%in%switchedEdges] <- "gray"
    barplot(estimatedFreq[boolFreqRand],las=2,ylim=c(0,1),ylab="Supportive Genes [%]",col=colTmp,main="Bootstrap Supported Edges")
    
    out1 <- estimatedFreq[boolFreqRand]
    out2 <- c("blue",NA,"red")[sign(expectedModel[boolFreqRand])+2]

    estimatedFreq <- estimatedFreq[estimatedFreq>=freqFilt&boolFreqRand]

    genesTmp <- lapply(names(estimatedFreq),function(conditionTmp)
    {
        ratios <- mat2/mat1
        ratios[ratios<(1+ratioTh)&ratios>(1-ratioTh)]=0 
    
        ratios[ratios>(1+ratioTh)&ratios!=0]=1
        ratios[ratios<(1-ratioTh)&ratios!=0]=(-1)

        i <- strsplit(conditionTmp,"_")[[1]][[1]]
        j <- strsplit(conditionTmp,"_")[[1]][[2]]
        if(conditionTmp%in%names(expectedModel))
        {
            if(expectedModel[conditionTmp]==1)
            {
                common_genes[which(ratios[,i]*ratios[,j]>0)]
            }else{
                common_genes[which(ratios[,i]*ratios[,j]<0)]
            }               
        }
    })
    names(genesTmp) <- names(estimatedFreq)
    
    list(Freq,common_genes,genesTmp,out1,out2)
}

Untreated_PlaB_0.25 <- GenesCorr(mat1=Untreated,mat2=PlaB,ratioTh=0.25,expectedModel=significantCleanedCorUntreated,freqFilt=0,iterations=1000,seedTmp=0,sfTh=0.001,redifineExpectedModel=FALSE);dev.off()
Untreated_LepB_0.25 <- GenesCorr(mat1=Untreated,mat2=LepB,ratioTh=0.25,expectedModel=significantCleanedCorUntreated,freqFilt=0,iterations=1000,seedTmp=0,sfTh=0.001,redifineExpectedModel=FALSE);dev.off()
Untreated_Harr_0.25 <- GenesCorr(mat1=Untreated_NoPoly,mat2=Harr,ratioTh=0.25,expectedModel=significantCleanedCorUntreated_NoPoly,freqFilt=0,iterations=1000,seedTmp=0,sfTh=0.001,redifineExpectedModel=FALSE);dev.off()

pdf("ratesDeltas_0.25.pdf",width=18,height=8)
par(mfrow=c(3,4))
Untreated_PlaB_0.25_ReDefined <- GenesCorr(mat1=Untreated,mat2=PlaB,ratioTh=0.25,expectedModel=significantCleanedCorUntreated,freqFilt=0,iterations=1000,seedTmp=0,sfTh=0.001,redifineExpectedModel=TRUE)
Untreated_LepB_0.25_ReDefined <- GenesCorr(mat1=Untreated,mat2=LepB,ratioTh=0.25,expectedModel=significantCleanedCorUntreated,freqFilt=0,iterations=1000,seedTmp=0,sfTh=0.001,redifineExpectedModel=TRUE)
Untreated_Harr_0.25_ReDefined <- GenesCorr(mat1=Untreated_NoPoly,mat2=Harr,ratioTh=0.25,expectedModel=significantCleanedCorUntreated_NoPoly,freqFilt=0,iterations=1000,seedTmp=0,sfTh=0.001,redifineExpectedModel=TRUE)
dev.off()

# Set switched edges to gray
Untreated_PlaB_0.25_Final <- Untreated_PlaB_0.25_ReDefined[4:5]
Untreated_LepB_0.25_Final <- Untreated_LepB_0.25_ReDefined[4:5]
Untreated_Harr_0.25_Final <- Untreated_Harr_0.25_ReDefined[4:5]

Untreated_PlaB_0.25_Final[[2]][!names(Untreated_PlaB_0.25_Final[[1]])%in%names(Untreated_PlaB_0.25[[4]])] <- "gray"
Untreated_LepB_0.25_Final[[2]][!names(Untreated_LepB_0.25_Final[[1]])%in%names(Untreated_LepB_0.25[[4]])] <- "gray"
Untreated_Harr_0.25_Final[[2]][!names(Untreated_Harr_0.25_Final[[1]])%in%names(Untreated_Harr_0.25[[4]])] <- "gray"

# Final networks - Figure 6 C-D and Figures S28 E
pdf("modulationsNetworks_0.25.pdf",width=10,height=4)
    par(mfrow=c(1,4))
    functionTmp <- function(EdgesWidth,EdgesCol,mainTmp="")
    {
        foe <- EdgesWidth
        qwe <- graph.data.frame(data.frame("from"=sapply(strsplit(names(foe),"_"),"[[",1)
                                          ,"to"=sapply(strsplit(names(foe),"_"),"[[",2))
                                ,directed=FALSE,vertices=paste0("k",1:9))
        igraph.options(plot.layout=layout.circle)
        colTmp <- EdgesCol
        plot.igraph(qwe,edge.width=abs(foe)*5,edge.color=colTmp,vertex.size=35,main=mainTmp)
        foe
    }

    significantCleanedCorPlaB <- functionTmp(Untreated_PlaB_0.25_Final[[1]],Untreated_PlaB_0.25_Final[[2]],mainTmp="Pladienolide B")
    significantCleanedCorLepB <- functionTmp(Untreated_LepB_0.25_Final[[1]],Untreated_LepB_0.25_Final[[2]],mainTmp="Leptomycin B")

    names(Untreated_PlaB_0.25_Final[[2]]) <- names(Untreated_PlaB_0.25_Final[[1]])
    names(Untreated_LepB_0.25_Final[[2]]) <- names(Untreated_LepB_0.25_Final[[1]])
    
    commonEdges <- intersect(names(Untreated_PlaB_0.25_Final[[1]]),names(Untreated_LepB_0.25_Final[[1]]))
    commonWeights <- (Untreated_PlaB_0.25_Final[[1]][commonEdges]+Untreated_LepB_0.25_Final[[1]][commonEdges])/2
    
    commonColors <- sapply(commonEdges,function(i)
    {
        if(Untreated_PlaB_0.25_Final[[2]][[i]]!=Untreated_LepB_0.25_Final[[2]][[i]])
        {
            "white"
        }else{Untreated_PlaB_0.25_Final[[2]][[i]]}
    })

    commonWeights <- commonWeights[names(commonColors[commonColors!="white"])]
    commonEdges <- names(commonWeights)

    significantCleanedCorInt <- functionTmp(commonWeights,commonColors[commonEdges],mainTmp="Intersection")

    significantCleanedCorHarr <- functionTmp(Untreated_Harr_0.25_Final[[1]],Untreated_Harr_0.25_Final[[2]],mainTmp="Harringtonine")
dev.off()

### Higher orders - Figure 6E and Figures S29
Untreated_PlaB_0.25_N <- table(table(unlist(Untreated_PlaB_0.25_ReDefined[[3]][commonEdges])))
Untreated_LepB_0.25_N <- table(table(unlist(Untreated_LepB_0.25_ReDefined[[3]][commonEdges])))

x_PlaB <- table(unlist(Untreated_PlaB_0.25_ReDefined[[3]][commonEdges]))
y_PlaB <- split(names(x_PlaB),x_PlaB)

z_PlaB <- unlist(sapply(seq_along(Untreated_PlaB_0.25_ReDefined[[3]][commonEdges]),function(i)rep(names(Untreated_PlaB_0.25_ReDefined[[3]][commonEdges])[[i]],length(Untreated_PlaB_0.25_ReDefined[[3]][commonEdges][[i]]))))
names(z_PlaB) <- unlist(Untreated_PlaB_0.25_ReDefined[[3]][commonEdges])
z_PlaB <- split(z_PlaB,names(z_PlaB))

k_PlaB <- sapply(sapply(y_PlaB,function(i)(unname(unlist(z_PlaB[i])))),function(j)
{
    c("coPro"=sum(j %in% c("k1_k2","k1_k3","k2_k3"))
     ,"postPro"=sum(j %in% c("k4_k5"))
     ,"crossPro"=sum(j %in% c("k2_k4","k3_k4","k3_k5"))
     ,"decay"=sum(j %in% c("k6_k9","k8_k7","k8_k9"))
     ,"bridge"=sum(j %in% c("k4_k9","k5_k7")))/length(j)
})
k_PlaB <- sapply(colnames(k_PlaB),function(i)Untreated_PlaB_0.25_N[[i]]*k_PlaB[,i])

x_LepB <- table(unlist(Untreated_LepB_0.25_ReDefined[[3]][commonEdges]))
y_LepB <- split(names(x_LepB),x_LepB)

z_LepB <- unlist(sapply(seq_along(Untreated_LepB_0.25_ReDefined[[3]][commonEdges]),function(i)rep(names(Untreated_LepB_0.25_ReDefined[[3]][commonEdges])[[i]],length(Untreated_LepB_0.25_ReDefined[[3]][commonEdges][[i]]))))
names(z_LepB) <- unlist(Untreated_LepB_0.25_ReDefined[[3]][commonEdges])
z_LepB <- split(z_LepB,names(z_LepB))

k_LepB <- sapply(sapply(y_LepB,function(i)(unname(unlist(z_LepB[i])))),function(j)
{
    c("coPro"=sum(j %in% c("k1_k2","k1_k3","k2_k3"))
     ,"postPro"=sum(j %in% c("k4_k5"))
     ,"crossPro"=sum(j %in% c("k2_k4","k3_k4","k3_k5"))
     ,"decay"=sum(j %in% c("k6_k9","k8_k7","k8_k9"))
     ,"bridge"=sum(j %in% c("k4_k9","k5_k7")))/length(j)
})
k_LepB <- sapply(colnames(k_LepB),function(i)Untreated_LepB_0.25_N[[i]]*k_LepB[,i])

k_PlaB <- cbind(k_PlaB,rep(NaN,4),c(3,1,3,3,2)/12)
k_LepB <- cbind(k_LepB,rep(NaN,4),c(3,1,3,3,2)/12)

colsTmp <- c("k1_k2"="black"
            ,"k1_k3"="black"
            ,"k2_k3"="black"
            ,"k4_k5"="brown"
            ,"k2_k4"="purple"
            ,"k3_k4"="purple"
            ,"k3_k5"="purple"
            ,"k6_k9"="green"
            ,"k8_k7"="green"
            ,"k8_k9"="green"
            ,"k4_k9"="orange"
            ,"k5_k7"="orange")

pdf("finalEdgesCoRegulation.pdf",width=9,height=6)
par(mfrow=c(2,3))
barplot(Untreated_PlaB_0.25_N,main="Pladienolide B",xlab="Number of Edges",ylab="Number of Genes",las=2)
barplot(Untreated_LepB_0.25_N,main="Leptomycin B",xlab="Number of Edges",ylab="Number of Genes",las=2)
plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
barplot(apply(k_PlaB,2,function(i)i/sum(i)),main="PlaB",las=2,col=unique(colsTmp),xlab="Number of Edges",ylab="")
barplot(apply(k_LepB,2,function(i)i/sum(i)),main="LepB",las=2,col=unique(colsTmp),xlab="Number of Edges",ylab="")
functionTmp(sign(commonWeights),colsTmp[names(commonWeights)],mainTmp="Legend")
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
