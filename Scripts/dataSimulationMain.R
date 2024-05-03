### External arguments extracted from simulateGene.sh
args = commandArgs(trailingOnly=TRUE)
 
for(v in args)
{
       vTmp <- strsplit(v,"=")[[1]]
       assign(vTmp[[1]],vTmp[[2]])
}

reps <- (unlist(strsplit(reps,",")))
TauFractions <- (unlist(strsplit(TauFractions,",")))
TauPoly <- (unlist(strsplit(TauPoly,",")))
Flags <- unlist(strsplit(Flags,","))
rates <- (unlist(strsplit(rates,",")))

reps <- as.numeric(reps)
translationCoeff <- as.numeric(translationCoeff)
TauFractions <- (as.numeric(TauFractions))
TauPoly <- (as.numeric(TauPoly))
nGenes <- as.numeric(nGenes)
cpus <- as.numeric(cpus)
ZeroThresh <- as.numeric(ZeroThresh)
MultFact <- as.numeric(MultFact)

### Libraries
library(parallel)
library(deSolve)
library(pheatmap)

source("/path/to/allInternalFunctions.R")
# source("path/to/allInternalFunctionsPatch_nuclearPrematureDecay.R") # Execute this line to simulate data with Nucleoplasmic Premature RNA degradation.

### Average gene
INSPEcT_nascent <- read.table("/path/to/regulated_genes_features.xls",sep="\t",header=TRUE)

### Median kinetic rates
synthesis=signif(median(INSPEcT_nascent[,"synthesis_0"]),2)
processing=signif(median(INSPEcT_nascent[,"processing_0"]),2)
degradation=signif(median(INSPEcT_nascent[,"degradation_0"]),2)

### Mean rates
exampleRates <- c(k1=synthesis
                 ,k2=processing
                 ,k3=processing
                 ,k4=0.1*processing
                 ,k5=0.25*processing
                 ,k6=1.5*processing
                 ,k7=translationCoeff*degradation
                 ,k8=degradation
                 ,k9=0.5*degradation
                 ,k10=degradation
                 ,k11=0.5*degradation)

initialRates <- list(10**ceiling(log10(exampleRates))*(1+seq_along(exampleRates)*1e-5)
                    ,10**round(log10(exampleRates))*(1+seq_along(exampleRates)*1e-5)
                    ,10**floor(log10(exampleRates))*(1+seq_along(exampleRates)*1e-5))

exampleRates <- exampleRates[names(exampleRates)%in%rates]
initialRates <- lapply(initialRates,function(i)i[names(i)%in%rates])

### Variation Coefficients
CVspecies = readRDS("/path/to/Files/CVs.rds")

## Add CVs for additional time points (1h and 2h)
CVspecies = lapply(CVspecies,function(i)
{
  addTmp1 <- addTmp2 <- i[grep("0.33",names(i))]
  names(addTmp1) <- gsub("0.33","1",names(addTmp1))
  names(addTmp2) <- gsub("0.33","2",names(addTmp2))

  c(i,addTmp1,addTmp2)
})

### Simulation (loop)
## Number of replicates
for(rep in reps)
{
    print(rep)
    
    simulatedDataset <- simulateData(exampleRates=exampleRates # Mean rates to sample.
                                    ,TauFractions=TauFractions # Time points with cellular fractionation.
                                    ,TauPoly=TauPoly # Time points with polysomal profiling.
                                    ,TauTotal=NULL # Time points with total RNA profiling (FIXED to null for our purposes).
                                    ,noise=TRUE # TRUE to add noise to data.
                                    ,CVs=CVspecies # Variation Coefficient for noise.
                                    ,Reps=rep # Number of replicates.
                                    ,nGenes=nGenes # Number of genes to be simulated.
                                    ,seed=1 # Seed for reproducibility.
                                    ,ZeroThresh=ZeroThresh # Minimum expression value.
                                    ,MultFact=MultFact # Number of nGenes to be modeled.
                                    ,cpus=cpus) # Number of CPUs.

    if(nrow(simulatedDataset$exampleData)<nGenes)
    {
        print("Multiple iterations needed.")
        print(2*ceiling(nGenes/(nrow(simulatedDataset$exampleData)-1)))
        for(i in 1:(2*ceiling(nGenes/(nrow(simulatedDataset$exampleData)-1))))
        {
        simulatedDatasetTmp <- simulateData(exampleRates=exampleRates # Mean rates to sample.
                                           ,TauFractions=TauFractions # Time points with cellular fractionation.
                                           ,TauPoly=TauPoly # Time points with polysomal profiling.
                                           ,TauTotal=NULL # Time points with total RNA profiling (FIXED to null for our purposes).
                                           ,noise=TRUE # TRUE to add noise to data.
                                           ,CVs=CVspecies # Variation Coefficient for noise.
                                           ,Reps=rep # Number of replicates.
                                           ,nGenes=nGenes # Number of genes to be simulated.
                                           ,seed=1+i # Seed for reproducibility.
                                           ,ZeroThresh=ZeroThresh # Minimum expression value.
                                           ,MultFact=MultFact # Number of nGenes to be modeled.
                                           ,cpus=cpus) # Number of CPUs.      

        simulatedDataset$exampleData <- rbind(simulatedDataset$exampleData,simulatedDatasetTmp$exampleData[-1,])
        simulatedDataset$DevDataTmp <- rbind(simulatedDataset$DevDataTmp,simulatedDatasetTmp$DevDataTmp[-1,])
        simulatedDataset$exampleRates <- rbind(simulatedDataset$exampleRates,simulatedDatasetTmp$exampleRates[-1,])
        simulatedDataset$singleReplicateExpressionData <- lapply(1:rep,function(i){rbind(simulatedDataset$singleReplicateExpressionData[[i]],simulatedDatasetTmp$singleReplicateExpressionData[[i]][-1,])})

        rownames(simulatedDataset$exampleData) <- 1:nrow(simulatedDataset$exampleData)
        rownames(simulatedDataset$DevDataTmp) <- 1:nrow(simulatedDataset$DevDataTmp)
        rownames(simulatedDataset$exampleRates) <- 1:nrow(simulatedDataset$exampleRates)
        for(i in 1:rep){rownames(simulatedDataset$singleReplicateExpressionData[[i]]) <- 1:nrow(simulatedDataset$singleReplicateExpressionData[[i]])}

        if(nrow(simulatedDataset$exampleData)>=nGenes){break}
        }
    }

    simulatedDataset$exampleData <- simulatedDataset$exampleData[1:nGenes,]
    simulatedDataset$DevDataTmp <- simulatedDataset$DevDataTmp[1:nGenes,]
    simulatedDataset$exampleRates <- simulatedDataset$exampleRates[1:nGenes,]
    simulatedDataset$singleReplicateExpressionData <- lapply(1:rep,function(i){simulatedDataset$singleReplicateExpressionData[[i]][1:nGenes,]})

    saveRDS(simulatedDataset,file=paste0("simulatedDataset_",simulatedDataset$name,".rds"))
    
    for(flag in Flags){
        inferedNumericalModel <- inferRates(expressionData=simulatedDataset$singleReplicateExpressionData # Expression data of the genes to be modeled.
                                           ,expressionDataDev=NULL # Standard deviations of the genes to be modeled.
                                           ,simulatedDataset=simulatedDataset # Simulated dataset if this is the case (just to produce real rates correlations).
                                           ,initialRates=initialRates # List of initial rates for optimization.
                                           ,TauFractions=TauFractions # Time points with cellular fractionation.
                                           ,TauPoly=TauPoly # Time points with polysomal profiling.
                                           ,TauTotal=NULL # Time points with total RNA profiling.
                                           ,cpus=cpus # Number of cpus.
                                           ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
                                           ,lowB=1e-6 # Lower boundary for the rates.
                                           ,upB=1e4 # Upper boundary for the rates.
                                           ,FlagDev=flag # Cost function.
                                           ,lambda=0.05 # Regularization strength.
                                           ,excludeSpecies=NULL # List of species to be excluded from the cost function.
                                           ,parFixed=NULL) # List of parameters to be excluded from the optimization.
        saveRDS(inferedNumericalModel,file=paste0("inferedNumericalModel_",simulatedDataset$name,".rds"))
    }
}
