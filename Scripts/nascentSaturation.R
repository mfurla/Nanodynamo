### Libraries
library("deSolve")
source("/path/to/allInternalFunctions.R")

### Custom functions
averagePlotFunction <- function(simulatedDataset,LOG=FALSE,abline=NULL)
{
	xTmp <- simulatedDataset$exampleData
	t_profiled <- simulatedDataset$TauGlobals

	yTmp <- cbind(apply(xTmp[,grep("^chpn",colnames(xTmp))]/xTmp[,"chpt"],2,mean)
				 ,apply(xTmp[,grep("^chmn",colnames(xTmp))]/xTmp[,"chmt"],2,mean)
				 ,apply(xTmp[,grep("^npn",colnames(xTmp))]/xTmp[,"npt"],2,mean)
				 ,apply(xTmp[,grep("^nmn",colnames(xTmp))]/xTmp[,"nmt"],2,mean)
				 ,apply(xTmp[,grep("^cyn",colnames(xTmp))]/xTmp[,"cyt"],2,mean)
				 ,apply(xTmp[,grep("^pyn",colnames(xTmp))]/xTmp[,"pyt"],2,mean))

	xTmp <- cbind(apply(xTmp[,grep("^chpn",colnames(xTmp))],2,mean)
				 ,apply(xTmp[,grep("^chmn",colnames(xTmp))],2,mean)
				 ,apply(xTmp[,grep("^npn",colnames(xTmp))],2,mean)
				 ,apply(xTmp[,grep("^nmn",colnames(xTmp))],2,mean)
				 ,apply(xTmp[,grep("^cyn",colnames(xTmp))],2,mean)
				 ,apply(xTmp[,grep("^pyn",colnames(xTmp))],2,mean))

	print(yTmp)

	if(LOG)
	{
		matplot(format(signif(log10(t_profiled)[-1],3),nsmall=2),yTmp[-1,],xlab="Log10(Time [h])",ylab="Nascent RNA [%]",ylim=c(0,1.05),las=2,col=c(1:6),type="l",lwd=2,lty=1)			
		abline(v=log10(abline),col="grey",lty=2,lwd=2)
	}else{
		matplot(format(signif(t_profiled,3),nsmall=2),yTmp,xlab="Time [h]",ylab="Nascent RNA [%]",ylim=c(0,1.05),las=2,col=c(1:6),type="l",lwd=2,lty=1)
		abline(v=abline,col="grey",lty=2,lwd=2)
	}
	legend("bottomright",col=c(1:6,"grey"),lwd=2,lty=c(rep(1,6),2),legend=c("Chp","Chm","Np","Nm","C","P","20 min."),bty="n")
}

### Time points
timePoints <- round(c(0,1,2,3,4,5,10,15,20,25,30,45,60,90,120)/60,2)

### CVs (set to 0)
CVspecies = readRDS("/path/to/Files/CVs.rds")

## Add CVs for additional time points (1h and 2h)
CVspecies = lapply(CVspecies,function(i)
{
	addTmp <- lapply(timePoints,function(j)
	{
		i[grep("0.33",names(i))]	
	})

	addTmp <- lapply(seq_along(addTmp),function(j)
	{
  		outTmp <- addTmp[[j]]
  		names(outTmp) <- gsub("0.33",timePoints[[j]],names(outTmp))
  		outTmp
	})

  outTmp <- c(i,unlist(addTmp))
  outTmp[outTmp!=0] <- 0
  outTmp[unique(names(outTmp))]
})

### Average gene
## Data from the supplemental material of de Pretis et al. Genome Research 2017
INSPEcT_nascent <- read.table("/path/to/Files/regulated_genes_features.xls",sep="\t",header=TRUE)

## Steady state rates
synthesis = signif(median(INSPEcT_nascent[,"synthesis_0"]),2)
processing = signif(median(INSPEcT_nascent[,"processing_0"]),2)
degradation = signif(median(INSPEcT_nascent[,"degradation_0"]),2)

## Rates distributions means
exampleRates <- c(k1 = synthesis # Synthesis
                 ,k2 = processing # Co-transcriptional processing
                 ,k3 = processing # Detachment of mature chromatin
                 ,k4 = 0.1*processing # Detachment of premature chromatin
                 ,k5 = 0.25*processing # Post-transcriptional processing 
                 ,k6 = 1.5*processing # Export
                 ,k7 = 0.5*degradation # Association with polysomes
                 ,k8 = degradation # Cytoplasmic degradation
                 ,k10 = degradation) # Polysomal degradation

### Data simulation
simulatedDataset <- simulateData(exampleRates = exampleRates # Mean rates to sample.
								,TauFractions = timePoints # Time points with cellular fractionation.
								,TauPoly = timePoints # Time points with polysomal profiling.
								,TauTotal = NULL # Time points with total RNA profiling (FIXED to null for our purposes).
								,noise = TRUE # TRUE to add noise to data.
								,CVs = CVspecies # Variation Coefficient for noise.
								,Reps = 1 # Number of replicates.
								,nGenes = 1000 # Number of genes to be simulated.
								,seed = 1 # Seed for reproducibility.
								,ZeroThresh = 1e-10 # Minimum expression value.
								,MultFact = 1e5) # Number of nGenes to be modeled.

### Plot - Figure S3
averagePlotFunction(simulatedDataset,LOG=FALSE,abline=c(20/60))

