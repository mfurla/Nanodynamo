### Libraries
library("INSPEcT")
library("pheatmap")

### Expression data
expressionLevelsUntreated <- readRDS("/path/to/Results/Untreated/expressionLevelsMerged.rds")

### NanoID results
nanoID <- readRDS("/path/to/Results/nanoIDRates.rds")

### Data aggregation
mergedExpressionLevelsUntreated <- sapply(c("normalizedCountsYesChpNpP"
										   ,"normalizedCountsYesChpPNoNp"
										   ,"normalizedCountsYesPNoChp"
										   ,"normalizedCountsYesChpNpNoP"
										   ,"normalizedCountsYesChpNoNpP"
										   ,"normalizedCountsNoChpP")
,function(i){
	print(i)

	j <- expressionLevelsUntreated[[i]][["Sample_Merged"]]
	j <- j[,!grepl("^py",colnames(j))]
	j <- j[,!grepl("0$",colnames(j))]
	j <- j[,!grepl("t$",colnames(j))]

	pTmp <- j[,grep("p0.33",colnames(j))]
	nTmp <- j[,grep("n0.33",colnames(j))]

	ppTmp <- tryCatch(apply(as.matrix(j[,grep("pp0.33",colnames(j))],nrow=nrow(j)),1,sum),error=function(e){rep(0,nrow(j))})
	pnTmp <- tryCatch(apply(as.matrix(j[,grep("pn0.33",colnames(j))],nrow=nrow(j)),1,sum),error=function(e){rep(0,nrow(j))})
	
	cbind("pt"=ppTmp+pnTmp,"pn"=pnTmp,"tt"=apply(pTmp,1,sum)+apply(nTmp,1,sum),"tn"=apply(nTmp,1,sum))
})
mergedExpressionLevelsUntreated <- do.call("rbind",mergedExpressionLevelsUntreated)

### INSPEcT
## Expression levels quantification
nascent <- quantifyExpressionsFromTrAbundance(trAbundaces = list(exonsAbundances=cbind("Rep1"=mergedExpressionLevelsUntreated[,"tn"]
																					  ,"Rep2"=1e-10*mergedExpressionLevelsUntreated[,"tn"] # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																					  ,"Rep3"=mergedExpressionLevelsUntreated[,"tn"]
																					  ,"Rep4"=1e-10*mergedExpressionLevelsUntreated[,"tn"]) # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																,intronsAbundances=cbind("Rep1"=mergedExpressionLevelsUntreated[,"pn"]
																						,"Rep2"=1e-10*mergedExpressionLevelsUntreated[,"pn"] # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																						,"Rep3"=mergedExpressionLevelsUntreated[,"pn"]
																						,"Rep4"=1e-10*mergedExpressionLevelsUntreated[,"pn"])) # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
											, experimentalDesign = c(Rep1="Untreated",Rep2="Untreated",Rep3="UntreatedBis",Rep4="UntreatedBis"))

total <- quantifyExpressionsFromTrAbundance(trAbundaces = list(exonsAbundances=cbind("Rep1"=mergedExpressionLevelsUntreated[,"tt"]
																					,"Rep2"=1e-10*mergedExpressionLevelsUntreated[,"tt"] # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																					,"Rep3"=mergedExpressionLevelsUntreated[,"tt"]
																					,"Rep4"=1e-10*mergedExpressionLevelsUntreated[,"tt"]) # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																,intronsAbundances=cbind("Rep1"=mergedExpressionLevelsUntreated[,"pt"]
																						,"Rep2"=1e-10*mergedExpressionLevelsUntreated[,"pt"] # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																						,"Rep3"=mergedExpressionLevelsUntreated[,"pt"]
																						,"Rep4"=1e-10*mergedExpressionLevelsUntreated[,"pt"])) # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
											, experimentalDesign = c(Rep1="Untreated",Rep2="Untreated",Rep3="UntreatedBis",Rep4="UntreatedBis"))

## Main INSPEcT function
nascentInspObj <- newINSPEcT(tpts=c("Untreated","UntreatedBis")
							,labeling_time=20/60
							,nascentExpressions=nascent
							,matureExpressions=total)

### Comparative analysis
## INSPEcT rates
INSPEcT_synthesis <- log10(ratesFirstGuess(nascentInspObj,"synthesis")[,"synthesis_Untreated"])
INSPEcT_processing <- log10(ratesFirstGuess(nascentInspObj,"processing")[,"processing_Untreated"])
INSPEcT_degradation <- log10(ratesFirstGuess(nascentInspObj,"degradation")[,"degradation_Untreated"])

## Nanodynamo rates
inferedRatesUntreated_NoChpP <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_noChpP.rds")
inferedRatesUntreated_YesChpNoNpP <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_yesChpNoNpP.rds")
inferedRatesUntreated_YesChpNpNoP <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_yesChpNpNoP.rds")
inferedRatesUntreated_YesChpNpP <- readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreatedMerged_yesChpNpP_multi.rds")
inferedRatesUntreated_YesChpPNoNp <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_yesChpPNoNp.rds")
inferedRatesUntreated_YesPNoChp <- readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_yesPNoChp.rds")

Nanodynamo_synthesis <- list("WithP"=log10(c(inferedRatesWT_YesChpNpP$inferedRates[,"k1"]
											,inferedRatesWT_YesChpPNoNp$inferedRates[,"k1"]
											,inferedRatesWT_YesPNoChp$inferedRates[,"k1"]))
							,"WithOutP"=log10(c(inferedRatesWT_YesChpNpNoP$inferedRates[,"k1"]
											,inferedRatesWT_YesChpNoNpP$inferedRates[,"k1"]
											,inferedRatesWT_NoChpP$inferedRates[,"k1"]))
							,"NanoID"=log10(nanoID[,"Synthesis"]))

Nanodynamo_degradation <- list("WithP"=log10(c(inferedRatesUntreated_YesChpNpP$inferedRates[,"k8"]
						  					  ,inferedRatesUntreated_YesChpPNoNp$inferedRates[,"k8"]
						  					  ,inferedRatesUntreated_YesPNoChp$inferedRates[,"k8"]))
							  ,"WithOutP"=log10(c(inferedRatesUntreated_YesChpNpNoP$inferedRates[,"k8"]
											     ,inferedRatesUntreated_YesChpNoNpP$inferedRates[,"k8"]
											     ,inferedRatesUntreated_NoChpP$inferedRates[,"k8"]))
							  ,"NanoID"=log10(nanoID[,"Degradation"]))

sapply(names(Nanodynamo_synthesis),function(j)
{
	i <- Nanodynamo_synthesis[[j]]
	commonGenesTmp <- intersect(names(i),names(INSPEcT_synthesis))
	print(length(commonGenesTmp))
	round(cor(INSPEcT_synthesis[commonGenesTmp],i[commonGenesTmp],method="s",use="c"),2)
})

# WithP WithOutP NanoID
#  0.95     0.81   0.99

sapply(names(Nanodynamo_degradation),function(j)
{
	i <- Nanodynamo_degradation[[j]]
	commonGenesTmp <- intersect(names(i),names(INSPEcT_degradation))
	print(length(commonGenesTmp))
	round(cor(INSPEcT_degradation[commonGenesTmp],i[commonGenesTmp],method="s",use="c"),2)
})

# WithP WithOutP NanoID 
#  0.11     0.66   0.97