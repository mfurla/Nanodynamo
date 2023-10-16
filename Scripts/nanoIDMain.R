### Expression data
expressionLevelsUntreated <- readRDS("/path/to/Results/Untreated/expressionLevelsMerged.rds")

### Data aggregation
mergedExpressionLevelsUntreated <- sapply(c("chromatin"
										   ,"nucleoplasmic"
										   ,"cytoplasmic")
,function(i){
	print(i)

	j <- expressionLevelsUntreated[[i]][[1]]
	j <- j[,!grepl("0$",colnames(j))]
	j <- j[,!grepl("t$",colnames(j))]

	pTmp <- j[,grep("p0.33",colnames(j))]
	nTmp <- j[,grep("n0.33",colnames(j))]

	if(is.data.frame(pTmp))
	{
		return(cbind("PreExisting"=apply(pTmp,1,sum),"Nascent"=apply(nTmp,1,sum)))
	}else{
		names(pTmp) <- names(nTmp) <- rownames(j)
		return(cbind("PreExisting"=pTmp,"Nascent"=nTmp))
	}
})

commonGenes <- names(which(table(unlist(sapply(mergedExpressionLevelsUntreated,rownames)))==3))
mergedExpressionLevelsUntreated <- mergedExpressionLevelsUntreated[[1]][commonGenes,]+mergedExpressionLevelsUntreated[[2]][commonGenes,]+mergedExpressionLevelsUntreated[[3]][commonGenes,]

### nano-ID modeling
alpha = log(2)/as.numeric(60*12)
labelingTime = 20/60

decay.rate.single.moleculeTmp = - alpha - (1/labelingTime)*log(1 - mergedExpressionLevelsUntreated[,"Nascent"]/apply(mergedExpressionLevelsUntreated,1,sum))
decay.rate.single.moleculeTmp[decay.rate.single.moleculeTmp <= 0] = NA
			
## Synthesis rate
synthesis.rate.single.moleculeTmp = apply(mergedExpressionLevelsUntreated,1,sum)*(alpha + decay.rate.single.moleculeTmp)

## Half life
half.lives.single.moleculeTmp = log(2)/decay.rate.single.moleculeTmp

## Results
nanoIDResults <- cbind(mergedExpressionLevelsUntreated
								 ,"Synthesis"=synthesis.rate.single.moleculeTmp
								 ,"Degradation"=decay.rate.single.moleculeTmp
								 ,"HalfLife"=half.lives.single.moleculeTmp)
saveRDS(nanoIDResults,"nanoIDResults.rds")