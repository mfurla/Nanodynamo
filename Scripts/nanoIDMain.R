#### Comparison against merged
### Expression data
expressionLevels <- readRDS("/path/to/Results/Untreated/fullModel/expressionLevels.rds")

### Data aggregation
expressionData <- expressionLevels$normalizedCountsStatistics$Mean$mean
expressionData <- expressionData[,!grepl("^py",colnames(expressionData))]
expressionData <- cbind("Nascent"=apply(expressionData[,grepl("n0.33$",colnames(expressionData))],1,sum)
					   ,"PreExisting"=apply(expressionData[,grepl("p0.33$",colnames(expressionData))],1,sum))

### nano-ID
alpha = log(2)/as.numeric(60*12)
labelingTime = 20/60

decay.rate.single.moleculeTmp = - alpha - (1/labelingTime)*log(1 - expressionData[,"Nascent"]/apply(expressionData,1,sum))
decay.rate.single.moleculeTmp[decay.rate.single.moleculeTmp <= 0] = NA
			
## Synthesis rate
synthesis.rate.single.moleculeTmp = apply(expressionData,1,sum)*(alpha + decay.rate.single.moleculeTmp)

## Half life
half.lives.single.moleculeTmp = log(2)/decay.rate.single.moleculeTmp

## Results
expressionData <- cbind(expressionData
					  ,"Synthesis"=synthesis.rate.single.moleculeTmp
					  ,"Degradation"=decay.rate.single.moleculeTmp
					  ,"HalfLife"=half.lives.single.moleculeTmp)
saveRDS(expressionData,"nanoIDResults.rds")