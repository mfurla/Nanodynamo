### Libraries
library("INSPEcT")
library("pheatmap")

#### Comparison against merged
### Expression data
expressionLevels <- readRDS("/path/to/Results/Untreated/expressionLevels.rds")

### NanoID results
nanoID <- readRDS("/path/to/Files/nanoIDResults.rds")

### INSPEcT
## Data aggregation
expressionData <- expressionLevels$normalizedCountsStatistics$Mean$mean
expressionData <- expressionData[,!grepl("^py",colnames(expressionData))]
expressionData <- expressionData[,!grepl("0$",colnames(expressionData))]
expressionData <- expressionData[,!grepl("t$",colnames(expressionData))]

pTmp <- expressionData[,grep("p0.33",colnames(expressionData))]
nTmp <- expressionData[,grep("n0.33",colnames(expressionData))]

ppTmp <- apply(as.matrix(expressionData[,grep("pp0.33",colnames(expressionData))],nrow=nrow(expressionData)),1,sum)
pnTmp <- apply(as.matrix(expressionData[,grep("pn0.33",colnames(expressionData))],nrow=nrow(expressionData)),1,sum)

expressionData <- cbind("pt"=ppTmp+pnTmp,"pn"=pnTmp,"tt"=apply(pTmp,1,sum)+apply(nTmp,1,sum),"tn"=apply(nTmp,1,sum))

## Modeling
nascent <- quantifyExpressionsFromTrAbundance(trAbundaces = list(exonsAbundances=cbind("Rep1"=expressionData[,"tn"]
																					  ,"Rep2"=1e-10*expressionData[,"tn"] # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																					  ,"Rep3"=expressionData[,"tn"]
																					  ,"Rep4"=1e-10*expressionData[,"tn"]) # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																,intronsAbundances=cbind("Rep1"=expressionData[,"pn"]
																						,"Rep2"=1e-10*expressionData[,"pn"] # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																						,"Rep3"=expressionData[,"pn"]
																						,"Rep4"=1e-10*expressionData[,"pn"])) # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
											, experimentalDesign = c(Rep1="WT",Rep2="WT",Rep3="WTBis",Rep4="WTBis"))

total <- quantifyExpressionsFromTrAbundance(trAbundaces = list(exonsAbundances=cbind("Rep1"=expressionData[,"tt"]
																					,"Rep2"=1e-10*expressionData[,"tt"] # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																					,"Rep3"=expressionData[,"tt"]
																					,"Rep4"=1e-10*expressionData[,"tt"]) # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																,intronsAbundances=cbind("Rep1"=expressionData[,"pt"]
																						,"Rep2"=1e-10*expressionData[,"pt"] # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
																						,"Rep3"=expressionData[,"pt"]
																						,"Rep4"=1e-10*expressionData[,"pt"])) # Simulated replicate - Required to estimate the variance which, however, is not involved in rates first guess estimation
											, experimentalDesign = c(Rep1="WT",Rep2="WT",Rep3="WTBis",Rep4="WTBis"))

nascentInspObj <- newINSPEcT(tpts=c("WT","WTBis")
							,labeling_time=20/60
							,nascentExpressions=nascent
							,matureExpressions=total)

INSPEcT_synthesis <- log10(ratesFirstGuess(nascentInspObj,"synthesis")[,"synthesis_WT"])
INSPEcT_processing <- log10(ratesFirstGuess(nascentInspObj,"processing")[,"processing_WT"])
INSPEcT_degradation <- log10(ratesFirstGuess(nascentInspObj,"degradation")[,"degradation_WT"])

### Comparisons with Nanodynamo
## Same genes different models
inferedRatesUntreatedMerged_yesChpNpP<-readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreatedMerged_YesChpNpP_multi.rds")
inferedRatesUntreatedMerged_yesChpNpNoP<-readRDS("/path/to/Results/Untreated/thirdRun_sameGenesAllModels/inferedRatesUntreatedMerged_yesChpNpNoP_multi.rds")
inferedRatesUntreatedMerged_yesChpPNoNp<-readRDS("/path/to/Results/Untreated/thirdRun_sameGenesAllModels/inferedRatesUntreatedMerged_yesChpPNoNp_multi.rds")
inferedRatesUntreatedMerged_yesChpNoNpP<-readRDS("/path/to/Results/Untreated/thirdRun_sameGenesAllModels/inferedRatesUntreatedMerged_yesChpNoNpP_multi.rds")
inferedRatesUntreatedMerged_yesPNoChp<-readRDS("/path/to/Results/Untreated/thirdRun_sameGenesAllModels/inferedRatesUntreatedMerged_yesPNoChp_multi.rds")
inferedRatesUntreatedMerged_noChpP<-readRDS("/path/to/Results/Untreated/thirdRun_sameGenesAllModels/inferedRatesUntreatedMerged_noChpP_multi.rds")

sameGenesSimplerModelsResults <- lapply(paste0("k",c(1,3,6,8)),function(i)
{
	cbind("Chp_Np_P"=inferedRatesUntreatedMerged_yesChpNpP$inferedRates[,i]
		 ,"Chp_Np"=inferedRatesUntreatedMerged_yesChpNpNoP$inferedRates[,i]
		 ,"Chp_P"=inferedRatesUntreatedMerged_yesChpPNoNp$inferedRates[,i]
		 ,"Chp"=inferedRatesUntreatedMerged_yesChpNoNpP$inferedRates[,i]
		 ,"P"=inferedRatesUntreatedMerged_yesPNoChp$inferedRates[,i]
		 ,"-"=inferedRatesUntreatedMerged_noChpP$inferedRates[,i])
})

names(sameGenesSimplerModelsResults) <- c(paste0("k",c(1,3,6,8)))

sameGenesSimplerModelsResults[[1]] <- cbind(sameGenesSimplerModelsResults[[1]]
										   ,"INSPEcT"=INSPEcT_synthesis[rownames(sameGenesSimplerModelsResults[[1]])]
										   ,"nano-ID"=nanoID[rownames(sameGenesSimplerModelsResults[[1]]),"Synthesis"])
sameGenesSimplerModelsResults[[4]] <- cbind(sameGenesSimplerModelsResults[[4]]
										   ,"INSPEcT"=INSPEcT_degradation[rownames(sameGenesSimplerModelsResults[[4]])]
										   ,"nano-ID"=nanoID[rownames(sameGenesSimplerModelsResults[[4]]),"Degradation"])

sameGenesSimplerModelsResultsCor <- lapply(sameGenesSimplerModelsResults,function(i)round(cor(i,method="s"),2))

### Figures S53, S54, S55, S56
for(i in names(sameGenesSimplerModelsResultsCor))
{
	pheatmap(sameGenesSimplerModelsResultsCor[[i]]
		    ,breaks=seq(0,1,length.out=20)
		    ,color=colorRampPalette(c("white","red"))(20)
		    ,width=3.5
		    ,height=3.5
		    ,filename=paste0(i,"Heatmap.pdf")
		    ,display_numbers=TRUE
		    ,cluster_rows=FALSE
		    ,cluster_cols=FALSE
		    ,main=i
		    ,legend_labels="Spearman correlation")
}

## Different genes
inferedRatesUntreatedMerged_yesChpNpP<-readRDS("/path/to/Results/Untreated/firstRun_FullModel/inferedRatesUntreatedMerged_YesChpNpP_multi.rds")
inferedRatesUntreatedMerged_yesChpNpNoP<-readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_YesChpNpNoP_multi.rds")
inferedRatesUntreatedMerged_yesChpPNoNp<-readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_YesChpPNoNp_multi.rds")
inferedRatesUntreatedMerged_yesChpNoNpP<-readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_YesChpNoNpP_multi.rds")
inferedRatesUntreatedMerged_yesPNoChp<-readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_YesPNoChp_multi.rds")
inferedRatesUntreatedMerged_noChpP<-readRDS("/path/to/Results/Untreated/secondRun_simplerModels/inferedRatesUntreatedMerged_NoChpP_multi.rds")

INSPEcT2_synthesis <- list("Nanodynamo"=log10(c(inferedRatesUntreatedMerged_yesChpNpP$inferedRates[,"k1"]
											   ,inferedRatesUntreatedMerged_yesChpPNoNp$inferedRates[,"k1"]
											   ,inferedRatesUntreatedMerged_yesPNoChp$inferedRates[,"k1"]
											   ,inferedRatesUntreatedMerged_yesChpNpNoP$inferedRates[,"k1"]
											   ,inferedRatesUntreatedMerged_yesChpNoNpP$inferedRates[,"k1"]
											   ,inferedRatesUntreatedMerged_noChpP$inferedRates[,"k1"]))
						  ,"NanoID"=log10(nanoID[,"Synthesis"]))

INSPEcT2_degradation <- list("Nanodynamo"=log10(c(inferedRatesUntreatedMerged_yesChpNpP$inferedRates[,"k8"]
											   ,inferedRatesUntreatedMerged_yesChpPNoNp$inferedRates[,"k8"]
											   ,inferedRatesUntreatedMerged_yesPNoChp$inferedRates[,"k8"]
											   ,inferedRatesUntreatedMerged_yesChpNpNoP$inferedRates[,"k8"]
											   ,inferedRatesUntreatedMerged_yesChpNoNpP$inferedRates[,"k8"]
											   ,inferedRatesUntreatedMerged_noChpP$inferedRates[,"k8"]))
							,"NanoID"=log10(nanoID[,"Degradation"]))

commonGenes <- data.frame(sapply(names(INSPEcT2_synthesis),function(j)
{
	i <- INSPEcT2_synthesis[[j]]
	commonGenesTmp <- intersect(names(i),names(INSPEcT_synthesis))
	log10(length(commonGenesTmp))
}))
colnames(commonGenes) <- "N. Genes"

synthesisCor <- sapply(names(INSPEcT2_synthesis),function(j)
{
	i <- INSPEcT2_synthesis[[j]]
	commonGenesTmp <- intersect(names(i),names(INSPEcT_synthesis))
	print(length(commonGenesTmp))
	round(cor(INSPEcT_synthesis[commonGenesTmp],i[commonGenesTmp],method="s",use="c"),2)
})

degradationCor <- sapply(names(INSPEcT2_degradation),function(j)
{
	i <- INSPEcT2_degradation[[j]]
	commonGenesTmp <- intersect(names(i),names(INSPEcT_degradation))
	print(length(commonGenesTmp))
	round(cor(INSPEcT_degradation[commonGenesTmp],i[commonGenesTmp],method="s",use="c"),2)
})

