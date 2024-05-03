###########################################################################################
### Code for the comparison of nascent and premature RNA among fractions and conditions ###
###########################################################################################

#### Untreated condition - nascent and premature across conditions
### Expression levels
expressionLevelsUntreated <- readRDS("/path/to/Results/Untreated/expressionLevels.rds")
fullModelGenes <- rownames(expressionLevelsUntreated$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$WT1$mean)

### Nascent RNA
## First replicate
chpn1 <- expressionLevelsUntreated$rawCounts$Chr_WT1[,"chpn0.33"]/expressionLevelsUntreated$rawCounts$Chr_WT1[,"chpt"]
names(chpn1) <- rownames(expressionLevelsUntreated$rawCounts$Chr_WT1)
chmn1 <- expressionLevelsUntreated$rawCounts$Chr_WT1[,"chmn0.33"]/expressionLevelsUntreated$rawCounts$Chr_WT1[,"chmt"]
names(chmn1) <- rownames(expressionLevelsUntreated$rawCounts$Chr_WT1)

npn1 <- expressionLevelsUntreated$rawCounts$Nuc_WT1[,"npn0.33"]/expressionLevelsUntreated$rawCounts$Nuc_WT1[,"npt"]
names(npn1) <- rownames(expressionLevelsUntreated$rawCounts$Nuc_WT1)
nmn1 <- expressionLevelsUntreated$rawCounts$Nuc_WT1[,"nmn0.33"]/expressionLevelsUntreated$rawCounts$Nuc_WT1[,"nmt"]
names(nmn1) <- rownames(expressionLevelsUntreated$rawCounts$Nuc_WT1)

cyn1 <- (expressionLevelsUntreated$rawCounts$Cyt_WT1[,"cypn0.33"]+expressionLevelsUntreated$rawCounts$Cyt_WT1[,"cymn0.33"])/(expressionLevelsUntreated$rawCounts$Cyt_WT1[,"cypt"]+expressionLevelsUntreated$rawCounts$Cyt_WT1[,"cymt"])
names(cyn1) <- rownames(expressionLevelsUntreated$rawCounts$Cyt_WT1)
pyn1 <- (expressionLevelsUntreated$rawCounts$Poly_WT1[,"pypn0.33"]+expressionLevelsUntreated$rawCounts$Poly_WT1[,"pymn0.33"])/(expressionLevelsUntreated$rawCounts$Poly_WT1[,"pypt"]+expressionLevelsUntreated$rawCounts$Poly_WT1[,"pymt"])
names(pyn1) <- rownames(expressionLevelsUntreated$rawCounts$Poly_WT1)

## Second replicate
chpn2 <- expressionLevelsUntreated$rawCounts$Chr_WT2[,"chpn0.33"]/expressionLevelsUntreated$rawCounts$Chr_WT2[,"chpt"]
names(chpn2) <- rownames(expressionLevelsUntreated$rawCounts$Chr_WT2)
chmn2 <- expressionLevelsUntreated$rawCounts$Chr_WT2[,"chmn0.33"]/expressionLevelsUntreated$rawCounts$Chr_WT2[,"chmt"]
names(chmn2) <- rownames(expressionLevelsUntreated$rawCounts$Chr_WT2)

npn2 <- expressionLevelsUntreated$rawCounts$Nuc_WT2[,"npn0.33"]/expressionLevelsUntreated$rawCounts$Nuc_WT2[,"npt"]
names(npn2) <- rownames(expressionLevelsUntreated$rawCounts$Nuc_WT2)
nmn2 <- expressionLevelsUntreated$rawCounts$Nuc_WT2[,"nmn0.33"]/expressionLevelsUntreated$rawCounts$Nuc_WT2[,"nmt"]
names(nmn2) <- rownames(expressionLevelsUntreated$rawCounts$Nuc_WT2)

cyn2 <- (expressionLevelsUntreated$rawCounts$Cyt_WT2[,"cypn0.33"]+expressionLevelsUntreated$rawCounts$Cyt_WT2[,"cymn0.33"])/(expressionLevelsUntreated$rawCounts$Cyt_WT2[,"cypt"]+expressionLevelsUntreated$rawCounts$Cyt_WT2[,"cymt"])
names(cyn2) <- rownames(expressionLevelsUntreated$rawCounts$Cyt_WT2)

pyn2 <- (expressionLevelsUntreated$rawCounts$Poly_WT2[,"pypn0.33"]+expressionLevelsUntreated$rawCounts$Poly_WT2[,"pymn0.33"])/(expressionLevelsUntreated$rawCounts$Poly_WT2[,"pypt"]+expressionLevelsUntreated$rawCounts$Poly_WT2[,"pymt"])
names(pyn2) <- rownames(expressionLevelsUntreated$rawCounts$Poly_WT2)

### Premature RNA
## First replicate
ch1 <- expressionLevelsUntreated$rawCounts$Chr_WT1[,"chpt"]/(expressionLevelsUntreated$rawCounts$Chr_WT1[,"chpt"]+expressionLevelsUntreated$rawCounts$Chr_WT1[,"chmt"])
names(ch1) <- rownames(expressionLevelsUntreated$rawCounts$Chr_WT1)

n1 <- expressionLevelsUntreated$rawCounts$Nuc_WT1[,"npt"]/(expressionLevelsUntreated$rawCounts$Nuc_WT1[,"npt"]+expressionLevelsUntreated$rawCounts$Nuc_WT1[,"nmt"])
names(n1) <- rownames(expressionLevelsUntreated$rawCounts$Nuc_WT1)

cy1 <- expressionLevelsUntreated$rawCounts$Cyt_WT1[,"cypt"]/(expressionLevelsUntreated$rawCounts$Cyt_WT1[,"cypt"]+expressionLevelsUntreated$rawCounts$Cyt_WT1[,"cymt"])
names(cy1) <- rownames(expressionLevelsUntreated$rawCounts$Cyt_WT1)

py1 <- expressionLevelsUntreated$rawCounts$Poly_WT1[,"pypt"]/(expressionLevelsUntreated$rawCounts$Poly_WT1[,"pypt"]+expressionLevelsUntreated$rawCounts$Poly_WT1[,"pymt"])
names(py1) <- rownames(expressionLevelsUntreated$rawCounts$Poly_WT1)

## Second replicate
ch2 <- expressionLevelsUntreated$rawCounts$Chr_WT2[,"chpt"]/(expressionLevelsUntreated$rawCounts$Chr_WT2[,"chpt"]+expressionLevelsUntreated$rawCounts$Chr_WT2[,"chmt"])
names(ch2) <- rownames(expressionLevelsUntreated$rawCounts$Chr_WT2)

n2 <- expressionLevelsUntreated$rawCounts$Nuc_WT2[,"npt"]/(expressionLevelsUntreated$rawCounts$Nuc_WT2[,"npt"]+expressionLevelsUntreated$rawCounts$Nuc_WT2[,"nmt"])
names(n2) <- rownames(expressionLevelsUntreated$rawCounts$Nuc_WT2)

cy2 <- expressionLevelsUntreated$rawCounts$Cyt_WT2[,"cypt"]/(expressionLevelsUntreated$rawCounts$Cyt_WT2[,"cypt"]+expressionLevelsUntreated$rawCounts$Cyt_WT2[,"cymt"])
names(cy2) <- rownames(expressionLevelsUntreated$rawCounts$Cyt_WT2)

py2 <- expressionLevelsUntreated$rawCounts$Poly_WT2[,"pypt"]/(expressionLevelsUntreated$rawCounts$Poly_WT2[,"pypt"]+expressionLevelsUntreated$rawCounts$Poly_WT2[,"pymt"])
names(py2) <- rownames(expressionLevelsUntreated$rawCounts$Poly_WT2)

### Removal of species under the minimim expression threshold
chpn1[chpn1<=1e-10] <- NaN;nmn1[nmn1<=1e-10] <- NaN;npn2[npn2<=1e-10] <- NaN;pyn1[pyn1<=1e-10] <- NaN;ch1[ch1<=1e-10] <- NaN;py1[py1<=1e-10] <- NaN;cy2[cy2<=1e-10] <- NaN
chmn1[chmn1<=1e-10] <- NaN;chpn2[chpn2<=1e-10] <- NaN;nmn2[nmn2<=1e-10] <- NaN;cyn2[cyn2<=1e-10] <- NaN;n1[n1<=1e-10] <- NaN;ch2[ch2<=1e-10] <- NaN;py2[py2<=1e-10] <- NaN
npn1[npn1<=1e-10] <- NaN;chmn2[chmn2<=1e-10] <- NaN;cyn1[cyn1<=1e-10] <- NaN;pyn2[pyn2<=1e-10] <- NaN;cy1[cy1<=1e-10] <- NaN;n2[n2<=1e-10] <- NaN

### Figure S10
par(mfrow=c(1,2))
boxplot(list(Chp=apply(cbind(chpn1[fullModelGenes],chpn2[fullModelGenes]),1,mean,na.rm=TRUE)
			,Np=apply(cbind(npn1[fullModelGenes],npn2[fullModelGenes]),1,mean,na.rm=TRUE)
			,Chm=apply(cbind(chmn1[fullModelGenes],chmn2[fullModelGenes]),1,mean,na.rm=TRUE)
			,Nm=apply(cbind(nmn1[fullModelGenes],nmn2[fullModelGenes]),1,mean,na.rm=TRUE)
			,C=apply(cbind(cyn1[fullModelGenes],cyn2[fullModelGenes]),1,mean,na.rm=TRUE)
			,P=apply(cbind(pyn1[fullModelGenes],pyn2[fullModelGenes]),1,mean,na.rm=TRUE))
		,outline=FALSE,ylab="Nascent RNA [%]",ylim=c(0,1),varwidth=TRUE,las=2)
boxplot(list(Ch=apply(cbind(ch1[fullModelGenes],ch2[fullModelGenes]),1,mean,na.rm=TRUE)
			,N=apply(cbind(n1[fullModelGenes],n2[fullModelGenes]),1,mean,na.rm=TRUE)
			,C=apply(cbind(cy1[fullModelGenes],cy2[fullModelGenes]),1,mean,na.rm=TRUE)
			,P=apply(cbind(py1[fullModelGenes],py2[fullModelGenes]),1,mean,na.rm=TRUE))
		,outline=FALSE,ylab="Premature RNA [%]",ylim=c(0,0.6),varwidth=TRUE,las=2)

### Numerical values
## Nascent
round(median(apply(cbind(chpn1[fullModelGenes],chpn2[fullModelGenes]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.45
round(median(apply(cbind(npn1[fullModelGenes],npn2[fullModelGenes]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.42
round(median(apply(cbind(chmn1[fullModelGenes],chmn2[fullModelGenes]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.33
round(median(apply(cbind(nmn1[fullModelGenes],nmn2[fullModelGenes]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.24
round(median(apply(cbind(cyn1[fullModelGenes],cyn2[fullModelGenes]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.27
round(median(apply(cbind(pyn1[fullModelGenes],pyn2[fullModelGenes]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.28

## Premature
round(median(apply(cbind(ch1[fullModelGenes],ch2[fullModelGenes]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.15
round(median(apply(cbind(n1[fullModelGenes],n2[fullModelGenes]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.07
round(median(apply(cbind(cy1[fullModelGenes],cy2[fullModelGenes]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.08
round(median(apply(cbind(py1[fullModelGenes],py2[fullModelGenes]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.07

#### PlaB vs Untreated - premature and nascent
## PlaB expression levels
expressionLevelsPlaB <- readRDS("/path/to/Results/PladienolideB/expressionLevels.rds")
fullModelGenes_PlaB <- rownames(expressionLevelsPlaB$normalizedCountsStatisticsSplitted$normalizedCountsTmpYesChpNpP$PlaB1$mean)

## Nascent RNA
# First replicate
chpn1_PlaB <- expressionLevelsPlaB$rawCounts$Chr_PlaB1[,"chpn0.33"]/expressionLevelsPlaB$rawCounts$Chr_PlaB1[,"chpt"]
names(chpn1_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Chr_PlaB1)
chmn1_PlaB <- expressionLevelsPlaB$rawCounts$Chr_PlaB1[,"chmn0.33"]/expressionLevelsPlaB$rawCounts$Chr_PlaB1[,"chmt"]
names(chmn1_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Chr_PlaB1)

npn1_PlaB <- expressionLevelsPlaB$rawCounts$Nuc_PlaB1[,"npn0.33"]/expressionLevelsPlaB$rawCounts$Nuc_PlaB1[,"npt"]
names(npn1_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Nuc_PlaB1)
nmn1_PlaB <- expressionLevelsPlaB$rawCounts$Nuc_PlaB1[,"nmn0.33"]/expressionLevelsPlaB$rawCounts$Nuc_PlaB1[,"nmt"]
names(nmn1_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Nuc_PlaB1)

cyn1_PlaB <- (expressionLevelsPlaB$rawCounts$Cyt_PlaB1[,"cypn0.33"]+expressionLevelsPlaB$rawCounts$Cyt_PlaB1[,"cymn0.33"])/(expressionLevelsPlaB$rawCounts$Cyt_PlaB1[,"cypt"]+expressionLevelsPlaB$rawCounts$Cyt_PlaB1[,"cymt"])
names(cyn1_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Cyt_PlaB1)
pyn1_PlaB <- (expressionLevelsPlaB$rawCounts$Poly_PlaB1[,"pypn0.33"]+expressionLevelsPlaB$rawCounts$Poly_PlaB1[,"pymn0.33"])/(expressionLevelsPlaB$rawCounts$Poly_PlaB1[,"pypt"]+expressionLevelsPlaB$rawCounts$Poly_PlaB1[,"pymt"])
names(pyn1_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Poly_PlaB1)

# Second replicate
chpn2_PlaB <- expressionLevelsPlaB$rawCounts$Chr_PlaB2[,"chpn0.33"]/expressionLevelsPlaB$rawCounts$Chr_PlaB2[,"chpt"]
names(chpn2_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Chr_PlaB2)
chmn2_PlaB <- expressionLevelsPlaB$rawCounts$Chr_PlaB2[,"chmn0.33"]/expressionLevelsPlaB$rawCounts$Chr_PlaB2[,"chmt"]
names(chmn2_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Chr_PlaB2)

npn2_PlaB <- expressionLevelsPlaB$rawCounts$Nuc_PlaB2[,"npn0.33"]/expressionLevelsPlaB$rawCounts$Nuc_PlaB2[,"npt"]
names(npn2_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Nuc_PlaB2)
nmn2_PlaB <- expressionLevelsPlaB$rawCounts$Nuc_PlaB2[,"nmn0.33"]/expressionLevelsPlaB$rawCounts$Nuc_PlaB2[,"nmt"]
names(nmn2_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Nuc_PlaB2)

cyn2_PlaB <- (expressionLevelsPlaB$rawCounts$Cyt_PlaB2[,"cypn0.33"]+expressionLevelsPlaB$rawCounts$Cyt_PlaB2[,"cymn0.33"])/(expressionLevelsPlaB$rawCounts$Cyt_PlaB2[,"cypt"]+expressionLevelsPlaB$rawCounts$Cyt_PlaB2[,"cymt"])
names(cyn2_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Cyt_PlaB2)

pyn2_PlaB <- (expressionLevelsPlaB$rawCounts$Poly_PlaB2[,"pypn0.33"]+expressionLevelsPlaB$rawCounts$Poly_PlaB2[,"pymn0.33"])/(expressionLevelsPlaB$rawCounts$Poly_PlaB2[,"pypt"]+expressionLevelsPlaB$rawCounts$Poly_PlaB2[,"pymt"])
names(pyn2_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Poly_PlaB2)

## Premature RNA
# First replicate
ch1_PlaB <- expressionLevelsPlaB$rawCounts$Chr_PlaB1[,"chpt"]/(expressionLevelsPlaB$rawCounts$Chr_PlaB1[,"chpt"]+expressionLevelsPlaB$rawCounts$Chr_PlaB1[,"chmt"])
names(ch1_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Chr_PlaB1)

n1_PlaB <- expressionLevelsPlaB$rawCounts$Nuc_PlaB1[,"npt"]/(expressionLevelsPlaB$rawCounts$Nuc_PlaB1[,"npt"]+expressionLevelsPlaB$rawCounts$Nuc_PlaB1[,"nmt"])
names(n1_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Nuc_PlaB1)

cy1_PlaB <- expressionLevelsPlaB$rawCounts$Cyt_PlaB1[,"cypt"]/(expressionLevelsPlaB$rawCounts$Cyt_PlaB1[,"cypt"]+expressionLevelsPlaB$rawCounts$Cyt_PlaB1[,"cymt"])
names(cy1_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Cyt_PlaB1)

py1_PlaB <- expressionLevelsPlaB$rawCounts$Poly_PlaB1[,"pypt"]/(expressionLevelsPlaB$rawCounts$Poly_PlaB1[,"pypt"]+expressionLevelsPlaB$rawCounts$Poly_PlaB1[,"pymt"])
names(py1_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Poly_PlaB1)

# Second replicate
ch2_PlaB <- expressionLevelsPlaB$rawCounts$Chr_PlaB2[,"chpt"]/(expressionLevelsPlaB$rawCounts$Chr_PlaB2[,"chpt"]+expressionLevelsPlaB$rawCounts$Chr_PlaB2[,"chmt"])
names(ch2_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Chr_PlaB2)

n2_PlaB <- expressionLevelsPlaB$rawCounts$Nuc_PlaB2[,"npt"]/(expressionLevelsPlaB$rawCounts$Nuc_PlaB2[,"npt"]+expressionLevelsPlaB$rawCounts$Nuc_PlaB2[,"nmt"])
names(n2_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Nuc_PlaB2)

cy2_PlaB <- expressionLevelsPlaB$rawCounts$Cyt_PlaB2[,"cypt"]/(expressionLevelsPlaB$rawCounts$Cyt_PlaB2[,"cypt"]+expressionLevelsPlaB$rawCounts$Cyt_PlaB2[,"cymt"])
names(cy2_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Cyt_PlaB2)

py2_PlaB <- expressionLevelsPlaB$rawCounts$Poly_PlaB2[,"pypt"]/(expressionLevelsPlaB$rawCounts$Poly_PlaB2[,"pypt"]+expressionLevelsPlaB$rawCounts$Poly_PlaB2[,"pymt"])
names(py2_PlaB) <- rownames(expressionLevelsPlaB$rawCounts$Poly_PlaB2)

chpn1_PlaB[chpn1_PlaB<=1e-10] <- NaN;nmn1_PlaB[nmn1_PlaB<=1e-10] <- NaN;npn2_PlaB[npn2_PlaB<=1e-10] <- NaN;pyn1_PlaB[pyn1_PlaB<=1e-10] <- NaN;ch1_PlaB[ch1_PlaB<=1e-10] <- NaN;py1_PlaB[py1_PlaB<=1e-10] <- NaN;cy2_PlaB[cy2_PlaB<=1e-10] <- NaN
chmn1_PlaB[chmn1_PlaB<=1e-10] <- NaN;chpn2_PlaB[chpn2_PlaB<=1e-10] <- NaN;nmn2_PlaB[nmn2_PlaB<=1e-10] <- NaN;cyn2_PlaB[cyn2_PlaB<=1e-10] <- NaN;n1_PlaB[n1_PlaB<=1e-10] <- NaN;ch2_PlaB[ch2_PlaB<=1e-10] <- NaN;py2_PlaB[py2_PlaB<=1e-10] <- NaN
npn1_PlaB[npn1_PlaB<=1e-10] <- NaN;chmn2_PlaB[chmn2_PlaB<=1e-10] <- NaN;cyn1_PlaB[cyn1_PlaB<=1e-10] <- NaN;pyn2_PlaB[pyn2_PlaB<=1e-10] <- NaN;cy1_PlaB[cy1_PlaB<=1e-10] <- NaN;n2_PlaB[n2_PlaB<=1e-10] <- NaN

## Premature Plot
commonGenes <- intersect(names(ch1),names(ch2))
ch <- apply(cbind(ch1[commonGenes],ch2[commonGenes]),1,mean,na.rm=TRUE)
ch <- ch[is.finite(ch)]

commonGenes_PlaB <- intersect(names(ch1_PlaB),names(ch2_PlaB))
ch_PlaB <- apply(cbind(ch1_PlaB[commonGenes_PlaB],ch2_PlaB[commonGenes_PlaB]),1,mean,na.rm=TRUE)
ch_PlaB <- ch_PlaB[is.finite(ch_PlaB)]

commonGenes <- intersect(names(ch),names(ch_PlaB))

### Figure 3C
smoothScatter(log2(ch[commonGenes]),log2(ch_PlaB[commonGenes]),xlab="Untreated",ylab="Pladienolide B",main="Log2 Proportion of\nPremature RNA reads",pch=".")
abline(0,1,lwd=2,col=1);abline(1,1,lwd=2,col=2);abline(2,1,lwd=2,col=3)

## Nascent
round(median(apply(cbind(chpn1_PlaB[fullModelGenes_PlaB],chpn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.45 vs 0.27
round(median(apply(cbind(npn1_PlaB[fullModelGenes_PlaB],npn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.42 vs 0.33
round(median(apply(cbind(chmn1_PlaB[fullModelGenes_PlaB],chmn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.33 vs 0.25
round(median(apply(cbind(nmn1_PlaB[fullModelGenes_PlaB],nmn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE),na.rm=TRUE),2) # 0.24 vs 0.17

ks.test(apply(cbind(chpn1_PlaB[fullModelGenes_PlaB],chpn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE),apply(cbind(chpn1[fullModelGenes],chpn2[fullModelGenes]),1,mean,na.rm=TRUE))$p.value
ks.test(apply(cbind(npn1_PlaB[fullModelGenes_PlaB],npn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE),apply(cbind(npn1[fullModelGenes],npn2[fullModelGenes]),1,mean,na.rm=TRUE))$p.value
ks.test(apply(cbind(chmn1_PlaB[fullModelGenes_PlaB],chmn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE),apply(cbind(chmn1[fullModelGenes],chmn2[fullModelGenes]),1,mean,na.rm=TRUE))$p.value
ks.test(apply(cbind(nmn1_PlaB[fullModelGenes_PlaB],nmn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE),apply(cbind(nmn1[fullModelGenes],nmn2[fullModelGenes]),1,mean,na.rm=TRUE))$p.value

### Figure S17
par(mfrow=c(1,2))
boxplot(list(Chp=apply(cbind(chpn1_PlaB[fullModelGenes_PlaB],chpn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE)
			,Np=apply(cbind(npn1_PlaB[fullModelGenes_PlaB],npn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE)
			,Chm=apply(cbind(chmn1_PlaB[fullModelGenes_PlaB],chmn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE)
			,Nm=apply(cbind(nmn1_PlaB[fullModelGenes_PlaB],nmn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE)
			,C=apply(cbind(cyn1_PlaB[fullModelGenes_PlaB],cyn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE)
			,P=apply(cbind(pyn1_PlaB[fullModelGenes_PlaB],pyn2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE))
		,outline=FALSE,ylab="Nascent RNA [%]",ylim=c(0,1),varwidth=TRUE,las=2)
boxplot(list(Ch=apply(cbind(ch1_PlaB[fullModelGenes_PlaB],ch2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE)
			,N=apply(cbind(n1_PlaB[fullModelGenes_PlaB],n2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE)
			,C=apply(cbind(cy1_PlaB[fullModelGenes_PlaB],cy2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE)
			,P=apply(cbind(py1_PlaB[fullModelGenes_PlaB],py2_PlaB[fullModelGenes_PlaB]),1,mean,na.rm=TRUE))
		,outline=FALSE,ylab="Premature RNA [%]",ylim=c(0,0.6),varwidth=TRUE,las=2)

