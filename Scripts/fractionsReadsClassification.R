#### Untreated
expressionLevels <- readRDS("/path/to/Results/Untreated/expressionLevelsMerged.rds")

### Nascent RNA
## First replicate
chpn1 <- expressionLevels$chromatin[[1]][,"chpn0.33"]/expressionLevels$chromatin[[1]][,"chpt"]
names(chpn1) <- rownames(expressionLevels$chromatin[[1]])
chmn1 <- expressionLevels$chromatin[[1]][,"chmn0.33"]/expressionLevels$chromatin[[1]][,"chmt"]
names(chmn1) <- rownames(expressionLevels$chromatin[[1]])

npn1 <- expressionLevels$nucleoplasmic[[1]][,"npn0.33"]/expressionLevels$nucleoplasmic[[1]][,"npt"]
names(npn1) <- rownames(expressionLevels$nucleoplasmic[[1]])
nmn1 <- expressionLevels$nucleoplasmic[[1]][,"nmn0.33"]/expressionLevels$nucleoplasmic[[1]][,"nmt"]
names(nmn1) <- rownames(expressionLevels$nucleoplasmic[[1]])

cyn1 <- expressionLevels$cytoplasmic[[1]][,"cyn0.33"]/expressionLevels$cytoplasmic[[1]][,"cyt"]
names(cyn1) <- rownames(expressionLevels$cytoplasmic[[1]])
pyn1 <- expressionLevels$polysomal[[1]][,"pyn0.33"]/expressionLevels$polysomal[[1]][,"pyt"]
names(pyn1) <- rownames(expressionLevels$polysomal[[1]])

### Premature RNA
ch1 <- expressionLevels$chromatin[[1]][,"chpt"]/(expressionLevels$chromatin[[1]][,"chpt"]+expressionLevels$chromatin[[1]][,"chmt"])
names(ch1) <- rownames(expressionLevels$chromatin[[1]])
n1 <- expressionLevels$nucleoplasmic[[1]][,"npt"]/(expressionLevels$nucleoplasmic[[1]][,"npt"]+expressionLevels$nucleoplasmic[[1]][,"nmt"])
names(n1) <- rownames(expressionLevels$nucleoplasmic[[1]])
cy1 <- expressionLevels$cytoplasmicP[[1]][,"cypt"]/(expressionLevels$cytoplasmicP[[1]][,"cypt"]+expressionLevels$cytoplasmicP[[1]][,"cymt"])
names(cy1) <- rownames(expressionLevels$cytoplasmicP[[1]])
py1 <- expressionLevels$polysomalP[[1]][,"pypt"]/(expressionLevels$polysomalP[[1]][,"pypt"]+expressionLevels$polysomalP[[1]][,"pymt"])
names(py1) <- rownames(expressionLevels$polysomalP[[1]])

### Significance
wilcox.test(cy1,py1)$p.value
 # [1] 2.403187e-38

### Plots - Figure S10
par(mfrow=c(1,2))
boxplot(list(Chp=chpn1[expressionLevels$expressedGenesYesChpNpP[[1]]],Np=npn1[expressionLevels$expressedGenesYesChpNpP[[1]]],Chm=chmn1[expressionLevels$expressedGenesYesChpNpP[[1]]],Nm=nmn1[expressionLevels$expressedGenesYesChpNpP[[1]]],C=cyn1[expressionLevels$expressedGenesYesChpNpP[[1]]],P=pyn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),outline=FALSE,ylab="Nascent RNA [%]",ylim=c(0,1),varwidth=TRUE,las=2)
boxplot(list(Ch=ch1[expressionLevels$expressedGenesYesChpNpP[[1]]],N=n1[expressionLevels$expressedGenesYesChpNpP[[1]]],C=cy1[expressionLevels$expressedGenesYesChpNpP[[1]]],P=py1[expressionLevels$expressedGenesYesChpNpP[[1]]]),outline=FALSE,ylab="Premature RNA [%]",ylim=c(0,0.6),varwidth=TRUE,las=2)

### Numerical values
## Premature
round(median(chpn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.39
round(median(chmn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.32
round(median(npn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.33
round(median(nmn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.24
round(median(cyn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.25
round(median(pyn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.29

round(median(ch1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.15
round(median(n1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.07
round(median(cy1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.06
round(median(py1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.07

#### Pladienolide B
expressionLevels <- readRDS("/path/to/Results/PladienolideB/expressionLevelsMerged.rds")

### Nascent RNA
## First replicate
chpn1 <- expressionLevels$chromatin[[1]][,"chpn0.33"]/expressionLevels$chromatin[[1]][,"chpt"]
names(chpn1) <- rownames(expressionLevels$chromatin[[1]])
chmn1 <- expressionLevels$chromatin[[1]][,"chmn0.33"]/expressionLevels$chromatin[[1]][,"chmt"]
names(chmn1) <- rownames(expressionLevels$chromatin[[1]])

npn1 <- expressionLevels$nucleoplasmic[[1]][,"npn0.33"]/expressionLevels$nucleoplasmic[[1]][,"npt"]
names(npn1) <- rownames(expressionLevels$nucleoplasmic[[1]])
nmn1 <- expressionLevels$nucleoplasmic[[1]][,"nmn0.33"]/expressionLevels$nucleoplasmic[[1]][,"nmt"]
names(nmn1) <- rownames(expressionLevels$nucleoplasmic[[1]])

cyn1 <- expressionLevels$cytoplasmic[[1]][,"cyn0.33"]/expressionLevels$cytoplasmic[[1]][,"cyt"]
names(cyn1) <- rownames(expressionLevels$cytoplasmic[[1]])
pyn1 <- expressionLevels$polysomal[[1]][,"pyn0.33"]/expressionLevels$polysomal[[1]][,"pyt"]
names(pyn1) <- rownames(expressionLevels$polysomal[[1]])

### Premature RNA
ch1 <- expressionLevels$chromatin[[1]][,"chpt"]/(expressionLevels$chromatin[[1]][,"chpt"]+expressionLevels$chromatin[[1]][,"chmt"])
names(ch1) <- rownames(expressionLevels$chromatin[[1]])
n1 <- expressionLevels$nucleoplasmic[[1]][,"npt"]/(expressionLevels$nucleoplasmic[[1]][,"npt"]+expressionLevels$nucleoplasmic[[1]][,"nmt"])
names(n1) <- rownames(expressionLevels$nucleoplasmic[[1]])
cy1 <- expressionLevels$cytoplasmicP[[1]][,"cypt"]/(expressionLevels$cytoplasmicP[[1]][,"cypt"]+expressionLevels$cytoplasmicP[[1]][,"cymt"])
names(cy1) <- rownames(expressionLevels$cytoplasmicP[[1]])
py1 <- expressionLevels$polysomalP[[1]][,"pypt"]/(expressionLevels$polysomalP[[1]][,"pypt"]+expressionLevels$polysomalP[[1]][,"pymt"])
names(py1) <- rownames(expressionLevels$polysomalP[[1]])

### Plots - Figure S16
pdf("fractionsClassificationPlaBMerged_fullModelGenes.pdf",width=5,height=4)
par(mfrow=c(1,2))
boxplot(list(Chp=chpn1[expressionLevels$expressedGenesYesChpNpP[[1]]],Np=npn1[expressionLevels$expressedGenesYesChpNpP[[1]]],Chm=chmn1[expressionLevels$expressedGenesYesChpNpP[[1]]],Nm=nmn1[expressionLevels$expressedGenesYesChpNpP[[1]]],C=cyn1[expressionLevels$expressedGenesYesChpNpP[[1]]],P=pyn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),outline=FALSE,ylab="Nascent RNA [%]",ylim=c(0,1),varwidth=TRUE,las=2)
boxplot(list(Ch=ch1[expressionLevels$expressedGenesYesChpNpP[[1]]],N=n1[expressionLevels$expressedGenesYesChpNpP[[1]]],C=cy1[expressionLevels$expressedGenesYesChpNpP[[1]]],P=py1[expressionLevels$expressedGenesYesChpNpP[[1]]]),outline=FALSE,ylab="Premature RNA [%]",ylim=c(0,0.6),varwidth=TRUE,las=2)
dev.off()

### Numerical values
## Premature
round(median(chpn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.22
round(median(chmn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.16
round(median(npn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.25
round(median(nmn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.16
round(median(cyn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.21
round(median(pyn1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.33

round(median(ch1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.25
round(median(n1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.10
round(median(cy1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.09
round(median(py1[expressionLevels$expressedGenesYesChpNpP[[1]]]),2) # 0.10

#### PlaB vs Untreated - Figure 3B
## Chp
expressionLevelsUntreated <- readRDS("/path/to/Results/Untreated/expressionLevelsMerged.rds")
expressionLevelsPlaB <- readRDS("/path/to/Results/PladienolideB/expressionLevelsMerged.rds")

chUntreated <- expressionLevelsUntreated$chromatin[[1]][,"chpt"]/(expressionLevelsUntreated$chromatin[[1]][,"chpt"]+expressionLevelsUntreated$chromatin[[1]][,"chmt"])
names(chUntreated) <- rownames(expressionLevelsUntreated$chromatin[[1]])

chPlaB <- expressionLevelsPlaB$chromatin[[1]][,"chpt"]/(expressionLevelsPlaB$chromatin[[1]][,"chpt"]+expressionLevelsPlaB$chromatin[[1]][,"chmt"])
names(chPlaB) <- rownames(expressionLevelsPlaB$chromatin[[1]])

ks.test(chUntreated[expressionLevelsUntreated$expressedGenesYesChpNpP[[1]]],chPlaB[expressionLevelsPlaB$expressedGenesYesChpNpP[[1]]])$p.value # 0
wilcox.test(chUntreated[expressionLevelsUntreated$expressedGenesYesChpNpP[[1]]],chPlaB[expressionLevelsPlaB$expressedGenesYesChpNpP[[1]]],alternative="less")$p.value # 5.467804e-68

allGenes <- intersect(names(chUntreated),names(chPlaB))
allGenesFull <- intersect(expressionLevelsUntreated$expressedGenesYesChpNpP[[1]],expressionLevelsPlaB$expressedGenesYesChpNpP[[1]])

smoothScatter(log2(chUntreated[allGenes]),log2(chPlaB[allGenes]),xlab="Untreated",ylab="Pladienolide B",main="Log2 Proportion of\nPremature RNA reads",pch=".")
abline(0,1,lwd=2,col=1);abline(1,1,lwd=2,col=2);abline(2,1,lwd=2,col=3)

## Nuclear nascent
chpnUntreated<- expressionLevelsUntreated$chromatin[[1]][,"chpn0.33"]/expressionLevelsUntreated$chromatin[[1]][,"chpt"]
names(chpnUntreated) <- rownames(expressionLevelsUntreated$chromatin[[1]])
chmnUntreated<- expressionLevelsUntreated$chromatin[[1]][,"chmn0.33"]/expressionLevelsUntreated$chromatin[[1]][,"chmt"]
names(chmnUntreated) <- rownames(expressionLevelsUntreated$chromatin[[1]])
npnUntreated<- expressionLevelsUntreated$nucleoplasmic[[1]][,"npn0.33"]/expressionLevelsUntreated$nucleoplasmic[[1]][,"npt"]
names(npnUntreated) <- rownames(expressionLevelsUntreated$nucleoplasmic[[1]])
nmnUntreated<- expressionLevelsUntreated$nucleoplasmic[[1]][,"nmn0.33"]/expressionLevelsUntreated$nucleoplasmic[[1]][,"nmt"]
names(nmnUntreated) <- rownames(expressionLevelsUntreated$nucleoplasmic[[1]])

chpnPlaB<- expressionLevelsPlaB$chromatin[[1]][,"chpn0.33"]/expressionLevelsPlaB$chromatin[[1]][,"chpt"]
names(chpnPlaB) <- rownames(expressionLevelsPlaB$chromatin[[1]])
chmnPlaB<- expressionLevelsPlaB$chromatin[[1]][,"chmn0.33"]/expressionLevelsPlaB$chromatin[[1]][,"chmt"]
names(chmnPlaB) <- rownames(expressionLevelsPlaB$chromatin[[1]])
npnPlaB<- expressionLevelsPlaB$nucleoplasmic[[1]][,"npn0.33"]/expressionLevelsPlaB$nucleoplasmic[[1]][,"npt"]
names(npnPlaB) <- rownames(expressionLevelsPlaB$nucleoplasmic[[1]])
nmnPlaB<- expressionLevelsPlaB$nucleoplasmic[[1]][,"nmn0.33"]/expressionLevelsPlaB$nucleoplasmic[[1]][,"nmt"]
names(nmnPlaB) <- rownames(expressionLevelsPlaB$nucleoplasmic[[1]])

ks.test(chpnUntreated[expressionLevelsUntreated$expressedGenesYesChpNpP[[1]]],chpnPlaB[expressionLevelsPlaB$expressedGenesYesChpNpP[[1]]])$p.value # 0
ks.test(chmnUntreated[expressionLevelsUntreated$expressedGenesYesChpNpP[[1]]],chmnPlaB[expressionLevelsPlaB$expressedGenesYesChpNpP[[1]]])$p.value # 0
ks.test(npnUntreated[expressionLevelsUntreated$expressedGenesYesChpNpP[[1]]],npnPlaB[expressionLevelsPlaB$expressedGenesYesChpNpP[[1]]])$p.value # 0
ks.test(nmnUntreated[expressionLevelsUntreated$expressedGenesYesChpNpP[[1]]],nmnPlaB[expressionLevelsPlaB$expressedGenesYesChpNpP[[1]]])$p.value # 0