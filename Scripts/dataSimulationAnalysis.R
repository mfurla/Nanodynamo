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
library("tidyr")
library("RColorBrewer")

source("/path/to/allInternalFunctions.R")

### Data loading
## Number of replicate
# Simulated datasets
simulatedDataset1 <- readRDS("path/to/Results/dataSimulation/FullModel_Replicates/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
simulatedDataset2 <- readRDS("path/to/Results/dataSimulation/FullModel_Replicates/simulatedDataset_Noise_TRUE_Reps_2_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
simulatedDataset3 <- readRDS("path/to/Results/dataSimulation/FullModel_Replicates/simulatedDataset_Noise_TRUE_Reps_3_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
simulatedDataset4 <- readRDS("path/to/Results/dataSimulation/FullModel_Replicates/simulatedDataset_Noise_TRUE_Reps_4_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
simulatedDataset5 <- readRDS("path/to/Results/dataSimulation/FullModel_Replicates/simulatedDataset_Noise_TRUE_Reps_5_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")

# Inferred rates
inferedNumericalModel1 <- readRDS("path/to/Results/dataSimulation/FullModel_Replicates/inferedNumericalModel_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
inferedNumericalModel2 <- readRDS("path/to/Results/dataSimulation/FullModel_Replicates/inferedNumericalModel_Noise_TRUE_Reps_2_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
inferedNumericalModel3 <- readRDS("path/to/Results/dataSimulation/FullModel_Replicates/inferedNumericalModel_Noise_TRUE_Reps_3_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
inferedNumericalModel4 <- readRDS("path/to/Results/dataSimulation/FullModel_Replicates/inferedNumericalModel_Noise_TRUE_Reps_4_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
inferedNumericalModel5 <- readRDS("path/to/Results/dataSimulation/FullModel_Replicates/inferedNumericalModel_Noise_TRUE_Reps_5_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")

## Temporal design
# Simulated datasets
simulatedDataset1tp <- readRDS("path/to/Results/dataSimulation/FullModel_Replicates/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
simulatedDataset2tp <- readRDS("path/to/Results/dataSimulation/FullModel_TemporalDesign/TimePoints_0_0.33_1/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_1_TauPoly_0_0.33_1_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
simulatedDataset3tp <- readRDS("path/to/Results/dataSimulation/FullModel_TemporalDesign/TimePoints_0_0.33_1_2/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_1_2_TauPoly_0_0.33_1_2_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")

# Inferred rates
inferedNumericalModel1tp <- readRDS("path/to/Results/dataSimulation/FullModel_Replicates/inferedNumericalModel_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
inferedNumericalModel2tp <- readRDS("path/to/Results/dataSimulation/FullModel_TemporalDesign/TimePoints_0_0.33_1/inferedNumericalModel_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_1_TauPoly_0_0.33_1_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
inferedNumericalModel3tp <- readRDS("path/to/Results/dataSimulation/FullModel_TemporalDesign/TimePoints_0_0.33_1_2/inferedNumericalModel_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_1_2_TauPoly_0_0.33_1_2_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")

## Models with nucleoplasmic RNA decay
# Simulated datasets
simulatedDatasetPND <- readRDS("path/to/Results/dataSimulation/PrematureNucleoplasmicDegradation/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3_k11_0.65.rds")
simulatedDatasetMND <- readRDS("path/to/Results/dataSimulation/MatureNucleoplasmicDegradation/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k9_0.65_k10_1.3.rds")

# Inferred rates
inferedNumericalPND <- readRDS("path/to/Results/dataSimulation/PrematureNucleoplasmicDegradation/inferedNumericalModel_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k10_1.3_k11_0.65.rds")
inferedNumericalMND <- readRDS("path/to/Results/dataSimulation/MatureNucleoplasmicDegradation/inferedNumericalModel_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k7_0.65_k8_1.3_k9_0.65_k10_1.3.rds")

### Correlations
## Number of replicates
replicatesCorrelations <- sapply(colnames(simulatedDataset1$exampleRates),function(i)
{
    c("1rep"=cor(simulatedDataset1$exampleRates[,i],inferedNumericalModel1$inferedRates[,i],method="s")
     ,"2rep"=cor(simulatedDataset2$exampleRates[,i],inferedNumericalModel2$inferedRates[,i],method="s")
     ,"3rep"=cor(simulatedDataset3$exampleRates[,i],inferedNumericalModel3$inferedRates[,i],method="s")
     ,"4rep"=cor(simulatedDataset4$exampleRates[,i],inferedNumericalModel4$inferedRates[,i],method="s")
     ,"5rep"=cor(simulatedDataset5$exampleRates[,i],inferedNumericalModel5$inferedRates[,i],method="s"))
})

## Temporal design
timePointsCorrelations <- sapply(colnames(simulatedDataset1$exampleRates),function(i)
{
    c("20 min"=cor(simulatedDataset1tp$exampleRates[,i],inferedNumericalModel1tp$inferedRates[,i],method="s")
     ,"20,60 min"=cor(simulatedDataset2tp$exampleRates[,i],inferedNumericalModel2tp$inferedRates[,i],method="s")
     ,"20,60,120 min"=cor(simulatedDataset3tp$exampleRates[,i],inferedNumericalModel3tp$inferedRates[,i],method="s"))
})

x <- rbind(replicatesCorrelations,rep(NaN,nrow(replicatesCorrelations)),timePointsCorrelations)
x <- x[,c("k1","k2","k3","k4","k5","k6","k8","k7","k10")]
colnames(x) <- paste0("k",1:9)

### Plots
## Replicates - Figure 1C
myPal <- brewer.pal(8,"Blues")

files <- list(inferedNumericalModel1$ratesCorrelations,
              inferedNumericalModel2$ratesCorrelations,
              inferedNumericalModel3$ratesCorrelations,
              inferedNumericalModel4$ratesCorrelations,
              inferedNumericalModel5$ratesCorrelations)
cor <- unname(data.frame(matrix(unlist(files),nrow=length(files)*ncol(files[[1]]),ncol=1)))
n <- c("k1","k2","k3","k4","k5","k6","k7","k8","k9")
r <- unlist(lapply(seq_along(c(1:5)),function(j){rep(j,9)}))
cor <- cbind("Reps" = as.factor(r),"Rates" = as.factor(rep(n,length(files))),"Correlations"=cor)

pl <- ggplot(cor,aes(x = Rates,y = Correlations,fill = Reps)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=myPal[3:7]) +
 theme(legend.position="top") + theme( axis.title = element_text(size = 12), axis.text = element_text(size = 8)) + 
  ylab("Spearman Correlation") + 
  labs(text= element_text("N° of reps")) 

## Labeling pulses - Figure 1D
files <- list(inferedNumericalModel1tp$ratesCorrelations,
              inferedNumericalModel2tp$ratesCorrelations,
              inferedNumericalModel3tp$ratesCorrelations)
cor <- unname(data.frame(matrix(unlist(files),nrow=length(files)*ncol(files[[1]]),ncol=1)))
n <- c("k1","k2","k3","k4","k5","k6","k7","k8","k9")
r <- unlist(lapply(c("0.33","0.33_1","0.33_1_2"),function(j){rep(j,9)}))
cor <- cbind("Pulses" = as.factor(r),"Rates" = as.factor(rep(n,3)),"Correlations"=cor)

pl <- ggplot(cor,aes(x = Rates,y = Correlations,fill = Pulses)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=myPal[4:6]) +
  theme(legend.position="top") + theme( axis.title = element_text(size = 12), axis.text = element_text(size = 8)) + 
  ylab("Spearman Correlation") + 
  labs(text= element_text("Labeling pulses")) 

# Expected vs Inferred data - Figure 1E
expected <- inferedNumericalModel1$expressionData
expected <- expected[,-grep("0$",colnames(expected))]
inferred <- inferedNumericalModel1$inferedData
inferred <- inferred[,-grep("0$",colnames(inferred))]
cor <- diag(round(cor(inferred,expected,method="s",use="c"),2))
times <- c(rep("Nascent",6),rep("Pre-existing",6),rep("Total",6))
Species <- rep(c("Chp","Chm","Np","Nm","C","P"),3)
ord <- c("chp","chm","np","nm","c","p")
ti <- c("n0.33","p0.33","t")
ord <- matrix(sapply(ord,function(i){sapply(ti,function(j){paste0(i,j,collapse="_")})}),nrow=1)
cor <- cor[order(match(names(cor),matrix(sapply(ord,function(i){sapply(ti,function(j){paste0(i,j,collapse="_")})}),nrow=1)))]
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

## Figure S4
foe <- unlist(GoodnessOfFit(inferedNumericalModel1tp,width=7,height=12,name="Sim",lowSat=0.025,upSat=0.975))

## Figure 1F
foe <- unlist(simulatedRatesCorrelations(object1=inferedNumericalModel1tp,simObject=simulatedDataset1tp,width=7,height=12,name="",lowSat=0.025,upSat=0.975))

## Figure S12A
RatesDistributions(object=inferedNumericalModel1tp,width=25,height=4,obj_name="Sim",xlimTmp=c(-2,6))

### Simpler models - Supplemental Material Figures
foe <- unlist(simulatedRatesCorrelations(object1=readRDS("/path/to/Results/dataSimulation/FullModel_NoPolysomal/inferedNumericalModel_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k8_1.3.rds")
                    ,simObject=readRDS("/path/to/Results/dataSimulation/FullModel_NoPolysomal/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k4_3_k5_7.5_k6_45_k8_1.3.rds")
                    ,width=7,height=12,name="FullModel_NoPolysomal",lowSat=0.025,upSat=0.975))
round(c(min(foe),max(foe),median(foe)),2)
# [1] 0.86 0.95 0.90
foe <- unlist(simulatedRatesCorrelations(object1=readRDS("/path/to/Results/dataSimulation/NoNucleoplasmicPremature/inferedNumericalModel_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
                    ,simObject=readRDS("/path/to/Results/dataSimulation/NoNucleoplasmicPremature/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
                    ,width=7,height=12,name="NoNucleoplasmicPremature",lowSat=0.025,upSat=0.975))
round(c(min(foe),max(foe),median(foe)),2)
# [1] 0.63 0.94 0.90
foe <- unlist(simulatedRatesCorrelations(object1=readRDS("/path/to/Results/dataSimulation/NoNucleoplasmicPremature_NoPolysomal/inferedNumericalModel_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k6_45_k8_1.3.rds")
                    ,simObject=readRDS("/path/to/Results/dataSimulation/NoNucleoplasmicPremature_NoPolysomal/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k2_30_k3_30_k6_45_k8_1.3.rds")
                    ,width=7,height=12,name="NoNucleoplasmicPremature_NoPolysomal",lowSat=0.025,upSat=0.975))
round(c(min(foe),max(foe),median(foe)),2)
# [1] 0.92 0.96 0.95
foe <- unlist(simulatedRatesCorrelations(object1=readRDS("/path/to/Results/dataSimulation/NoPremature/inferedNumericalModel_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k3_30_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
                    ,simObject=readRDS("/path/to/Results/dataSimulation/NoPremature/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k3_30_k6_45_k7_0.65_k8_1.3_k10_1.3.rds")
                    ,width=7,height=12,name="NoP",lowSat=0.025,upSat=0.975))
round(c(min(foe),max(foe),median(foe)),2)
# [1] 0.60 0.93 0.84
foe <- unlist(simulatedRatesCorrelations(object1=readRDS("/path/to/Results/dataSimulation/NoPremature_NoPolysomal/inferedNumericalModel_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k3_30_k6_45_k8_1.3.rds")
                    ,simObject=readRDS("/path/to/Results/dataSimulation/NoPremature_NoPolysomal/simulatedDataset_Noise_TRUE_Reps_1_CV_0.35_TauFractions_0_0.33_TauPoly_0_0.33_TauTotal__k1_12_k3_30_k6_45_k8_1.3.rds")
                    ,width=7,height=12,name="NoPremature_NoPolysomal",lowSat=0.025,upSat=0.975))
round(c(min(foe),max(foe),median(foe)),2)
# [1] 0.90 0.95 0.93

### Nuclear RNA decay - Figures S30, S31
foe <- unlist(simulatedRatesCorrelations(object1=inferedNumericalPND,simObject=simulatedDatasetPND,width=7,height=12,name="",lowSat=0.025,upSat=0.975))
round(c(min(foe),max(foe),median(foe)),2)

foe <- unlist(simulatedRatesCorrelations(object1=inferedNumericalMND,simObject=simulatedDatasetMND,width=7,height=12,name="",lowSat=0.025,upSat=0.975))
round(c(min(foe),max(foe),median(foe)),2)