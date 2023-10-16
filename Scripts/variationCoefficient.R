### Libraries
library(ggplot2)

### Data loading - from the supplemental material of de Pretis et al. Genome Research 2017
load("path/to/3T9MycER_rpkms.RData")

### CV profiling
## Total RNA
totalRNAMeanTmp <- apply(totalRpkmsExons,1,mean)
totalRNASdTmp <- apply(totalRpkmsExons,1,sd)
totalRNACvTmp <- totalRNASdTmp/totalRNAMeanTmp

## Premature RNA
prematureRNAMeanTmp <- apply(totalRpkmsIntrons,1,mean)
prematureRNASdTmp <- apply(totalRpkmsIntrons,1,sd)
prematureRNACvTmp <- prematureRNASdTmp/prematureRNAMeanTmp

## Total Nascent RNA
totalNascentRNAMeanTmp <- apply(labeledRpkmsExons,1,mean)
totalNascentRNASdTmp <- apply(labeledRpkmsExons,1,sd)
totalNascentRNACvTmp <- totalNascentRNASdTmp/totalNascentRNAMeanTmp

## Premature Nascent RNA
prematureNascentRNAMeanTmp <- apply(labeledRpkmsIntrons,1,mean)
prematureNascentRNASdTmp <- apply(labeledRpkmsIntrons,1,sd)
prematureNascentRNACvTmp <- prematureNascentRNASdTmp/prematureNascentRNAMeanTmp

### Plot - Figure S2
pdf("CV_distribution.pdf",width=4,height=4)
ggplot(data.frame(CV=c(totalRNACvTmp,prematureRNACvTmp,totalNascentRNACvTmp,prematureNascentRNACvTmp)
				 ,RNA_species=rep(c("Total_RNA","Premature_RNA","Total_Nascent_RNA","Premature_Nascent_RNA")
				 ,each=length(totalRNACvTmp))), aes(x = CV, fill = RNA_species)) + geom_density(alpha = 0.5) + theme(legend.position = c(0.65, 0.65))
dev.off()