# Functions that estimates CVs at gene-level
CVEstimate <- function(Fraction){
	sapply(colnames(Fraction[[1]]),function(i){
		vec <- cbind(Fraction[[1]][,i],Fraction[[2]][,i])

		m <- apply(vec,1,mean)
		std <- apply(vec,1,sd)
		CV <- std/m
	})
}

# Read untreated expression levels
WT <- readRDS("/path/to/Results/Untreated/expressionLevels.rds")
#Picking common genes between replicates
CVall <- list(WT$normalizedCountsStatistics$WT1$mean,WT$normalizedCountsStatistics$WT2$mean)
CVall <- lapply(CVall,function(i){i[intersect(rownames(CVall[[1]]),rownames(CVall[[2]])),]})

#Dividing CVs per fraction
Chromatin <- list(CVall[[1]][,grep("ch",colnames(CVall[[1]]))],
				  CVall[[2]][,grep("ch",colnames(CVall[[2]]))])

Nucleoplasmic <- list(CVall[[1]][,grep("^n",colnames(CVall[[1]]))],
				  CVall[[2]][,grep("^n",colnames(CVall[[2]]))])

Cytoplasmic <- list(CVall[[1]][,grep("cyp",colnames(CVall[[1]]))]+CVall[[1]][,grep("cym",colnames(CVall[[1]]))],
					CVall[[2]][,grep("cyp",colnames(CVall[[2]]))]+CVall[[2]][,grep("cym",colnames(CVall[[2]]))])
Cytoplasmic <- lapply(Cytoplasmic,function(i){colnames(i) <- gsub("cyp","cy",colnames(i));i})

Polysomal <- list(CVall[[1]][,grep("pyp",colnames(CVall[[1]]))]+CVall[[1]][,grep("pym",colnames(CVall[[1]]))],
					CVall[[2]][,grep("pyp",colnames(CVall[[2]]))]+CVall[[2]][,grep("pym",colnames(CVall[[2]]))])
Polysomal <- lapply(Polysomal,function(i){colnames(i) <- gsub("pyp","py",colnames(i));i})

# Adding total 
Chromatin <- lapply(Chromatin,function(i)
{
	i <- cbind(i
			  ,chp0.33=rowSums(i[,c("chpp0.33","chmp0.33")])
			  ,chn0.33=rowSums(i[,c("chpn0.33","chmn0.33")])
			  ,chp0=rowSums(i[,c("chpp0","chmp0")])
			  ,chn0=rowSums(i[,c("chpn0","chmn0")])
			  ,cht=rowSums(i[,c("chpt","chmt")]))
})
Nucleoplasmic <- lapply(Nucleoplasmic,function(i)
{
	i <- cbind(i
			  ,np0.33=rowSums(i[,c("npp0.33","nmp0.33")])
			  ,nn0.33=rowSums(i[,c("npn0.33","nmn0.33")])
			  ,np0=rowSums(i[,c("npp0","nmp0")])
			  ,nn0=rowSums(i[,c("npn0","nmn0")])
			  ,nt=rowSums(i[,c("npt","nmt")]))
})

#Computing CVs for each gene
CV_chr <- CVEstimate(Chromatin)
CV_nucleo <- CVEstimate(Nucleoplasmic)
CV_cyto <- CVEstimate(Cytoplasmic)
CV_poly <- CVEstimate(Polysomal)

### Boxplot
renameFunctionTmp <- function(mat)
{
	  colnames(mat) <- gsub("chpp","Premature_PreExisting_",colnames(mat))
  	  colnames(mat) <- gsub("chpn","Premature_Nascent_",colnames(mat))
  	  colnames(mat) <- gsub("chpt","Premature_Total",colnames(mat))
  	  colnames(mat) <- gsub("chmp","Mature_PreExisting_",colnames(mat))
  	  colnames(mat) <- gsub("chmn","Mature_Nascent_",colnames(mat))
  	  colnames(mat) <- gsub("chmt","Mature_Total",colnames(mat))
  	  colnames(mat) <- gsub("npp","Premature_PreExisting_",colnames(mat))
  	  colnames(mat) <- gsub("npn","Premature_Nascent_",colnames(mat))
  	  colnames(mat) <- gsub("npt","Premature_Total",colnames(mat))
  	  colnames(mat) <- gsub("nmp","Mature_PreExisting_",colnames(mat))
  	  colnames(mat) <- gsub("nmn","Mature_Nascent_",colnames(mat))
  	  colnames(mat) <- gsub("nmt","Mature_Total",colnames(mat))
  	  colnames(mat) <- gsub("cyp","PreExisting_",colnames(mat))
  	  colnames(mat) <- gsub("cyn","Nascent_",colnames(mat))
  	  colnames(mat) <- gsub("cyt","Total",colnames(mat))
  	  colnames(mat) <- gsub("pyp","PreExisting_",colnames(mat))
  	  colnames(mat) <- gsub("pyn","Nascent_",colnames(mat))
  	  colnames(mat) <- gsub("pyt","Total",colnames(mat))
	  colnames(mat) <- gsub("chp","PreExisting_",colnames(mat))
	  colnames(mat) <- gsub("chn","Nascent_",colnames(mat))
	  colnames(mat) <- gsub("cht","Total",colnames(mat))
	  colnames(mat) <- gsub("np","PreExisting_",colnames(mat))
	  colnames(mat) <- gsub("nn","Nascent_",colnames(mat))
	  colnames(mat) <- gsub("nt","Total",colnames(mat))
  	  colnames(mat) <- gsub("0.33","20min",colnames(mat))
  	  mat
}

par(mfrow=c(2,4))
boxplot(renameFunctionTmp(CV_chr[,!grepl("n0$",colnames(CV_chr))]),main="Chromatin",outline=FALSE,notch=TRUE,las=2,ylab="CV")
boxplot(renameFunctionTmp(CV_nucleo[,!grepl("n0$",colnames(CV_nucleo))]),main="Nucleoplasmic",outline=FALSE,notch=TRUE,las=2,ylab="CV")
boxplot(renameFunctionTmp(CV_cyto[,!grepl("n0$",colnames(CV_cyto))]),main="Cytoplasmic",outline=FALSE,notch=TRUE,las=2,ylab="CV")
boxplot(renameFunctionTmp(CV_poly[,!grepl("n0$",colnames(CV_poly))]),main="Polysomal",outline=FALSE,notch=TRUE,las=2,ylab="CV")

# Computing mean CVs
Cv_chr <-  apply(CV_chr,2,median,na.rm=TRUE)
Cv_nucleo <-  apply(CV_nucleo,2,median,na.rm=TRUE)
Cv_cyto <-  apply(CV_cyto,2,median,na.rm=TRUE)
Cv_poly <-  apply(CV_poly,2,median,na.rm=TRUE)

### Setting not defined CVs to 0.5 which is higher than the median CV of each fraction.
Cv_chr[is.na(Cv_chr)] <- 0.5
Cv_nucleo[is.na(Cv_nucleo)] <- 0.5
Cv_cyto[is.na(Cv_cyto)] <- 0.5
Cv_poly[is.na(Cv_poly)] <- 0.5

CVs <- list("Chromatin"=Cv_chr
			,"Nucleoplasmic"=Cv_nucleo
			,"Cytoplasmic"=Cv_cyto
			,"Polysomal"=Cv_poly)
saveRDS(file="CVs.rds",Cvs)
