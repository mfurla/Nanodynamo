inferRates <- function(expressionData # Expression data of the genes to be modeled.
                      ,expressionDataDev # Standard deviations of the genes to be modeled.
                      ,simulatedDataset=NULL # Simulated dataset if this is the case (just to produce real rates correlations).
                      ,initialRates # List of initial rates for optimization.
                      ,TauFractions # Time points with cellular fractionation.
                      ,TauPoly # Time points with polysomal profiling.
                      ,TauTotal # Time points with total nascent and pre-existing RNA profiling.
                      ,cpus # Number of cpus.
                      ,logOptim=TRUE # TRUE to optimize the model parameters in the Log space.
                      ,lowB=1e-6 # Lower boundary for the rates.
                      ,upB=1e4 # Upper boundary for the rates.
                      ,FlagDev="FC" # Cost function.
                      ,lambda=0.05 # Regularization strength.
                      ,excludeSpecies=NULL # List of species to be excluded from the cost function.
                      ,parFixed=NULL) # List of parameters to be excluded from the optimization.
{
  # Definition of the dataset as simulated or not.
  if(!is.null(simulatedDataset))
  {
    simulatedDatasetFlag <- TRUE
    initialRates <- lapply(initialRates,function(i)i[colnames(simulatedDataset$exampleRates)])
  }else{
    simulatedDatasetFlag <- FALSE
  }

  # Sorting of expression levels.
  colOrderTmp <- sort(colnames(expressionData))
  expressionData <- expressionData[,colOrderTmp]
  expressionDataDev <- expressionDataDev[,colOrderTmp]

  # If standard deviation are missing assume a CV of 1.
  if(is.null(expressionDataDev)){expressionDataDev<-expressionData}
  
  # If the optimization in the logarithmic space is required transform the initial rates.
  if(logOptim){initialRates <- lapply(initialRates,log)}

  # Selection of the configuration to be simulated according to the provided example rates or the RNA species.
  # This means definition of the modelSimulation function (ode implementation), of dataGeneration function (ode solution),
  # and in case of polysomal RNA profiling exclusion of nascent RNA from the species for the cost function optimization.

  if("k9"%in%names(unlist(initialRates))
     &any(grepl("^chpp",colnames(expressionData)))
     &any(grepl("^npp",colnames(expressionData)))
     &any(grepl("^py",colnames(expressionData))))
  {
    dataGeneration <- dataGenerationFullND
    print("Inference mode: Nuclear decay.")
  }else if("k11"%in%names(unlist(initialRates))
     &any(grepl("^chpp",colnames(expressionData)))
     &any(grepl("^npp",colnames(expressionData)))
     &any(grepl("^py",colnames(expressionData))))
  {
    dataGeneration <- dataGenerationFullND_prem
    print("Inference mode: Nuclear decay prem.")
  }else if(any(grepl("^chpp",colnames(expressionData)))
          &any(grepl("^npp",colnames(expressionData)))
          &!any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: No polysomal RNA, Chromatin premature RNA, Nucleoplasmic premature RNA.")
    dataGeneration <- dataGenerationNoPoly
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k7","k9","k10")]))
  }else if(any(grepl("^chpp",colnames(expressionData)))
          &any(grepl("^npp",colnames(expressionData)))
          &any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: Full model")
    dataGeneration <- dataGenerationFull
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k9")]))
  }else if(any(grepl("^chpp",colnames(expressionData)))
         &!any(grepl("^npp",colnames(expressionData)))
         &!any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: No polysomal RNA, Chromatin premature RNA.")
    dataGeneration <- dataGenerationNoPolyNoNucP
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k4","k5","k7","k9","k10")]))
  }else if(any(grepl("^chpp",colnames(expressionData)))
         &!any(grepl("^npp",colnames(expressionData)))
         &any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: Yes polysomal RNA, Chromatin premature RNA.")
    dataGeneration <- dataGenerationNoNucP
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k4","k5","k9")]))
  }else if(!any(grepl("^chpp",colnames(expressionData)))
          &!any(grepl("^npp",colnames(expressionData)))
          &!any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: No polysomal RNA, No premature RNA.")
    dataGeneration <- dataGenerationNoPolyNoP
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k2","k4","k5","k7","k9","k10")]))
  }else if(!any(grepl("^chpp",colnames(expressionData)))
          &!any(grepl("^npp",colnames(expressionData)))
          &any(grepl("^py",colnames(expressionData))))
  {
    print("Inference mode: Yes polysomal RNA, No premature RNA.")
    dataGeneration <- dataGenerationNoP
    initialRates <- lapply(initialRates,function(i)c(i[!names(i)%in%c("k2","k4","k5","k9")]))
  }else{
    print("Issues in the input data!")
    return(NULL)
  }

  # Inference loop for any gene (row) in expressionData.
  inferedRates <- t(mcsapply(1:nrow(expressionData),function(i)
  {
    print(i)

    # Expression levels isolation.
    exampleDataTmp <- unlist(expressionData[i,])
    exampleDataDevTmp <- unlist(expressionDataDev[i,])

    # Loop over three different algorithms.
    inferedRatesTmp <- sapply(c("L-BFGS-B","BFGS","Nelder-Mead"),function(m)
    {
      # Loop over all the initial rates configurations.
      foe <- sapply(initialRates,function(initialRatesTmp)
      {
        tryCatch(optim(par=initialRatesTmp # Initial rates.
                      ,FlagDev=FlagDev # Cost function selection.
                      ,fn=costFunction # Error function to be minimized.
                      ,dataGeneration=dataGeneration # Data generation function.
                      ,parFixed=parFixed # Fixed parameters.
                      ,data=exampleDataTmp # Experimental data.
                      ,dev=exampleDataDevTmp # Experimental data standard deviation.
                      ,lowB=lowB # Lower bound for rates.
                      ,upB=upB # Upper bound for rates.
                      ,TauFractions=TauFractions # Time points with cellular fractionation.
                      ,TauPoly=TauPoly # Time points with polysomal profiling.
                      ,TauTotal=TauTotal # Time points with total RNA profiling.
                      ,logOptim=logOptim # TRUE to optimize the model parameters in the Log space.
                      ,method=m # Optimization algorithm.
                      ,excludeSpecies=excludeSpecies # List of species to be excluded from the cost function.
                      ,lambda=lambda # Regularization strength.
                      ,control=list(maxit=1e9,reltol=1e-50)) # Optimization algorithm FIXED parameter.
          ,error=function(e){list("par"=rep(NaN,length(initialRatesTmp))
                                                                                   ,"value"=Inf
                                                                                   ,"counts"=NaN
                                                                                   ,"convergence"=NaN
                                                                                   ,"message"="Error!")})
      })
      # Selection of the best model among the different initial conditions.
      foe[,which.min(foe["value",])]
    })

    # Selection of the best model among the different algorithms.
    winnerTmp <- names(which.min(inferedRatesTmp["value",]))
    inferedRatesTmp <- inferedRatesTmp[,winnerTmp]

    return(inferedRatesTmp$par)
  },mc.cores=cpus))

  # From logarithmic to linear rates.
  if(logOptim){inferedRates <- exp(inferedRates)}
  rownames(inferedRates) <- rownames(expressionData)

  # Generation of the expression levels for the optimal set of rates.
  inferedData <- t(sapply(1:nrow(inferedRates),function(i)dataGeneration(exampleRates=inferedRates[i,] # Rates.
                                                                        ,TauFractions=TauFractions # Time points with cellular fractionation.
                                                                        ,TauPoly=TauPoly # Time points with polysomal profiling.
                                                                        ,TauTotal=TauTotal))) # Time points with total nascent and pre-existing RNA profiling.
  rownames(inferedData) <- rownames(inferedRates)

  # Loop to estimate for each gene the model Log Likelihood as a specific cost function.
  LogL <- t(sapply(1:nrow(expressionData),function(x){costFunction(par=inferedRates[x,] # Set of rates.
                                                                  ,data=expressionData[x,] # Experimental data.
                                                                  ,dev=expressionDataDev[x,] # Experimental data standard deviation.
                                                                  ,TauFractions = TauFractions # Time points with cellular fractionation.
                                                                  ,TauPoly = TauPoly # Time points with polysomal profiling.
                                                                  ,TauTotal = TauTotal # Time points with total RNA profiling.
                                                                  ,dataGeneration=dataGeneration # Data generation function.
                                                                  ,parFixed=NULL # List of rates to be excluded from the optimization.
                                                                  ,logOptim=FALSE # Always FALSE because we previously converted (potential) logarithmic rates.
                                                                  ,excludeSpecies=NULL # List of species to be excluded from the cost function.
                                                                  ,lowB=1e-6 # Lower bound for rates.
                                                                  ,upB=1e4 # Upper bound for rates.
                                                                  ,FlagDev="LL" # We are interested in LogLikelihood.
                                                                  ,lambda=0.05) # Regularization strength.
                                                      }))

  # AIC and BIC goodness of fit metrics.
  AIC <- 2*ncol(inferedRates) - 2*LogL
  BIC <- ncol(inferedRates)*log(ncol(expressionData))- 2*LogL

  # All goodness of fit metrics.
  metrics <- rbind(LL=LogL,AIC=AIC,BIC=BIC)
  rownames(metrics) <- c("LL","AIC","BIC")
  colnames(metrics) <- rownames(expressionData)

  # If data are simulated real rates correlations estimation.
  if(simulatedDatasetFlag)
  {
    realRates <- simulatedDataset$exampleRates
    
    print("Rates correlations:")
    print(round(cor(inferedRates,realRates,method="s"),2))

    cc <- round(cor(inferedRates,realRates,method="s"),2)

    relativeErrors <- cbind(inferedRates,realRates,(1-inferedRates/realRates))
    
    colnames(relativeErrors)[1:ncol(realRates)] <- paste0(colnames(realRates),".inf")
    colnames(relativeErrors)[(ncol(realRates)+1):(2*ncol(realRates))] <- paste0(colnames(realRates),".real")
    colnames(relativeErrors)[((2*ncol(realRates))+1):(3*ncol(realRates))] <- paste0(colnames(realRates),".err")

    idxs <- list(grep("inf",colnames(relativeErrors))
               ,grep("real",colnames(relativeErrors))
               ,grep("err",colnames(relativeErrors)))
    ratesCorrelations=round(sapply(idxs[[1]],function(i)cor(relativeErrors[,i],relativeErrors[,(max(idxs[[1]])+i)],method="s",use="c")),2)
  }else{relativeErrors=ratesCorrelations=NULL}

  # Output.
  list(inferedRates=inferedRates, inferedData=inferedData, expressionData=expressionData, metrics=t(metrics), relativeErrors=relativeErrors, ratesCorrelations=ratesCorrelations)
}

fullModelSimulationND_prem <- function(t,y,parms)
{
  k1=parms[["k1"]] # Synthesis.
  k2=parms[["k2"]] # Co-transcriptional splicing.
  k3=parms[["k3"]] # Detachment of spliced chromatin associated RNA.
  k4=parms[["k4"]] # Detachment of unspliced chromatin associated RNA.
  k5=parms[["k5"]] # Post-transcriptional splicing.
  k6=parms[["k6"]] # Export.
  k7=parms[["k7"]] # Translation.
  k8=parms[["k8"]] # Cytoplasmic decay
  k11=parms[["k11"]] # Nuclear decay.
  k10=parms[["k10"]] # Polysomal decay

  # Implementation of the ODE model.
  # System of equations for NASCENT RNA
  dy1=k1 - (k2+k4)*y[1]                     #y1=Chp
  dy2=k2*y[1] - k3*y[2]                     #y2=Chm
  dy3=k4*y[1] - (k5+k11)*y[3]               #y3=Np
  dy4=k3*y[2] + k5*y[3] - k6*y[4]           #y4=Nm
  dy5=k6*y[4] - (k7+k8)*y[5]                #y5=C
  dy6=k7*y[5] - k10*y[6]                    #y6=P             

  nascentTmp <- c(dy1,dy2,dy3,dy4,dy5,dy6)

  # System of equations for PRE-EXISTING RNA
  dy7=0 - (k2+k4)*y[7]                       #y7=Chp
  dy8=k2*y[7] - k3*y[8]                      #y8=Chm
  dy9=k4*y[7] - (k5+k11)*y[9]                #y9=Np
  dy10=k3*y[8] + k5*y[9] - k6*y[10]          #y10=Nm
  dy11=k6*y[10] - (k7+k8)*y[11]              #y11=C
  dy12=k7*y[11] - k10*y[12]                   

  preExistingTmp <- c(dy7,dy8,dy9,dy10,dy11,dy12)

  # Output.
  list(c(nascentTmp,preExistingTmp))
}

dataGenerationFullND_prem <- function(exampleRates # Set of rate.
                                ,TauFractions # Time points with cellular fractionation.
                                ,TauPoly # Time points with polysomal profiling.
                                ,TauTotal) # Time points with total RNA profiling.
{
  k1=exampleRates[["k1"]] # Synthesis.
  k2=exampleRates[["k2"]] # Co-transcriptional splicing.
  k3=exampleRates[["k3"]] # Detachment of spliced chromatin associated RNA.
  k4=exampleRates[["k4"]] # Detachment of unspliced chromatin associated RNA.
  k5=exampleRates[["k5"]] # Post-transcriptional splicing.
  k6=exampleRates[["k6"]] # Export.
  k7=exampleRates[["k7"]] # Translation.
  k8=exampleRates[["k8"]] # Cytoplasmic decay
  k11=exampleRates[["k11"]] # Nuclear decay.
  k10=exampleRates[["k10"]] # Polysomal decay

  # Initial conditions.
  yini <- c(y1=0,y2=0,y3=0,y4=0,y5=0,y6=0
           ,y7= k1/(k2+k4)
           ,y8= k1*k2/(k3*(k2+k4))
           ,y9= k4*k1/((k2+k4)*(k5+k11))
           ,y10= (k1/(k6*(k2+k4)))*(k2+(k4*k5)/(k5+k11))
           ,y11= (k1/((k7+k8)*(k2+k4)))*(k2+(k4*k5)/(k5+k11))
           ,y12= (k7/k10)*(k1/((k7+k8)*(k2+k4)))*(k2+(k4*k5)/(k5+k11)))

  # All time-points to be simulated
  TauGlobals <- unique(c(TauFractions,TauPoly,TauTotal))

  # Numerical solution of the ODE model.
  exampleData <- ode(times=TauGlobals # Time points.
                    ,y=yini # Initial conditions.
                    ,func=fullModelSimulationND_prem # Model.
                    ,parms=exampleRates) # Rates.

  # Removal of time points (first column).
  exampleData <- c(t(exampleData[,-1]))
  
  # Set of colnames.
  names(exampleData) <- c(sapply(TauGlobals,function(i)paste0(c("chpn","chmn","npn","nmn","cyn","pyn"
                                                               ,"chpp","chmp","npp","nmp","cyp","pyp"),i)))

  # Steady state solution (initial conditions).
  exampleData <- c(tail(yini,(length(yini)/2)),exampleData)
  names(exampleData)[1:(length(yini)/2)] <- c("chpt","chmt","npt","nmt","cyt","pyt")

  # Cytoplasmic RNA is re-defined as the sum of polysomal RNA and cytoplasmic RNA from the model
  # because this is closer to the experimental data.
  # Here we could add a flag to avoid this step if this assumption is not true for a specific experiment.
  exampleData[sort(names(exampleData[grep("^cyn",names(exampleData))]))] <- exampleData[sort(names(exampleData[grep("^cyn",names(exampleData))]))]+exampleData[sort(names(exampleData[grep("^pyn",names(exampleData))]))] 
  exampleData[sort(names(exampleData[grep("^cyp",names(exampleData))]))]  <- exampleData[sort(names(exampleData[grep("^cyp",names(exampleData))]))]+exampleData[sort(names(exampleData[grep("^pyp",names(exampleData))]))] 
  exampleData["cyt"] <- exampleData["cyt"]+exampleData["pyt"]

  exampleData <- c(exampleData[grep("t$",names(exampleData))] # We keep the total fractions
                   ,unlist(lapply(c(TauFractions,TauTotal),function(i)exampleData[grepl(paste0(i,"$"),names(exampleData))
                                                                          &!grepl(paste0("^pyn",i,"$"),names(exampleData))
                                                                          &!grepl(paste0("^pyp",i,"$"),names(exampleData))])) # We keep all the data for the profiled fraction (no polysomal which can be profiled at different time points).
                   ,unlist(lapply(TauPoly,function(i)exampleData[grepl(paste0("^pyn",i,"$"),names(exampleData))
                                                                |grepl(paste0("^pyp",i,"$"),names(exampleData))])))

  # If some fractionation time-points have been profiled only for total RNA we compute the total nascent and total
  # pre-existing, and we remove the corresponding fractions.
  if(!is.null(TauTotal))
  {
    exampleData <- c(exampleData
                    ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("n",i,"$"),names(exampleData))
                                                                 &!grepl(paste0("^pyn",i,"$"),names(exampleData))])}) # For certain time-points we just have nascent total RNA.
                    ,sapply(TauTotal,function(i){sum(exampleData[grepl(paste0("p",i,"$"),names(exampleData))
                                                                 &!grepl(paste0("^pyp",i,"$"),names(exampleData))])})) # For certain time-points we just have pre-existing total RNA.
    
    names(exampleData)[names(exampleData)==""] <- c(paste0("n",TauTotal),paste0("p",TauTotal))
    
    # Removal of fractions for time-points profiled without cellular fractionation.
    exampleData <- exampleData[!(names(exampleData)%in%c(sapply(TauTotal,function(i)paste0(c("chpn","chmn","npn","nmn","cyn","chpp","chmp","npp","nmp","cyp"),i))))]
    
    # Removal of eventual duplicated conditions
    exampleData <- exampleData[!duplicated(names(exampleData))]
  }
  
  exampleData <- exampleData[!duplicated(names(exampleData))]
  exampleData <- exampleData[sort(names(exampleData))]
}

simulateData <- function(exampleRates # Mean rates to sample.
                        ,TauFractions # Time points with cellular fractionation.
                        ,TauPoly # Time points with polysomal profiling.
                        ,TauTotal # Time points with total nascent and pre-existing RNA profiling.
                        ,noise # TRUE to add noise to data.
                        ,CV # Variation Coefficient for noise.
                        ,Reps # Number of replicates.
                        ,nGenes # Number of genes to be simulated.
                        ,seed # Seed for reproducibility.
                        ,ZeroThresh # Minimum expression value.
                        ,MultFact) # Number of nGenes to be modeled.
{
  if(length(intersect(TauFractions,TauTotal))>0)
  {
    print("Define independent Total and Fractions RNA temporal profiles.")
    stop()
  }

  # All time-points to be simulated
  TauGlobals <- unique(c(TauFractions,TauTotal,TauPoly))

  # Selection of the configuration to be simulated according to the provided example rates.
  # This means definition of the modelSimulation function (ode implementation) and of dataGeneration function (ode solution).
  
  if("k9"%in%names(exampleRates))
  { #Full Model with also nuclear decay
    modelSimulation <- fullModelSimulationND
    dataGeneration <- dataGenerationFullND
    print("Full Model with mature nuclear decay")
  }else if("k11"%in%names(exampleRates))
  { #Full Model with nuclear decay on premature
    modelSimulation <- fullModelSimulationND_prem
    dataGeneration <- dataGenerationFullND_prem
    print("Full Model with premature nuclear decay")
  }else if(all(c("k1","k2","k3","k4","k5","k6","k7","k8","k10")%in%names(exampleRates)))
  { #Full Model
    modelSimulation <- fullModelSimulation
    dataGeneration <- dataGenerationFull
    print("Full Model")
  }else if(all(c("k1","k2","k3","k4","k5","k6","k8")%in%names(exampleRates)))
  { #Full Model without translation
    modelSimulation <- modelSimulationNoPoly
    dataGeneration <- dataGenerationNoPoly
    print("Full Model without translation")
  }else if(all(c("k1","k2","k3","k6","k7","k8","k10")%in%names(exampleRates)))
  { #No Nucleoplasmic Premature - yes tranlation
    modelSimulation <- modelSimulationNoNucP
    dataGeneration <- dataGenerationNoNucP
    print("No Nucleoplasmic Premature - yes tranlation")
  }else if(all(c("k1","k2","k3","k6","k8")%in%names(exampleRates)))
  { #No Nucleoplasmic Premature - no tranlation
    modelSimulation <- modelSimulationNoPolyNoNucP
    dataGeneration <- dataGenerationNoPolyNoNucP
    print("No Nucleoplasmic Premature - no tranlation")
  }else if(all(c("k1","k3","k6","k7","k8","k10")%in%names(exampleRates)))
  { #No premature - yes translation
    modelSimulation <- modelSimulationNoP
    dataGeneration <- dataGenerationNoP
    print("No premature - yes translation")
  }else if(all(c("k1","k3","k6","k8")%in%names(exampleRates)))
  { #No premature - no translation
    modelSimulation <- modelSimulationNoPolyNoP
    dataGeneration <- dataGenerationNoPolyNoP
    print("No premature - no translation")
  }else{
    print("Error in configuration!")
  }

  # Definition of a specific name for the configuration for output purposes.
  elements <- c("Noise","Reps","CV","TauFractions","TauPoly","TauTotal",names(exampleRates))
  values <- c(noise,Reps,CV,paste0(TauFractions,collapse="_"),paste0(TauPoly,collapse = "_"),paste0(TauTotal,collapse = "_"),unname(exampleRates))
  nameTmp <- character(0)
  for (i in 1:length(elements)){
    name <- paste0(c(elements[i],values[i]),collapse = "_")
    nameTmp <- paste0(c(nameTmp,name),collapse ="_")
  }

  # Sampling of the rates.
  set.seed(seed)
  multipleExampleRates <- sapply(exampleRates,function(x)rnorm(nGenes*MultFact,mean=x,sd=5*x))

  #Selection of positive rate-sets.
  multipleExampleRates <- multipleExampleRates[apply(multipleExampleRates,1,function(row)all(row>=0)),]
  
  # Inclusion of the original rates set.
  multipleExampleRates <- rbind(exampleRates,multipleExampleRates)
  rownames(multipleExampleRates) <- 1:nrow(multipleExampleRates)

  # Data generation and addition of noise at the expression level.
  set.seed(seed)
  Obj <- t(apply(multipleExampleRates,1,function(par)
  {
    # Data generation.
    exampleDataTmp <- dataGeneration(exampleRates=par # Rates.
                                    ,TauFractions=TauFractions # Time points with cellular fractionation.
                                    ,TauPoly=TauPoly # Time points with polysomal profiling.
                                    ,TauTotal=TauTotal) # Time points with total nascent and pre-existing RNA profiling.

    if(noise)
    {
      # Noise addition and averages computation.
      ListDataTmp <- dataAverageFunction(v=exampleDataTmp # Expression levels.
                                        ,nRep=Reps # Number of replicates.
                                        ,CV=CV) # Variation Coefficient.
      return(ListDataTmp)
    }else{
      return(exampleDataTmp)
    }
  }))
  
  # Extraction of the quantities of interest.
  exampleData <- t(sapply(Obj,"[[",1))          #means 
  rownames(exampleData) <- 1:nrow(exampleData)
  DevDataTmp <- t(sapply(Obj,"[[",2))           #deviations
  rownames(DevDataTmp) <- 1:nrow(DevDataTmp)
  
  # Conditions with low coverage are set to the minimum coverage value.
  exampleData[abs(exampleData)<ZeroThresh] <- ZeroThresh

  # Selection of suitable examples without negative expression levels (not biologically reasonable).
  exampleData[which(exampleData<0)] <- NaN
  exampleData <- exampleData[apply(exampleData,1,function(r)all(is.finite(r))),]

  # Selection of at last nGenes examples.
  if(nrow(exampleData)<nGenes){print(paste0("Warning: less than ",nGenes," genes provided."))}
  exampleData <- exampleData[1:min(nrow(exampleData),nGenes),]
  multipleExampleRates <- multipleExampleRates[rownames(exampleData),]
  DevDataTmp <- DevDataTmp[rownames(exampleData),]

  rownames(exampleData) <- 1:nrow(exampleData)
  rownames(DevDataTmp) <- 1:nrow(DevDataTmp)
  rownames(multipleExampleRates) <- 1:nrow(multipleExampleRates)

  # Output.
  return(list("exampleData"=exampleData,"DevDataTmp" = DevDataTmp
              ,"exampleRates"=multipleExampleRates,"name"=nameTmp
              ,"TauGlobals"=TauGlobals,"TauFractions"=TauFractions
              ,"TauPoly"=TauPoly,"TauTotal"=TauTotal
              ,"modelSimulation"=modelSimulation,"dataGeneration"=dataGeneration))
}