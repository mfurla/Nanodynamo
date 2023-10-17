# Nanodynamo
**Nanodynamo: studying the dynamics of RNA metabolism through mathematical modelling and native RNA profiling.**

This repository contains files and scripts required to reproduce the analyses presented in the manuscript.

The _condaEnvironment.yml_ file allows to reproduce the Anaconda environment used for our analyses.

In the _Scripts_ folder we have collected:

- A Julia script used to evaluate parameters identifiability (parametersIdentifiability.ipynb),
- An R script to plot polysome profiling data (polysomalProfiling.R),
- An R script to investigate nascent RNA saturation with simulated data (nascentSaturation.R),
- A bash script to simulate genes with Nanodynamo and analyze the resulting dataset (dataSimulationBash.sh),
- An R script to simulate genes with Nanodynamo and analyze the resulting dataset (dataSimulationMain.R),
- An R script to reproduce our analyses on simulated data (dataSimulationAnalysis.R),
- An R script to characterize nano-ID accuracy (nanoIDAccuracy.R),
- An R script to characterize nascent and premature reads fractions (fractionsReadsClassification.R),
- An R script to estimate rates according to the nano-ID model (nanoIDMain.R),
- An R script to estimate rates with INSPEcT (INSPEcTAnalysisAllGenes.R),
- Three R scripts with core Nanodynamo's funcitons (allInternalFunctions.R, allInternalFunctionsPatch_nuclearPrematureDecay.R, allInternalFunctionsPatch_nucleoplasmicPrematureExport.R),
- An R script to model genes with simplified models (sameGenesSimplifiedModels.R),
- An R script to reproduce our analyses on rates couplings (couplings.R),
- An R script to model the experimental datasets (experimentalDataAnalysis.R),
- An R script to reproduce our analyses on experimental data (mainModelingResults.R).

~ _mainModelingResults.R_: Figures 2A-D, 2E-F, 3D-G, 4A-B, 4D-F, 5A, 5C-F, S4, S7-9, S11-14, S17-27A, S32.  
~ _nanoIDAccuracy.R_: Figure S1.  
~ _nascentSaturation.R_: Figure S3.  
~ _variationCoefficient.R_: Figure S2.  
~ _polysomalProfiling.R_: Figure S6.  
~ _fractionsReadsClassification.R_: Figures 3B, S10, S16.  
~ _couplings.R_: Figures 6A-E, S27-29.  
~ _dataSimulationAnalysis.R_: Figures 1C-F, S4, S12A, S30-31, Supplemental Methods Figures.  
~ _sameGenesSimplifiedModels.R_: Supplemental Methods Figures.  

In the _Results_ folder we have collected:

- The rates inferred with the nano-ID model for the Untreated condition (nanoIDRates.rds),
- Rates and gene expression levels estimated with Nanodynamo for the Untreated condition (Untreated),
- Rates and gene expression levels estimated with Nanodynamo for the Pladienolide B treated condition (PladienolideB),
- Rates and gene expression levels estimated with Nanodynamo for the Leptomycin B treated condition (LeptomycinB),
- Rates and gene expression levels estimated with Nanodynamo for the Harringtonin condition (Harringtonin),

In the _Files_ folder we have collected:

- Polysome profiling results (polysomalProfiling),
- PolyA tail length quantifications (PolyATails),
- Data produced during nano-ID training (nanoID.modification.probabilities.labeled.RData,nanoID.modification.probabilities.unlabeled.RData
,nanoIDResults.rds),
- UTRs sequences (utr3.fa,utr5.fa),
- Data required to estimate gene expression levels variation coefficients (3T9MycER_rpkms.RData),
- Data required to estimate rates initial conditions for simulations and models fit (regulated_genes_features.xls).
