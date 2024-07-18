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
- Two R scripts with core Nanodynamo's funcitons (allInternalFunctions.R, allInternalFunctionsPatch_nuclearPrematureDecay.R),
- An R script to model genes with simplified models (sameGenesSimplifiedModels.R),
- An R script to model the experimental datasets (experimentalDataAnalysis.R),
- An R script to reproduce our analyses on experimental data (mainModelingResults.R),
- An R script to estimate RNA species CVs from experimental data (CoefficientOfVariations.R).

~ _mainModelingResults.R_: Figures 2B-E, 3E-I, 4B-E, 4G, 5B-F, 6A-D, S7-9, S11, S12B, S13-14, S18-24, S26, S28-36, S39-41.
~ _nanoIDAccuracy.R_: Figure S2.  
~ _nascentSaturation.R_: Figure S3.  
~ _polysomalProfiling.R_: Figures 2A, 3D, 4A, 5A.  
~ _fractionsReadsClassification.R_: Figures 3C, S10, S17.  
~ _couplings.R_: Figures 6A-E, S27-29.  
~ _dataSimulationAnalysis.R_: Figures 1D, 1E, S4-5, S12A, S37-38, S43-52.
~_INSPEcTAnalysisAllGenes.R_: Figures S53-56.

In the _Results_ folder we have collected:

- Rates and gene expression levels estimated with Nanodynamo for the Untreated condition (Untreated),
- Rates and gene expression levels estimated with Nanodynamo for the Pladienolide B treated condition (PladienolideB),
- Rates and gene expression levels estimated with Nanodynamo for the Leptomycin B treated condition (LeptomycinB),
- Rates and gene expression levels estimated with Nanodynamo for the Harringtonine condition (Harringtonine),
- Simulated data and corresponding inferred rates.

In the _Files_ folder we have collected:

- The rates inferred with the nano-ID model for the Untreated condition (nanoIDResults.rds),
- Polysome profiling results (polysomalProfiling),
- PolyA tail length quantifications (PolyATails),
- Data produced during nano-ID training (nanoID.modification.probabilities.labeled.RData,nanoID.modification.probabilities.unlabeled.RData
,nanoIDResults.rds),
- UTRs sequences (utr3.fa,utr5.fa),
- Data required to estimate rates initial conditions for simulations and models fit (regulated_genes_features.xls),
- Encode samples IDs used for the identification of RBPs and TFs targets (encodeFilesRBPs.txt, encodeFilesTFs.txt
),
- RNA species CVs (CVs.rds).
