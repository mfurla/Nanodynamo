#!/bin/bash

# Relevant features to replicate the simulated configurations:
#
#     reps=1,2,3,4,5 List of configurations with different number of replicates to be modeled (e.g. from 1 to 5).
#     TauFractions=0,0.33 List of labeling times to be modeled for cellular fractions (e.g. 0' and 20').
#     TauPoly=0,0.33 List of labeling times to be modeled for polysomal RNA (e.g. 0' and 20').
#     nGenes=1000 Number of genes to be modeled.
#     rates=k1,k2,k3,k4,k5,k6,k7,k8,k10 List of rates to be included in the model, this selects the model to be simulated (e.g. all the rates of the full model).
#
# Important technical parameters:
#     cpus=15 Number of threads to be used for parallelization (e.g. 15 CPUs).

# Example 1: Simulation of the full model, 1000 genes, 1 to 5 replicates, one labeling time of 20'
Rscript /path/to/dataSimulationMain.R reps=1,2,3,4,5 TauFractions=0,0.33 TauPoly=0,0.33 translationCoeff=0.5 Flags=FC nGenes=1000 cpus=15 ZeroThresh=0.0000000001 MultFact=100000 rates=k1,k2,k3,k4,k5,k6,k7,k8,k10

# Example 2: Simulation of a model missing polysomal RNA, 1000 genes, 1 replicate, three labeling times of 20' 60' and 120'
Rscript /path/to/dataSimulationMain.R reps=1,2,3,4,5 TauFractions=0,0.33,1,2 TauPoly=0,0.33,1,2 translationCoeff=0.5 Flags=FC nGenes=1000 cpus=15 ZeroThresh=0.0000000001 MultFact=100000 rates=k1,k2,k3,k4,k5,k6,k8

conda deactivate