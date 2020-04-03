#!/bin/bash

#PBS -l walltime=36:00:00,mem=250gb,nodes=1:ppn=1
#PBS -N met_mixed_model_GrainProtein
# #PBS -N met_mixed_model_MaltExtract
# #PBS -N met_mixed_model_GrainYield
# #PBS -N met_mixed_model_HeadingDate
# #PBS -N met_mixed_model_PlantHeight
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n
#PBS -q amdsmall

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Projects/BarleyNurseryAnalysis/Scripts/PhenotypicModeling

module load R/3.5.2_mkl


## Run the R script
# GrainProtein
Rscript met_mixed_models.R GrainProtein

# MaltExtract
# Rscript met_mixed_models.R MaltExtract

# GrainYield
# Rscript met_mixed_models.R GrainYield

# HeadingDate
# Rscript met_mixed_models.R HeadingDate

# PlantHeight
# Rscript met_mixed_models.R PlantHeight


