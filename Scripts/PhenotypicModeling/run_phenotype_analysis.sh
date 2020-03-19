#!/bin/bash

#PBS -l walltime=24:00:00,mem=150gb,nodes=1:ppn=1
#PBS -N met_mixed_model
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n
#PBS -q amdsmall

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Projects/BarleyNurseryAnalysis/Scripts/PhenotypicModeling

module load R/3.5.2_mkl

# Met analysis
Rscript met_mixed_models.R

