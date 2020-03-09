#!/bin/bash

#PBS -l walltime=48:00:00,mem=62gb,nodes=1:ppn=1
#PBS -N met_mixed_model
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Projects/BarleyNurseryAnalysis

module load R/3.5.0
# module load R/3.5.2_mkl

# Met analysis
Rscript met_mixed_models.R

