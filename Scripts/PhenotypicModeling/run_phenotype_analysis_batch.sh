#!/bin/bash

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Projects/BarleyNurseryAnalysis/Scripts/PhenotypicModeling


# Array of traits
TRAITS=(BetaGlucan DiastaticPower FreeAminoNitrogen GrainProtein GrainYield \
HeadingDate KernelWeight MaltExtract PlantHeight PlumpGrain SolubleProteinTotalProtein TestWeight)

# Set the queue settings
QUEUE_SETTINGS='-l walltime=36:00:00,mem=250gb,nodes=1:ppn=1'

# Set the node for job submission
QUEUE='amdsmall'


################################
## Do not edit below
################################


set -e
set -u
set -o pipefail

# Iterate over traits and launch the job

for tr in ${TRAITS[@]}; do

  echo "module load R/3.5.2_mkl && Rscript met_mixed_models.R $tr" \
  | qsub ${QUEUE_SETTINGS} -M neyha001@umn.edu -m abe -N met_mixed_model_$tr -r n -q $QUEUE

done