## Barley Nursery Analysis
## 
## Summarize phenotypic data in the dataset
## 

# Base script
proj_dir <- getwd()
source(file.path(proj_dir, "startup.R"))

# Additional packages
library(GGally)
library(paletteer)
library(ggrepel)
library(cowplot)


# A list of traits to analyze
traits <- unique(pheno_dat$trait)



# Compare simple model with or without pedigree ---------------------------


# Iterate over nurseries, management, and traits





## Convert the pedigree mat to A
A <- as.matrix(pedigree_Amat[genotypes, genotypes]); rm(pedigree_Amat)
