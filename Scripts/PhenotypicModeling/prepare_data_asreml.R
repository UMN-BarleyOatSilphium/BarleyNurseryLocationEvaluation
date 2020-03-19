## Barley Nursery Analysis
## 
## Prepare data for use in ASREML
## - Calculate matrix inverses
## - Collect phenotype data to model
## 
## 

# Run on a local machine
proj_dir <- getwd()
source(file.path(proj_dir, "startup.R"))


## Subset more relevant traits
traits_tokeep <- c("GrainProtein", "MaltExtract", "GrainYield", "HeadingDate", "PlantHeight")

## Nest phenotypic data by trait, nursery, and management
pheno_dat_subset <- pheno_dat %>%
  filter(trait %in% traits_tokeep) %>%
  select(-trial) %>%
  left_join(., select(trial_metadata, trial, environment, nursery, management)) %>%
  mutate_at(vars(trial, line_name, environment), as.factor) %>%
  group_by(trait, nursery, management) %>%
  # Filter out low-observation combinations
  filter(n() >= 100) %>%
  ungroup() %>%
  droplevels()
  
  
  
## Nest
pheno_to_model <- pheno_dat_subset %>%
  group_by(trait, nursery, management) %>%
  nest() %>%
  mutate(out = list(NULL))


# Vector of genotype names
genotypes <- levels(pheno_dat_subset$line_name)

## Convert the pedigree mat to A; remove the larger matrix
A <- as.matrix(pedigree_Amat[genotypes, genotypes]); rm(pedigree_Amat)
# Convert to sparse
A <- Matrix(A, sparse = TRUE)
# Convert to matrix of row, column, value in row-major order
A_df <- summary(A) %>%
  arrange(i, j)



## Do the same thing for the inverse
Ainv <- solve(A)
Ainv_df <- summary(Ainv) %>%
  arrange(i, j) %>%
  as.matrix() %>%
  # Add rowNames attribute
  structure(., "INVERSE" = TRUE, rowNames = row.names(A))



## Save data for use by ASREML
save("pheno_to_model", "Ainv_df", file = file.path(data_dir, "data_for_asreml_modeling.RData"))
save("pheno_to_model", "Ainv_df", 
     file = "G:/AGRO/BARLEY_LAB/Jeff/Projects/BarleyNurseryAnalysis/Data/data_for_asreml_modeling.RData")




## Replicate over environments
Ainv_df1 <- bind_rows(replicate(450, Ainv_df, simplify = FALSE))









