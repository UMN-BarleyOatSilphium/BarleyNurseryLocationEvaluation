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


# Trait-trial metadata to use
trial_trait_metadata1 <- trial_trait_metadata %>%
  mutate(replications = ifelse(trait %in% quality_traits, 1, replications)) %>%
  distinct_at(vars(trial, trait, replications))


## Nest phenotypic data by trait, nursery, and management
pheno_dat_subset <- pheno_dat %>%
  # Add nursery and management variables
  mutate(nursery = tolower(str_sub(trial, 1, 3)),
         management = tolower(str_extract(trial, "Rainfed|Irrigated"))) %>%
  # Add replication information
  left_join(., trial_trait_metadata1) %>%
  mutate_at(vars(trial, line_name, environment), as.factor)


## Plot histograms of trait values
pheno_dat_subset %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~ trait, scales = "free") +
  theme_acs()


pheno_dat_subset1 <- pheno_dat_subset %>%
  # Keep only the important traits
  filter(trait %in% traits_keep) %>%
  # Logit transform plump/thin grains
  mutate(value = case_when(
    trait %in% c("PlumpGrain", "ThinGrains") ~ logitTransform(value / 100),
    TRUE ~ value
  )) %>%
  # Center and scale data within trials
  group_by(trial, trait) %>%
  mutate(value_scale = as.numeric(scale(value))) %>%
  ungroup()

  
## Nest
pheno_to_model <- pheno_dat_subset1 %>%
  group_by(trait, nursery, management) %>%
  nest() %>%
  mutate(out = list(NULL))


# Vector of genotype names
genotypes <- levels(pheno_dat_subset1$line_name)

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
# save("pheno_to_model", "Ainv_df", "A", file = file.path(data_dir, "data_for_asreml_modeling.RData"))
save("pheno_to_model", "Ainv_df", "A",
     file = "Z:/BARLEY_LAB/Jeff/Projects/BarleyNurseryAnalysis/Data/data_for_asreml_modeling.RData")



