## Barley Nursery Analysis
## 
## Fit models to analyze the nursery data
## 
## Author: Jeff Neyhart
## Last modified: 25 Feb. 2020
## 

## This script will fit mixed-models to analyze data across the barley nursery set. One model will be
## fitted per trait/nursery/management combination.
## 
## Mixed-model breakdown:
## - Environments are fixed
## - Genotypes and GxE are random with constant variance
## - Heterogenous error per environment
## 


# # Run on a local machine
# proj_dir <- getwd()
# source(file.path(proj_dir, "startup.R"))

# Run the source script
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Projects/BarleyNurseryAnalysis"
source(file.path(proj_dir, "startup_MSI.R"))



# Other packages
library(lmerTest)
library(broom)
library(Matrix)
library(modelr)
library(nlme)

## Convert the pedigree mat to A
A <- as.matrix(pedigree_Amat[levels(pheno_dat$line_name), levels(pheno_dat$line_name)])


## Nest phenotypic data by trait, nursery, and management
pheno_to_model <- pheno_dat %>%
  left_join(., select(trial_metadata, trial, nursery, management)) %>%
  group_by(trait, nursery, management) %>%
  nest() %>%
  mutate(out = list(NULL))

## Iterate over rows in this data.frame
for (i in seq_len(nrow(pheno_to_model))) {
 
  df1 <- droplevels(pheno_to_model$data[[i]])
  
  # ## Subset data to demo
  # df1 <- subset(df1, year1 < 2008) %>%
  #   droplevels()
  
  # Subset the A matrix
  A_use <- A %>% 
    subset(row.names(.) %in% levels(df1$line_name), colnames(.) %in% levels(df1$line_name))
  
  ## Build a relationship matrix for GxE
  # E matrix for environments
  E <- diag(nlevels(df1$environment)); dimnames(E) <- replicate(2, levels(df1$environment), simplify = F)
  AE <- kronecker(A_use, E, make.dimnames = TRUE)
  
  # Fit the model
  fit <- mmer(fixed = value ~ 1 + environment,
              random = ~vs(line_name, Gu = A_use) + vs(line_name:environment, Gu = AE),
              rcov = ~vs(ds(environment), units),
              data = df1, date.warning = FALSE)
  
  # Heritability
  H2 <- herit.mmer(x = fit)
  
  # Extract blups and variance components
  BLUPs <- fit$U %>% 
    map("value") %>% 
    map(~tibble(line_name = names(.x), blup = .x)) %>% 
    modify_at(2, .f = ~separate(., line_name, c("line_name", "environment"), sep = ":")) %>%
    reduce(full_join, by = "line_name") %>% 
    mutate(pred_value = blup.x + blup.y) %>% 
    select(line_name, environment, pred_value)
  
  # Extract variance components
  var_comp <- fit$sigma %>% 
    tibble(component = names(.), variance = .) %>% 
    unnest(variance)
  
  
  # Assemble into a df
  pheno_to_model$out[[i]] <- tibble(H2 = H2, var_comp = list(var_comp), BLUPs = list(BLUPs))
  
}


## Unnest matrix
met_mm_out <- pheno_to_model %>%
  select(-data) %>% 
  unnest(out)

# Save
save("met_mm_out", file = file.path(result_dir, "met_mixed_model_output.RData"))






## Example output for multiple environments and multiple traits
example_traits <- c("GrainYield", "GrainProtein")

# ## Extract phenotypic data for these traits
# ## Find trials where at least one trait was observed; then find
# ## trials where both traits were observed
# trait_pairs <- bind_rows(
#   tibble(trait1 = example_traits, trait2 = trait1),
#   rename_all(as.data.frame(t(combn(example_traits, 2)), stringsAsFactors = FALSE), ~c("trait1", "trait2"))
# )
# 
# # Phenotype data to use
# pheno_dat_use <- left_join(pheno_dat, select(trial_metadata, trial, nursery, management), by = "trial")
# 
# pheno_to_model_example <- trait_pairs %>%
#   mutate(data = map2(trait1, trait2, ~{
#     pheno_dat_use %>%
#       filter(trait %in% c(.x, .y)) %>%
#       # Filter trials where both traits were observed
#       group_by(trial) %>%
#       filter(n_distinct(trait) == n_distinct(c(.x, .y))) %>%
#       spread(trait, value) %>%
#       group_by(nursery, management) %>%
#       nest()
#   })) %>% 
#   unnest(data)



### Example analysis to determine the best model ####

# For each of 5 trait, sample 5 locations with at least 5 years of data
nLocations <- 10
nYears <- 5


## Find traits with the required number of locations and years
pheno_dat_example <- pheno_dat %>% 
  filter(year %in% 2000:2010) %>%
  group_by(location) %>% 
  filter(n_distinct(year) >= nYears) %>% 
  group_by(trait) %>% 
  filter(n_distinct(location) >= nLocations) %>%
  ungroup()
  
set.seed(1254)
sample_traits <- c(example_traits, sample(x = setdiff(unique(pheno_dat_example$trait), example_traits), 
                                          size = 5 - length(example_traits)))




## Nest phenotypic data by trait
set.seed(1254)
pheno_to_model_example <- pheno_dat_example %>%
  filter(trait %in% sample_traits) %>% 
  split(.$trait) %>% 
  map_df(~filter(., location %in% sample(unique(location), nLocations))) %>%
  left_join(., select(trial_metadata, trial, nursery, management)) %>%
  group_by(trait) %>%
  nest() %>%
  mutate(out = list(NULL))


## Iterate over rows in this data.frame
for (i in seq_len(nrow(pheno_to_model_example))) {

  df1 <- droplevels(pheno_to_model_example$data[[i]])
  # Trait or traits
  tr <- str_subset(names(df1), "^[A-Z]")

  # ## Subset data to demo
  # df1 <- subset(df1, year1 < 2008) %>%
  #   droplevels()

  # Subset the A matrix
  A_use <- A %>%
    subset(row.names(.) %in% levels(df1$line_name), colnames(.) %in% levels(df1$line_name))

  ## Build a relationship matrix for GxE
  # E matrix for environments
  # E <- diag(nlevels(df1$environment)); dimnames(E) <- replicate(2, levels(df1$environment), simplify = F)
  E <- Diagonal(nlevels(df1$environment)); dimnames(E) <- replicate(2, levels(df1$environment), simplify = F)
  
  # E <- select(df1, line_name, environment, value) %>%
  #   spread(environment, value) %>%
  #   select(levels(df1$environment)) %>%
  #   cor(., use = "pairwise.complete.obs")
  AE <- as.matrix(kronecker(A_use, E, make.dimnames = TRUE))
  dimnames(AE) <- replicate(2, paste(rep(row.names(A_use), each = ncol(E)), row.names(E), sep = ":"), simplify = F)

  
  # Fit the model
  fit <- mmer(fixed = value ~ 1 + environment,
              random = ~vs(line_name, Gu = A_use) + vs(line_name:environment, Gu = AE),
              rcov = ~vs(ds(environment), units),
              data = df1, date.warning = FALSE)

  # Heritability
  H2 <- herit.mmer(x = fit, method = "Cullis")

  ## Get BLUEs (environment main effects/means)
  BLUEs <- coef(fit) %>%
    select(environment = Effect, estimate = Estimate) %>%
    mutate(environment = levels(df1$environment),
           estimate = c(estimate[1], estimate[1] + estimate[-1]))

  # Extract blups and variance components
  BLUPs <- fit$U %>%
    map("value") %>%
    map(~tibble(line_name = names(.x), blup = .x)) %>%
    modify_at(2, .f = ~separate(., line_name, c("line_name", "environment"), sep = ":")) %>%
    reduce(full_join, by = "line_name") %>%
    mutate(pgv_g_nonzero = blup.x != 0,
           pgv_ge_nonzero = blup.y != 0,
           pgv = blup.x + blup.y) %>%
    select(line_name, environment, pgv_g_nonzero, pgv_ge_nonzero, pgv)

  # Generate predicted phenotypic values
  PPV <- full_join(BLUPs, BLUEs, by = "environment") %>%
    mutate(ppv = pgv + estimate) %>%
    left_join(., select(df1, line_name, environment, pheno = value)) %>%
    select(-estimate)


  # Extract variance components
  var_comp <- fit$sigma %>%
    tibble(component = names(.), variance = .) %>%
    unnest(variance)

  # Extract model summary
  model_summ <- tibble(AIC = fit$AIC, BIC = fit$BIC, LL = last(fit$monitor[1,]))


  # ## Extract useful components from the model fit:
  # ## PEV
  # ## Model suitability estimates
  # ##
  # ##
  # fit1 <- fit[c("VarU", "PevU")]


  # Assemble into a df
  pheno_to_model_example$out[[i]] <- tibble(H2 = H2, var_comp = list(var_comp), pred_pheno = list(PPV),
                                            model_summ = list(model_summ))

}

## Unnest df
met_mm_out_example <- pheno_to_model_example %>%
  select(-data) %>%
  unnest(out) %>%
  select(-fit_comp)

# Save
save("met_mm_out_example", file = file.path(result_dir, "met_mixed_model_output_example.RData"))



