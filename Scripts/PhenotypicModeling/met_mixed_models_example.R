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


# Run on a local machine
proj_dir <- getwd()
source(file.path(proj_dir, "startup.R"))


# Other packages
library(lmerTest)
library(broom)
library(Matrix)
library(modelr)
library(nlme)


## Subset more relevant traits
traits_tokeep <- c("GrainProtein", "MaltExtract", "GrainYield", "HeadingDate", "PlantHeight")


## Nest phenotypic data by trait, nursery, and management
pheno_dat_subset <- pheno_dat %>%
  filter(trait %in% traits_tokeep) %>%
  select(-trial) %>%
  left_join(., select(trial_metadata, trial, environment, nursery, management)) %>%
  mutate_at(vars(trial, line_name, environment), as.factor)

## Minimum number of years
min_years <- 5

## Example output for multiple environments and multiple traits
example_traits <- c("GrainYield", "GrainProtein")


## Apply some filters to the dataset
## 1. locations must be observed in >= 5 years
pheno_dat_subset1 <- pheno_dat_subset %>% 
  filter(trait %in% example_traits) %>%
  group_by(trait, nursery, management, location) %>% 
  filter(n_distinct(year1) >= min_years) %>%
  ungroup() %>%
  droplevels()

pheno_tabs <- distinct(pheno_dat_subset1, location, year, trait, nursery, management) %>%
  unite(group, nursery, management) %>%
  xtabs(~ location + year + group, data = .)

## Sum years for the numbers of 2s
year_score <- apply(X = pheno_tabs, MARGIN = c(2,3), function(x) sum(x == n_distinct(pheno_dat_subset1$trait)))

# List of sequencies
year_seq_list <- range(pheno_dat_subset1$year1) %>% 
  {seq(from = first(.), to = last(.), by = 1)} %>% 
  accumulate(., ~tail(c(.x, .y), min_years)) %>%
  subset(., map_lgl(., ~length(.) == min_years)) %>%
  map(as.character)

## Find the sequence with the highest score
year_score_summ <- map(year_seq_list, ~apply(X = year_score[.,], MARGIN = 2, FUN = sum))
best_seq_index <- apply(X = do.call("rbind", year_score_summ), MARGIN = 2, FUN = which.max)
best_seq <- map(best_seq_index, ~year_seq_list[[.]])


## Subset those years for all traits
pheno_dat_example <- best_seq %>%
  imap_dfr(~filter(pheno_dat_subset1, year %in% .x, nursery == str_split(.y, pattern = "_")[[1]][1], 
               management == str_split(.y, pattern = "_")[[1]][2])) %>%
  droplevels()
       



# Get a list of genotypes
genotypes <- levels(pheno_dat_example$line_name)

## Convert the pedigree mat to A
A <- as.matrix(pedigree_Amat[genotypes, genotypes])


## Nest
pheno_to_model <- pheno_dat_example %>%
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
  # E <- diag(nlevels(df1$environment)); dimnames(E) <- replicate(2, levels(df1$environment), simplify = F)
  E <- Diagonal(nlevels(df1$environment)); dimnames(E) <- replicate(2, levels(df1$environment), simplify = F)
  
  # E <- select(df1, line_name, environment, value) %>%
  #   spread(environment, value) %>%
  #   select(levels(df1$environment)) %>%
  #   cor(., use = "pairwise.complete.obs")
  AE <- as.matrix(kronecker(A_use, E, make.dimnames = TRUE))
  dimnames(AE) <- replicate(2, paste(rep(row.names(A_use), each = ncol(E)), row.names(E), sep = ":"), simplify = F)
  
  ### Fit the model ###
  
  fit <- mmer(fixed = value ~ 1 + environment,
              random = ~vs(line_name, Gu = A_use) + vs(line_name:environment, Gu = AE),
              rcov = ~vs(ds(environment), units),
              data = df1, 
              date.warning = FALSE)

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
  pheno_to_model$out[[i]] <- tibble(H2 = H2, var_comp = list(var_comp), pred_pheno = list(PPV),
                                    model_summ = list(model_summ))
  
  
}


## Unnest df
met_mm_out_example <- pheno_to_model %>%
  select(-data) %>%
  unnest(out)

# Save
save("met_mm_out_example", file = file.path(result_dir, "met_mixed_model_output_example.RData"))



