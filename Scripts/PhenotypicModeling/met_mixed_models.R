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
library(broom)
library(Matrix)
library(modelr)

# Get the script arguments 
# arg1: traits
# arg2: nursery
# arg3: management
args <- commandArgs(trailingOnly = TRUE)
# Create an argument filter expression
arg1 <- if (is.na(args[1])) unique(pheno_dat$trait) else args[1]
arg2 <- if (is.na(args[2])) unique(trial_metadata$nursery) else args[2]
arg3 <- if (is.na(args[3])) unique(trial_metadata$management) else args[3]

## Subset more relevant traits
traits_tokeep <- c("GrainProtein", "MaltExtract", "GrainYield", "HeadingDate", "PlantHeight",
                   "BetaGlucan", "DiastaticPower", "FreeAminoNitrogen", "KernelWeight",
                   "PlumpGrain", "SolubleProteinTotalProtein", "TestWeight")

# Minimum number of years in which a location was observed
min_year <- 5
# Minimum number of locations in which a trait was observed
min_loc <- 5


## Nest phenotypic data by trait, nursery, and management
pheno_to_model <- pheno_dat %>%
  select(-trial) %>%
  left_join(., select(trial_metadata, trial, environment, nursery, management)) %>%
  filter(trait %in% traits_tokeep) %>%
  # Filter from args
  filter(trait %in% arg1, nursery %in% arg2, management %in% arg3) %>%
  mutate_at(vars(trial, line_name, environment), as.factor)


## Apply some filters to the dataset
## 1. locations must be observed in >= 5 years
pheno_to_model <- pheno_to_model %>% 
  group_by(trait, nursery, management, location) %>% 
  filter(n_distinct(year1) >= min_year) %>%
  group_by(trait, nursery, management) %>% 
  filter(n_distinct(location) >= min_loc) %>%
  ungroup() %>%
  droplevels()

# # # ## Assess observations
# pheno_to_model %>%
#   group_by(trait, nursery, management) %>%
#   summarize_at(vars(line_name, environment, location, year), n_distinct) %>%
#   mutate(n_pred_obs = line_name * environment) %>%
#   View


## Set the model expression to evaluatate
model_exp <- expression({
  fit <- mmer(fixed = value ~ 1 + environment,
              random = ~vs(line_name, Gu = A) + vs(line_name:environment, Gu = AE),
              rcov = ~vs(ds(environment), units),
              data = df1, 
              date.warning = FALSE)
})


## Nest
pheno_to_model <- pheno_to_model %>%
  group_by(trait, nursery, management) %>%
  nest() %>%
  mutate(out = list(NULL))


## Iterate over rows in this data.frame
for (i in seq_len(nrow(pheno_to_model))) {
 
  df1 <- droplevels(pheno_to_model$data[[i]])
  
  # ## Subset data to demo
  # df1 <- subset(df1, year1 < 2004) %>%
  #   droplevels()
  
  # Get a list of genotypes
  genotypes <- levels(df1$line_name)
  
  
  # Subset the pedigree relationship matrix
  A <- as.matrix(pedigree_Amat[genotypes, genotypes])
  
  ## Build a relationship matrix for GxE
  # E matrix for environments
  # E <- diag(nlevels(df1$environment)); dimnames(E) <- replicate(2, levels(df1$environment), simplify = F)
  E <- Diagonal(nlevels(df1$environment)); dimnames(E) <- replicate(2, levels(df1$environment), simplify = F)
  
  # E <- select(df1, line_name, environment, value) %>%
  #   spread(environment, value) %>%
  #   select(levels(df1$environment)) %>%
  #   cor(., use = "pairwise.complete.obs")
  AE <- as.matrix(kronecker(A, E, make.dimnames = TRUE))
  dimnames(AE) <- replicate(2, paste(rep(row.names(A), each = ncol(E)), row.names(E), sep = ":"), simplify = F)
  
  ### Fit the model ###
  
  # # Example
  # fit <- mmer(fixed = value ~ 1 + environment,
  #             random = ~vs(ds(environment), line_name, Gu = A),
  #             rcov = ~vs(ds(environment), units),
  #             data = df1, 
  #             date.warning = FALSE)
  
  
  
  ## Try to fit the model; capture the output
  model_try <- try( model_stdout <- capture.output({ eval(model_exp) }), silent = TRUE )
  # Print the output
  print(model_stdout)
  
  ## If model try is an error, skip
  if (class(model_try) == "try-error") {
    
    H2 <- var_comp <- PPV <- NULL
    model_summ <- model_try
    
    
  } else {
  
    # If model fit is empty, try using a smaller number of iterations; for instance find
    # the maximum logLik and use those iterations
    itry <- 1
    while (is_empty(fit) & itry == 1) {
      
      # Find the number of iterations that maximized the logLik
      best_iter <- model_stdout %>% 
        subset(., str_detect(., "singular", negate = T)) %>% 
        read_table(.) %>%
        subset(., LogLik == max(LogLik), iteration, drop = TRUE)
      
      # Refit
      eval(model_exp)
                  
      # Increase the counter
      itry = itry + 1
      
    }
  
  
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
    
  }
  
  
  # Assemble into a df
  pheno_to_model$out[[i]] <- tibble(H2 = list(H2), var_comp = list(var_comp), pred_pheno = list(PPV),
                                    model_summ = list(model_summ))
  
  ## Print a message
  cat(paste0("\n\nModel fitted for dataset: ", paste0(unlist(pheno_to_model[i,1:3]), collapse = ", ")))
  
}


## Unnest matrix
met_mm_out <- pheno_to_model %>%
  select(-data)

# Save
filename <- file.path(result_dir, paste0(c(arg1, "met_mixed_model_output.RData"), collapse = "_"))
save("met_mm_out", file = filename)


