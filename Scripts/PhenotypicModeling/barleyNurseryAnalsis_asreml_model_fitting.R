## ASREML modeling

library(asreml)
library(tibble)
library(dplyr)
library(stringr)
library(tidyr)

# Set the working directory
setwd("G:/AGRO/BARLEY_LAB/Jeff/Projects/BarleyNurseryAnalysis/")
proj_dir <- getwd()

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")

# Read or load data
load(file.path(data_dir, "data_for_asreml_modeling.RData"))

## Function to set factor to sum-to-zero contrasts
fct_contr_sum <- function(x, drop.levels = FALSE) {
  stopifnot(is.factor(x))
  stopifnot(is.logical(drop.levels))
  
  # Drop levels, if called for 
  x1 <- if (drop.levels) droplevels(x = x) else x
  
  # Redefine contrasts as sum-to-zero
  x1_contrasts <- contr.sum(levels(x1))
  colnames(x1_contrasts) <- head(levels(x1), -1)
  contrasts(x1) <- x1_contrasts
  return(x1)
}


## Sort by nursery and trait
asreml_out <- pheno_to_model %>%
  arrange(nursery, management, trait)



# Model data using FA covariance structure --------------------------------

# Then implement the algorithm for identifying clusters of non-COI
# environment as described by Burgeno et al 2008

# Options for asreml
asreml.options(fail = "soft", rotate.fa = FALSE, ai.sing = TRUE, maxit = 150, workspace = "1gb")

first_null <- min(which(sapply(asreml_out$out, is.null)))

## Testing
i = 1
##

## Iterate over rows in pheno_tomodel df
for (i in seq(first_null, nrow(asreml_out))) {
  
  df1 <- droplevels(asreml_out$data[[i]]) %>%
    # subset(location %in% sample_locs) %>%
    droplevels() %>%
    arrange(environment, line_name) %>%
    mutate_at(vars(environment, location, year), fct_contr_sum) %>%
    mutate(wts = 1 / replications)
  
  # Calculate the total phenotypic variance
  varP <- var(df1$value)

  ## Subset the A matrix ##
  A1 <- A[levels(df1$line_name), levels(df1$line_name), drop = FALSE]
  # Take the inverse
  A1_inv <- as.matrix(solve(A1)); dimnames(A1_inv) <- dimnames(A1)
  A1_inv <- structure(A1_inv, INVERSE = TRUE)

  
  # Fit simple models with fixed environment and random genotypes (with or without A mat)
  modelI <- asreml(fixed = value ~ 1 + environment,
                   random = ~line_name,
                   residual = ~ dsum(~ id(units) | environment),
                   weights = replications,
                   data = df1)
  
  modelA <- asreml(fixed = value ~ 1 + environment,
                   random = ~vm(line_name, A1_inv),
                   residual = ~ dsum(~ id(units) | environment),
                   weights = replications,
                   data = df1)
  
  
  ## Fit FA(2) models using environmental heterogeneous variances
  asreml_fit <- asreml(fixed = value ~ 1 + environment,
                       random = ~line_name:fa(environment, 2),
                       residual = ~ dsum(~ id(units) | environment),
                       weights = replications,
                       data = df1)

  ## Get the summary - contained parameter estimates
  fit_summary <- summary(asreml_fit)
  
  ## Get the BLUPs for genotypes
  fit_blups <- asreml_fit$coefficients$random %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    tbl_df() %>%
    separate(term, c("line_name", "term"), sep = ":") %>%
    mutate(line_name = str_remove(line_name, "line_name_")) %>%
    rename(estimate = bu)
  
  ## Separate into environment-specific blups and the genotype scores for the FA model
  fit_blups_pgv <- fit_blups %>%
    filter(str_detect(term, "Comp", negate = TRUE)) %>%
    mutate(term1 = str_extract(term, "location|year|environment"),
           level = lapply(term, FUN = str_split_fixed, pattern = "_", n = 2),
           level = sapply(level, "[[", 2)) %>%
    select(line_name, term = term1, level, estimate)
  
  fit_blups_scores <- fit_blups %>%
    filter(str_detect(term, "Comp")) %>%
    mutate(term1 = str_extract(term, "location|year|environment"),
           component = str_extract(term, "Comp[0-9]*") %>% str_replace(., "Comp", "fa")) %>%
    select(line_name, term = term1, component, estimate)
  
  
  ## Return the summary
  asreml_out$out[[i]] <- tibble(varP = varP, modelI_aic = summary(modelI)$aic, modelA_aic = summary(modelA)$aic,
                                summary = list(fit_summary), pgv_blups = list(fit_blups_pgv),
                                geno_scores = list(fit_blups_scores))
  
  
}
  
  
## Edit the df
analysis_models_fa <- asreml_out %>%
  select(-data) %>%
  unnest(out)

## Save
save("analysis_models_fa", file = file.path(result_dir, "phenotypic_analysis_asreml_fa.RData"))

