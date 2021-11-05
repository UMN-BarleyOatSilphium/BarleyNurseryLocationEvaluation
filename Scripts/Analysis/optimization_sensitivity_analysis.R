## Barley Nursery Location Analysis
## 
## Location ranking and optimization
## 
## This script will use the summarized data from the analyze_locations.R script
## to perform location ranking and optimization.
## 
## 


# Base script
proj_dir <- getwd()
source(file.path(proj_dir, "startup.R"))

# Additional packages
library(gramEvol)

# A list of traits to analyze
traits <- unique(pheno_dat$trait)

# Minimum number of environments to include a location
min_env <- 3

# Load the location analysis data
load(file.path(result_dir, "location_analysis_results.RData"))



 
# Optimize the selection of locations for all traits per nursery
# 
# For each location penalty level, select a random set of locations and measure
# fitness components

# Vector of location penalties
# loc_penalties <- c(0, 0.01, 0.05, 0.1, 0.5)
loc_penalties <- 0.01 # Use only the penalty deemed optimal

# List of trait weights
trait_weights <- list(
  equal = setNames(rep(1, length(traits_keep)), traits_keep)
)



# Perform a sensitivity analysis of the metric weights --------------------


# First create a df defining the parameter range
sensitivity_param_df <- crossing(varY_weight = seq(0, 1, by = 0.25), repAvg_weight = varY_weight, 
                                 reprAvg_weight = varY_weight) %>%
  # Remove rows where all the weights are zero
  filter(apply(X = ., MARGIN = 1, FUN = function(row) any(row != 0))) %>%
  # # # Remove rows where all the weights are the same
  # # filter(apply(X = ., MARGIN = 1, FUN = function(row) (all(row != 0) & n_distinct(row) > 1) | (any(row == 0) & n_distinct(row) == 3))) %>%
  # filter(apply(X = ., MARGIN = 1, FUN = function(row) all(row == 0) | all(row != 0) | (any(row == 0) & n_distinct(row) == 3))) %>%
  # # Add equal-weight rows when 1 weight is 0
  # add_row(varY_weight = c(0, 1, 1), repAvg_weight = c(1, 0, 1), reprAvg_weight = c(1, 1, 0)) %>%
  # Add location penalties; add nursery and management
  crossing(., loc_penalty = loc_penalties, distinct(location_opimization_input, nursery, management)) %>%
  # Add output list
  mutate(optimization = list(NULL)) # %>% filter(loc_penalty == 0.01, nursery == "mvn")


# Iterate over these permutations
first_null <- min(which(sapply(X = sensitivity_param_df$optimization, FUN = is.null)))
index <- seq(first_null, nrow(sensitivity_param_df))
pb <- progress::progress_bar$new(total = length(index))
for (i in index) {
  
  # Get the parameters
  pen <- sensitivity_param_df$loc_penalty[i]
  comp_weight <- select(sensitivity_param_df, contains("weight"))[i,] %>% 
    rename_all(~str_remove(., "_weight")) %>% 
    unlist()
  
  # Subset the location optimization input
  df <- left_join(sensitivity_param_df[i,], location_opimization_input, by = c("nursery", "management"))
  
  # Extract lists of matrices
  GL <- df$GL
  RL <- df$RL
  corL <- df$repeatability
  repL <- df$representativeness2
  
  # Scale the matrices to 0-1
  RL1 <- mapply(GL, RL, FUN = function(.x, .y) scale2(x = list(.x, .y), range = "01"), SIMPLIFY = FALSE)
  GL1 <- lapply(RL1, "[[", 1)
  RL1 <- lapply(RL1, "[[", 2)
  # Scale everything to 0-1
  corL1 <- lapply(X = corL, FUN = scale2, range = "01")
  repL1 <- lapply(X = repL, FUN = scale2, range = "01")

  # List of traits
  tr_list <- df$trait
  # Trait weights are all equal
  tr_weights <- trait_weights$equal
  
  # Vector of unique agronomic trait / maltq trait locations
  agro_locs <- reduce(map(GL[df$trait %in% agro_traits], rownames), union)
  core_agro_locs <- reduce(map(GL[df$trait %in% agro_traits], rownames), intersect)
  maltq_locs <- intersect(agro_locs, reduce(map(GL[df$trait %in% quality_traits], rownames), union))
  core_maltq_locs <- intersect(agro_locs, reduce(map(GL[df$trait %in% quality_traits], rownames), intersect))
  
  # Vector of all possible locations
  all_locations <-  sort(union(agro_locs, maltq_locs))
  # Max values of the decision variable for each location (depends on whether malt quality data was collected)
  loc_values <- setNames(object = ifelse(all_locations %in% maltq_locs, 2, 1), nm = all_locations)
  # List of core locations - at least one from each vector must be non-zero
  non_zero_loc_list <- list(core_agro_locs = core_agro_locs, core_maltq_locs = core_maltq_locs)
  
  
  # Implement the genetic algorithm
  optim_out <- GeneticAlg.int(genomeLen = length(all_locations), genomeMin = rep(0, length(loc_values)), genomeMax = loc_values, 
                              popSize = 100, iterations = 150, verbose = FALSE, mutationChance = 0.8, 
                              evalFunc = function(x) -fitness_int(x = x, locs = all_locations, traits = tr_list, G.list = GL1,
                                                                  R.list = RL1, P.list = corL1, M.list = repL1,
                                                                  component.weights = comp_weight, trait.weights = tr_weights,
                                                                  non.zero.env.list = non_zero_loc_list, env.penalty = pen))
                                                                  
  # # For testing
  # x = x; locs = all_locations; traits = tr_list; G.list = GL
  # R.list = RL; P.list = corL; M.list = repL; varY.scaling = varYscaling
  # component.weights = comp_weight; trait.weights = tr_weights
  # non.zero.env.list = non_zero_loc_list; env.penalty = pen; return.all = FALSE
  
  
  
  # Return the fitness components
  optim_x <- optim_out$best$genome
  agro_optim_loc <- all_locations[optim_x != 0]
  maltq_optim_locs <- all_locations[optim_x == 2]
  
  # Calculate scaled fitness and unscaled fitness
  fit_out <- fitness_int(x = optim_x, locs = all_locations, traits = tr_list, G.list = GL, R.list = RL, 
                         P.list = corL, M.list = repL, return.all = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("trait")
  
  fit_out_scaled <- fitness_int(x = optim_x, locs = all_locations, traits = tr_list, G.list = GL1, R.list = RL1, 
                                P.list = corL1, M.list = repL1, return.all = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("trait")
  
  # Scale and average the fitness components
  fit_out_summarize <- fit_out %>%
    mutate(varY = 1 - (varY / df$varY_scaling)) %>%
    summarize_at(vars(-trait), weighted.mean, w = tr_weights, na.rm = TRUE)
  
  # Summarize 
  fit_out_scaled_summarize <- fit_out_scaled %>%
    summarize_at(vars(-trait), weighted.mean, w = tr_weights, na.rm = TRUE)
  
  # Create a tibble that summarizes the optimization results
  optim_results <- tibble(method = "optimization", fitness = list(fit_out), fitness_scaled = list(fit_out_scaled), 
                          fitness_summarize = list(fit_out_summarize), fitness_scaled_summarize = list(fit_out_scaled_summarize),
                          optim_loc = list(list(agro_loc = agro_optim_loc, maltq_loc = maltq_optim_locs)))
  
  # Save this to the list
  sensitivity_param_df$optimization[[i]] <- optim_results
  
  pb$tick()
  
}


sensitivity_test_output <- sensitivity_param_df %>% 
  unnest()



## Model the results of the sensitivity test

# Fit a model to estimate the effect of weights
sensitivity_test_fitness_summary <- sensitivity_test_output %>%
  unnest(fitness_scaled_summarize) %>%
  mutate(varY = 1 - varY) %>%
  gather(metric, metric_value, varY, repAvg, reprAvg)

sensitivity_test_fitness_summary_fitted_models <- sensitivity_test_fitness_summary %>%
  group_by(metric, nursery) %>%
  do(fit = lm(metric_value ~ varY_weight + repAvg_weight + reprAvg_weight, data = .)) %>%
  ungroup()

# Effect plots
sensitivity_test_fitness_summary_fitted_models_effects <- sensitivity_test_fitness_summary_fitted_models %>%
  mutate(effects_df = list(NULL))

for (i in seq_len(nrow(sensitivity_test_fitness_summary_fitted_models_effects))) {
  mod <- sensitivity_test_fitness_summary_fitted_models_effects$fit[[i]]
  effects_df <- effects::allEffects(mod) %>%
    as.data.frame() %>%
    imap_dfr(~rename(.x, weight = 1) %>% mutate(metricWeight = .y)) %>%
    as_tibble() %>%
    left_join(., rename(broom::tidy(summary(mod)), metricWeight = term), by = "metricWeight")
  # Add to the df
  sensitivity_test_fitness_summary_fitted_models_effects$effects_df[[i]] <- effects_df

}


sensitivity_test_fitness_summary_fitted_models_effects1 <- sensitivity_test_fitness_summary_fitted_models_effects %>%
  select(-fit) %>%
  unnest(effects_df) %>%
  mutate(signif_asterick = case_when(p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", TRUE ~ ""),
         annotation = paste0("beta[", metricWeight, "]==", format_numbers(estimate, 2), "*'", signif_asterick, "'"))



# Repeat the same for each trait

# Fit a model to estimate the effect of weights
sensitivity_test_fitness <- sensitivity_test_output %>%
  unnest(fitness_scaled) %>%
  mutate(varY = 1 - varY) %>%
  gather(metric, metric_value, varY, repAvg, reprAvg)

sensitivity_test_fitness_fitted_models <- sensitivity_test_fitness %>%
  group_by(metric, nursery, trait) %>%
  do(fit = lm(metric_value ~ varY_weight + repAvg_weight + reprAvg_weight, data = .)) %>%
  ungroup()

# Effect plots
sensitivity_test_fitness_fitted_models_effects <- sensitivity_test_fitness_fitted_models %>%
  mutate(effects_df = list(NULL))

for (i in seq_len(nrow(sensitivity_test_fitness_fitted_models_effects))) {
  mod <- sensitivity_test_fitness_fitted_models_effects$fit[[i]]
  effects_df <- effects::allEffects(mod) %>%
    as.data.frame() %>%
    imap_dfr(~rename(.x, weight = 1) %>% mutate(metricWeight = .y)) %>%
    as_tibble() %>%
    left_join(., rename(broom::tidy(summary(mod)), metricWeight = term), by = "metricWeight")
  # Add to the df
  sensitivity_test_fitness_fitted_models_effects$effects_df[[i]] <- effects_df

}

sensitivity_test_fitness_fitted_models_effects1 <- sensitivity_test_fitness_fitted_models_effects %>%
  select(-fit) %>%
  unnest(effects_df) %>%
  mutate(signif_asterick = case_when(p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", TRUE ~ ""),
         annotation = paste0("beta[", metricWeight, "]==", format_numbers(estimate, 2), "*'", signif_asterick, "'"))



# Save the results
save("sensitivity_test_output", "sensitivity_test_fitness_fitted_models_effects1",
     "sensitivity_test_fitness_summary_fitted_models_effects1",
     file = file.path(result_dir, "sensitivity_test_results.RData"))




