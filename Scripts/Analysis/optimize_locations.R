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
library(GGally)
library(paletteer)
library(ggrepel)
library(cowplot)
library(gramEvol)


# Color pallete for environments and line names
term_colors <- set_names(paletteer_d("wesanderson::Zissou1")[c(1,5)], "line_name", "environment")

# Function to rename terms
f_term_replace <- function(x) str_replace_all(x, c("line_name" = "Genotype", "environment" = "Environment"))

# A list of traits to analyze
traits <- unique(pheno_dat$trait)

# Load the location analysis data
load(file.path(result_dir, "location_analysis_results.RData"))



# Optimize selection of locations -----------------------------------------

# Optimize the selection of locations for all traits per nursery
# 
# For each location penalty level, select a random set of locations and measure
# fitness components
# 
# Number of resamples
nSamples <- 50

# Vector of location penalties
loc_penalties <- c(0, 0.01, 0.05, 0.1, 0.5)

# List of trait weights
trait_weights <- list(
  equal = setNames(rep(1, length(traits_keep)), traits_keep),
  weighted = c('BetaGlucan' = 1.0, 'DiastaticPower' = 0.5, 'FreeAminoNitrogen' = 0.25, 'GrainProtein' = 1.0, 
               'GrainYield' = 0.75, 'HeadingDate' = 0.25, 'MaltExtract' = 1.0, 'PlantHeight' = 0.25, 
               'PlumpGrain' = 1.0, 'SolubleProteinTotalProtein' = 0.5, 'TestWeight' = 0.25)
)

# List of component weights
comp_weights <- list(
  equal = c(varY = 1, repAvg = 1, reprAvg = 1),
  weighted = c(varY = 1, repAvg = 0.75, reprAvg = 0.25)
)


# Run the optimization procedure 
optimized_locations_all_traits <- location_opimization_input %>%
  crossing(loc_penalty = loc_penalties, tibble(trait_weight_group = names(trait_weights), trait_weights),
           tibble(comp_weight_group = names(comp_weights), comp_weights), scale_components = c(TRUE, FALSE)) %>%
  arrange(nursery, management, loc_penalty, comp_weight_group, trait_weight_group, scale_components) %>%
  group_by(nursery, management, loc_penalty, comp_weight_group, trait_weight_group, scale_components) %>% 
  # Add indices
  cbind(., group_id = group_indices(.)) %>%
  do({
    df <- .
    print(unique(df$group_id))
    
    # Location penalty
    pen <- df$loc_penalty[[1]]
    # Fitness component weights
    comp_weight <- df$comp_weights[[1]]
    # List of traits
    tr_list <- df$trait
    # Trait weights
    tr_weights <- df$trait_weights[[1]]
    
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
    
    # Switch to scale components or not
    # 
    # First, do not scale
    if (!all(df$scale_components)) {

      # Adjust the varY scaling
      varYscaling <- df$varY_scaling
      
      optim_out <- GeneticAlg.int(genomeLen = length(all_locations), genomeMin = rep(0, length(loc_values)), genomeMax = loc_values, 
                                  popSize = 100, iterations = 150, verbose = FALSE, mutationChance = 0.8, 
                                  evalFunc = function(x) -fitness_int(x = x, locs = all_locations, traits = tr_list, G.list = GL,
                                                                      R.list = RL, P.list = corL, M.list = repL, varY.scaling = varYscaling,
                                                                      component.weights = comp_weight, trait.weights = tr_weights,
                                                                      non.zero.env.list = non_zero_loc_list, env.penalty = pen))
      
    # Yes scale
    } else {
      
      # Adjust the varY scaling
      varYscaling <- 1
      
      optim_out <- GeneticAlg.int(genomeLen = length(all_locations), genomeMin = rep(0, length(loc_values)), genomeMax = loc_values, 
                                   popSize = 100, iterations = 150, verbose = FALSE, mutationChance = 0.8, 
                                   evalFunc = function(x) -fitness_int(x = x, locs = all_locations, traits = tr_list, G.list = GL1,
                                                                       R.list = RL1, P.list = corL1, M.list = repL1, varY.scaling = varYscaling,
                                                                       component.weights = comp_weight, trait.weights = tr_weights,
                                                                       non.zero.env.list = non_zero_loc_list, env.penalty = pen))

      
    }
    
    
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
    
    ## Now randomly sample the same number of selelcted locations and calculate fitness
    # Determine the number of agronomic locations without maltq
    n_just_agro_locs <- length(setdiff(agro_optim_loc, maltq_optim_locs))
    n_agro_maltq_locs <- length(maltq_optim_locs)
    
    # Generate samples
    blank_optim_loc <- rep(0, length(optim_x))
    random_loc_samples <- replicate(n = nSamples, simplify = FALSE, expr = {
      # First sample locations for malting quality
      maltq_loc_sample <- sample(x = which(loc_values == 2), size = n_agro_maltq_locs)
      # Then sample for agronomic traits from the remainder
      just_agro_loc_sample <- sample(x = setdiff(seq_along(loc_values), maltq_loc_sample), size = n_just_agro_locs)
      
      blank_optim_loc[just_agro_loc_sample] <- 1
      blank_optim_loc[maltq_loc_sample] <- 2
      blank_optim_loc
    })
    
    # Calculate fitness for the random samples
    random_loc_fitness <- random_loc_samples %>%
      map(~fitness_int(x = .x, locs = all_locations, traits = tr_list, G.list = GL, R.list = RL, P.list = corL, M.list = repL, 
                       return.all = TRUE) %>% as.data.frame() %>% rownames_to_column("trait") )
    
    random_loc_fitness_scaled <- random_loc_samples %>%
      map(~fitness_int(x = .x, locs = all_locations, traits = tr_list, G.list = GL1, R.list = RL1, P.list = corL1, M.list = repL1, 
                       return.all = TRUE) %>% as.data.frame() %>% rownames_to_column("trait") )
    
    # Summarize
    random_loc_fitness_summarize <- random_loc_fitness %>%
      map_df(~{
        mutate(.x, varY = 1 - (varY / df$varY_scaling)) %>%
          summarize_at(vars(-trait), weighted.mean, w = tr_weights, na.rm = TRUE)
      })

    random_loc_fitness_scaled_summarize <- random_loc_fitness_scaled %>%
      map_df(~summarize_at(.x, vars(-trait), weighted.mean, w = tr_weights, na.rm = TRUE))
    
    
    # Save random fitness results
    random_results <- tibble(method = "random", fitness = list(random_loc_fitness), fitness_scaled = list(random_loc_fitness_scaled), 
                             fitness_summarize = list(random_loc_fitness_summarize), fitness_scaled_summarize = list(random_loc_fitness_scaled_summarize))
    
    
    ## Calculate the fitness of the ranked choice environments ##
    best_rank_locations <- subset(nursery_best_rank_locations, nursery %in% df$nursery, location, drop = TRUE)
    # Convert this to 0, 1, 2
    best_rank_x <- loc_values
    best_rank_x[!names(best_rank_x) %in% best_rank_locations] <- 0
    
    best_rank_fitness <- fitness_int(x = best_rank_x, locs = all_locations, traits = tr_list, G.list = GL, R.list = RL, 
                           P.list = corL, M.list = repL, return.all = TRUE) %>%
      as.data.frame() %>%
      rownames_to_column("trait")
    
    best_rank_fitness_scaled <- fitness_int(x = best_rank_x, locs = all_locations, traits = tr_list, G.list = GL1, R.list = RL1, 
                                  P.list = corL1, M.list = repL1, return.all = TRUE) %>%
      as.data.frame() %>%
      rownames_to_column("trait")
    
    # Scale and average the fitness components
    best_rank_fitness_summarize <- best_rank_fitness %>%
      mutate(varY = 1 - (varY / df$varY_scaling)) %>%
      summarize_at(vars(-trait), weighted.mean, w = tr_weights, na.rm = TRUE)
    
    # Summarize 
    best_rank_fitness_scaled_summarize <- best_rank_fitness_scaled %>%
      summarize_at(vars(-trait), weighted.mean, w = tr_weights, na.rm = TRUE)
    
    # Save ranked locations fitness results
    best_rank_results <- tibble(method = "best_rank", fitness = list(best_rank_fitness), fitness_scaled = list(best_rank_fitness_scaled), 
                                fitness_summarize = list(best_rank_fitness_summarize), fitness_scaled_summarize = list(best_rank_fitness_scaled_summarize),
                                optim_loc = list(list(agro_loc = all_locations[best_rank_x != 0],  maltq_loc = all_locations[best_rank_x == 2])))
    
    # Combine and return results
    bind_rows(optim_results, random_results, best_rank_results)
    
  }) %>% ungroup()


## Save the optimization steps from an example
set.seed(1023)

# Run the optimization procedure 
optimized_locations_example <- location_opimization_input %>%
  crossing(loc_penalty = loc_penalties, tibble(trait_weight_group = names(trait_weights), trait_weights),
           tibble(comp_weight_group = names(comp_weights), comp_weights), scale_components = c(TRUE, FALSE)) %>%
  filter(loc_penalty == 0.01, comp_weight_group == "equal", trait_weight_group == "equal", scale_components) %>%
  group_by(nursery, management, loc_penalty, comp_weight_group, trait_weight_group, scale_components) %>% 
  do({
    df <- .

    # Location penalty
    pen <- df$loc_penalty[[1]]
    # Fitness component weights
    comp_weight <- df$comp_weights[[1]]
    # List of traits
    tr_list <- df$trait
    # Trait weights
    tr_weights <- df$trait_weights[[1]]
    
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
    # Adjust the varY scaling
    varYscaling <- 1
    # Empty list to store the best population at each generation
    best_pop_gen <- list()
    
    # Monitor function
    monitor <- function(optim) {
      best_pop_gen <- append(best_pop_gen, values = list(optim$best$genome))
      assign(x = "best_pop_gen", value = best_pop_gen, envir = .GlobalEnv)
    }
    
    
    optim_out <- GeneticAlg.int(genomeLen = length(all_locations), genomeMin = rep(0, length(loc_values)), genomeMax = loc_values, 
                                popSize = 100, iterations = 250, verbose = FALSE, mutationChance = 0.8, 
                                evalFunc = function(x) -fitness_int(x = x, locs = all_locations, traits = tr_list, G.list = GL1,
                                                                    R.list = RL1, P.list = corL1, M.list = repL1, varY.scaling = varYscaling,
                                                                    component.weights = comp_weight, trait.weights = tr_weights,
                                                                    non.zero.env.list = non_zero_loc_list, env.penalty = pen), 
                                monitorFunc = monitor) 
                                
    
    # Determine the selected locations for each of the best genomes
    selected_locations_gen <- map(best_pop_gen, ~{
      list(agro_optim_loc = all_locations[.x != 0],
           maltq_optim_locs = all_locations[.x == 2])
    }) %>% 
      # Convert to df
      transpose() %>%
      as_tibble() %>%
      mutate(iter = seq_len(nrow(.)))
    
    # Calculate the fitness components
    
    # Calculate scaled fitness and unscaled fitness
    fit_out <- map(best_pop_gen, ~{
      fitness_int(x = .x, locs = all_locations, traits = tr_list, G.list = GL, R.list = RL, 
                  P.list = corL, M.list = repL, return.all = TRUE) %>%
        as.data.frame() %>%
        rownames_to_column("trait")
    })
    
    fit_out_scaled <- map(best_pop_gen, ~{
      fitness_int(x = .x, locs = all_locations, traits = tr_list, G.list = GL1, R.list = RL1, 
                  P.list = corL1, M.list = repL1, return.all = TRUE) %>%
        as.data.frame() %>%
        rownames_to_column("trait")
    })
    
    # Scale and average the fitness components
    fit_out_summarize <- fit_out %>%
      map(~{
        mutate(.x, varY = 1 - (varY / df$varY_scaling)) %>% 
          summarize_at(vars(-trait), weighted.mean, w = tr_weights, na.rm = TRUE)
      }) %>%
      imap_dfr(~mutate(.x, iter = .y))
    
    # Summarize 
    fit_out_scaled_summarize <- fit_out_scaled %>%
      map(~{
        mutate(.x, varY = 1 - (varY / df$varY_scaling)) %>% 
          summarize_at(vars(-trait), weighted.mean, w = tr_weights, na.rm = TRUE)
      }) %>%
      imap_dfr(~mutate(.x, iter = .y))
    
    # Create a tibble that summarizes the optimization results
    optim_results <- tibble(method = "optimization", fitness = list(fit_out), fitness_scaled = list(fit_out_scaled), 
                            fitness_summarize = list(fit_out_summarize), fitness_scaled_summarize = list(fit_out_scaled_summarize),
                            selected_locations_gen = list(selected_locations_gen))
    
  }) %>% ungroup()


# Save everything
save("optimized_locations_all_traits", "optimized_locations_example", 
     file = file.path(result_dir, "location_optimization_results.RData"))


