## Barley Nursery Analysis
## 
## Analyze stability from MET BLUPs
## 
## Author: Jeff Neyhart
## 


# Base script
proj_dir <- getwd()
source(file.path(proj_dir, "startup.R"))

# Additional packages
library(GGally)
library(paletteer)
library(ggrepel)
library(cowplot)
library(lme4)
library(modelr)



# Color pallete for environments and line names
term_colors <- set_names(paletteer_d(package = "wesanderson", palette = "Zissou1")[c(1,5)], "line_name", "environment")

# Function to rename terms
f_term_replace <- function(x) str_replace_all(x, c("line_name" = "Genotype", "environment" = "Environment"))




#######################
# Load data
#######################

# List data files
file_list <- list.files(path = result_dir, pattern = "met_mixed_model_output.RData", full.names = TRUE)
# Load the data and return the object
data_list <- lapply(X = file_list, FUN = function(f) { load(f); get("met_mm_out") })

# Bind
met_mm_out <- bind_rows(data_list) %>%
  ungroup() %>%
  mutate(out = map(out, ~mutate(., H2 = unlist(H2)))) %>%
  unnest(out) %>%
  filter(! map_lgl(pred_pheno, is.null))

## Unnest the predicted phenotypic data
pheno_blup <- met_mm_out %>%
  unnest(pred_pheno) %>%
  # filter_at(vars(contains("nonzero")), all_vars(.)) %>%
  # Filter for non-zero genotypic main effect predictions
  filter(pgv_g_nonzero) %>%
  left_join(., select(trial_metadata, environment, trial)) %>%
  mutate_at(vars(line_name, trial, environment), as.factor)

pheno_dat_tomodel <- pheno_dat %>%
  select(-trial) %>%
  left_join(., select(trial_metadata, environment, trial, nursery, management)) %>%
  mutate_at(vars(line_name, trial, environment), as.factor)



#######################
# Determine the most appropriate time window for calculating stability
#######################

# First, summarize the number of environments per trait and year
pheno_dat_tomodel %>% 
  group_by(nursery, management, trait, year) %>% 
  summarize(nEnv = n_distinct(environment)) %>%
  summarize_at(vars(nEnv), list(~min, ~max, ~median)) %>%
  as.data.frame()

# A list of traits to analyze
traits_tokeep <- c("GrainProtein", "MaltExtract", "GrainYield", "HeadingDate", "PlantHeight",
                   "BetaGlucan", "DiastaticPower", "FreeAminoNitrogen", "KernelWeight",
                   "PlumpGrain", "SolubleProteinTotalProtein", "TestWeight")

# Data frame of scenarios to keep
scenarios_tokeep <- pheno_dat_tomodel %>% 
  filter(trait %in% traits_tokeep) %>%
  group_by(nursery, management, trait) %>% 
  summarize(nEnv = n_distinct(environment)) %>% 
  ungroup() %>% 
  filter(nEnv > 1)





# For each nursery, management, trait, find number of common lines across various time windows
# Set possible year windows
year_windows <- 2:5

pheno_dat_tomodel_window <- pheno_dat_tomodel %>% 
  inner_join(., scenarios_tokeep) %>%
  group_by(nursery, management, trait) %>%
  do({
    df <- droplevels(.)
    
    # Sequence of years
    year_seq <- sort(unique(df$year1))
    # Generate lists of windows
    window_lists <- year_windows %>%
      set_names(., paste0("window", .)) %>%
      map(~{
        # Vector of start years
        year_start <- seq_range(x = year_seq, by = .x)
        # Remove the last start if it is the last year
        if (last(year_start) == last(year_seq)) year_start <- head(year_start, -1)
        # Vector of end years
        year_end <- c(tail(year_start - 1, -1), last(year_seq))
        
        # Generate year sequence to use
        year_seq_use <- map2(.x = year_start, .y = year_end, ~year_seq[between(year_seq, .x, .y)]) %>% 
          subset(., map_lgl(., ~length(.) > 0))  
        # If any are of length 1, merge with the next sequence
        which_length_one <- which(map_lgl(year_seq_use, ~length(.) == 1))
        if (length(which_length_one) > 0 & any(which_length_one < length(year_seq_use))) {
          year_seq_use1 <- year_seq_use
          year_seq_use1[[which_length_one]] <- unlist(year_seq_use1[c(which_length_one, which_length_one+1)])
          year_seq_use <- year_seq_use1[-(which_length_one + 1)]
        }
        year_seq_use
        
      }) %>% 
      map(~`names<-`(., seq_along(.)))
    
    ## Find the number of environments and common genotypes per window
    window_summary <- window_lists %>%
      map(~map_df(., ~{
        common_geno <- filter(df, year %in% .x) %>%
          droplevels() %>%
          split(.$year) %>% 
          map("line_name")
        
        if (is_empty(common_geno)) {
          common_geno <- NULL
        } else {
          common_geno <- reduce(common_geno, intersect)
        }
        
        
        filter(df, year %in% .x, line_name %in% common_geno) %>%
          droplevels() %>%
          group_by(line_name) %>%
          summarize(nEnv = n_distinct(environment), nYear = n_distinct(year), n = n(),
                    min_year = min(year1), max_year = max(year1))
        
      })) %>% imap_dfr(~mutate(.x, window_size = .y))
    
    ## Range in number of lines and environments per window
    window_summary_range <- window_summary %>%
      group_by(window_size, min_year, max_year) %>%
      summarize(nGeno = n_distinct(line_name), min_env = min(nEnv), max_env = max(nEnv)) %>%
      group_by(window_size) %>%
      summarize(min_geno = min(nGeno), max_geno = max(nGeno), min_env = min(min_env), max_env = max(max_env))
    
    ## Find the common checks across all windows for each window length
    window_common_checks <- window_summary %>% 
      group_by(window_size) %>% 
      mutate(nWindow = n_distinct(min_year)) %>% 
      group_by(window_size, line_name) %>% 
      filter(n() == nWindow) %>% 
      ungroup() %>% 
      distinct(window_size, line_name)
    
    ## Return all of this information
    tibble(summary = list(window_summary), range = list(window_summary_range), common_checks = list(window_common_checks))
    
  }) %>% ungroup()


## Plot number of lines/environments for each subset scenario
pheno_dat_tomodel_window_summary <- pheno_dat_tomodel_window %>%
  unnest(summary) %>%
  group_by(nursery, management, trait, window_size, min_year, max_year) %>%
  summarize(nGeno = n_distinct(line_name), min_nEnv = min(nEnv), max_nEnv = max(nEnv)) %>%
  ungroup()

g <- pheno_dat_tomodel_window_summary %>%
  filter(trait == "GrainYield") %>%
  mutate(min_year = as.factor(min_year)) %>%
  ggplot(aes(x = nGeno, y = min_nEnv, color = min_year)) +
  geom_point() +
  scale_color_discrete(guide = FALSE) +
  facet_grid(window_size ~ nursery + management) +
  theme_presentation2(base_size = 12)
plotly::ggplotly(g)





#######################
# Calculate the stability of each line
#######################

## First fit full models to calculate environmental effects
overall_environmental_means <- pheno_dat_tomodel %>% 
  inner_join(., scenarios_tokeep) %>%
  group_by(nursery, management, trait) %>%
  do({
    df <- droplevels(.)
    
    # Set environment to factor
    df1 <- mutate_at(df, vars(environment, line_name), fct_contr_sum)
    
    # Fit a model with random genotype
    fit_mm <- lmer(value ~ 1 + (1|line_name) + environment, data = df1)
    
    ## Return the environmental effects ##
    # Get the grand mean
    mu <- fixef(fit_mm)[1]
    
    # Get the environmental means and return
    tidy(fit_mm) %>%
      filter(str_detect(term, "environment")) %>% 
      mutate(term = str_remove(term, "environment")) %>% 
      add_row(term = last(levels(df1$environment)), estimate = -sum(.$estimate)) %>% 
      select(environment = term, effect = estimate)
  
  }) %>% ungroup()


# Year window to use
year_window <- "window2"


## Set up data to model
pheno_stability_tomodel <- pheno_dat_tomodel_window %>% 
  unnest(summary) %>% 
  filter(window_size == year_window) %>%
  # Add phenotypic data
  left_join(., pheno_dat_tomodel) %>%
  # Filter for years
  filter(year1 >= min_year & year1 <= max_year) %>%
  # Add environmental means
  left_join(., overall_environmental_means) %>%
  # Create names for subsets
  unite("subset", window_size, min_year, max_year) %>%
  select(nursery, management, trait, subset, environment, line_name, value, effect)
  
  
# Normal FW model
pheno_stability_subsets <- pheno_stability_tomodel %>%
  group_by(nursery, management, trait, subset, line_name) %>%
  do(tibble(dat = list(.), fit = list(lm(value ~ effect, data = .)), nObs = nrow(.))) %>%
  ungroup() %>%
  mutate(geno_mean = map_dbl(fit, ~coef(.)[1]),
         b_i = map_dbl(fit, ~coef(.)[2]),
         d_i = map_dbl(fit, sigma)) %>%
  select(-fit)

## Stability of long-term checks overall
pheno_stability_checks <- pheno_stability_tomodel %>%
  inner_join(., select(filter(unnest(pheno_dat_tomodel_window, common_checks), window_size == year_window), -window_size)) %>%
  group_by(nursery, management, trait, line_name) %>%
  do(tibble(dat = list(.), fit = list(lm(value ~ effect, data = .)), nObs = nrow(.))) %>%
  ungroup() %>%
  mutate(geno_mean = map_dbl(fit, ~coef(.)[1]),
         b_i = map_dbl(fit, ~coef(.)[2]),
         d_i = map_dbl(fit, sigma)) %>%
  select(-fit)


## Model common checks
pheno_stability_checks_tomodel <- pheno_stability_subsets %>%
  inner_join(., select(pheno_stability_checks, nursery:line_name, geno_mean:d_i), 
             by = c("nursery", "management", "trait", "line_name")) %>%
  rename_at(vars(ends_with(".x")), ~str_replace_all(., ".x", "_subset")) %>%
  rename_at(vars(ends_with(".y")), ~str_replace_all(., ".y", "_overall"))

## Test for subset and line_name effects
pheno_stability_checks_fits <- pheno_stability_checks_tomodel %>%
  group_by(nursery, management, trait) %>%
  # Remove instances with less than 2 lines
  filter(n_distinct(line_name) >= 2) %>%
  do({
    df <- .
    b_i_fit <- lm(b_i_subset ~ line_name + subset, df)
    d_i_fit <- lm(d_i_subset ~ line_name + subset, df)
    g_i_fit <- lm(geno_mean_subset ~ line_name + subset, df)
    ## return
    tibble(b_i_fit = list(b_i_fit), d_i_fit = list(d_i_fit), g_i_fit = list(g_i_fit))
  })

## Anovas




## Plot common checks
pheno_stability_checks_plot <- pheno_stability_subsets %>%
  inner_join(., select(pheno_stability_checks, nursery:line_name)) %>%
  unnest(dat) %>%
  select(-ends_with("1"))

pheno_stability_checks_plot %>%
  filter(trait == "GrainYield", nursery == "mvn") %>%
  ggplot(aes(x = effect, y = value, color = subset)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_smooth(aes(group = line_name), method = "lm", se = FALSE) +
  facet_grid(trait ~ line_name, switch = "y") +
  theme_presentation2(12)



## Model stability over time
pheno_stability_subsets_fits <- pheno_stability_subsets %>%
  group_by(nursery, management, trait) %>%
  do({
    df <- .
    b_i_fit <- lm(b_i ~ line_name + subset, df)
    d_i_fit <- lm(d_i ~ line_name + subset, df)
    g_i_fit <- lm(geno_mean ~ line_name + subset, df)  
    tibble(b_i_fit = list(b_i_fit), d_i_fit = list(d_i_fit), g_i_fit = list(g_i_fit))
  }) %>% ungroup()




## Plot stability slopes over subsets
pheno_stability_subsets %>%
  ggplot(aes(x = b_i, fill = subset)) +
  geom_density() +
  scale_fill_discrete(guide = FALSE) +
  facet_wrap(~ trait, scales = "free_x") +
  theme_presentation2(12)












