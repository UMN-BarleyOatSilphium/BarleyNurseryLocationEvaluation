## Barley Nursery Analysis
## 
## Analyze MET from barley nursery
## 
## Author: Jeff Neyhart
## Last modified: 26 Feb. 2020
## 



# Base script
proj_dir <- getwd()
source(file.path(proj_dir, "startup.R"))

# Additional packages
library(GGally)
library(Bilinear)





#######################
# Load data
#######################

load(file.path(result_dir, "met_mixed_model_output_example.RData"))




#######################
# Visualize phenotypic balance
#######################


## Unnest the predicted phenotypic data
pheno_blup <- met_mm_out_example %>%
  unnest(pred_pheno) %>%
  filter_at(vars(contains("nonzero")), all_vars(.))


## Visualize balance of lines in trials
ge_tab_pred <- pheno_blup %>% 
  left_join(., select(trial_metadata, environment, trial)) %>%
  # Complete trial-line_name cases
  mutate_at(vars(line_name, trial), as.factor) %>%
  complete(line_name, trial, trait, fill = list(pheno = NA, ppv = NA)) %>%
  mutate(obs = case_when(
    !is.na(pheno) ~ "phenotype",
    !is.na(ppv) ~ "prediction",
    TRUE ~ as.character(NA)
  )) %>%
  group_by(line_name, trait) %>%
  mutate(nPheno = sum(obs == "phenotype", na.rm = T) / n_distinct(trial),
         nMissing = sum(is.na(obs)) / n_distinct(trial)) %>%
  ungroup() %>% 
  arrange(trait, trial, desc(nPheno), nMissing) %>%
  mutate_at(vars(line_name, trial), fct_inorder)



freq_df_plot <- ge_tab_pred %>%
  # Convert line_name and location to numeric
  mutate_at(vars(line_name, trial), list(num = ~fct_inseq(as.factor(as.numeric(.)))))

## Plot tiles of years of genotype-location observations
g_pred_freq <- freq_df_plot %>%
  filter(trait == first(trait)) %>%
  ggplot(aes(x = trial_num, y = line_name_num, fill = obs)) +
  geom_tile() +
  # facet_grid(~ trait, drop = TRUE, scales = "free_x", space = "free_x") +
  scale_y_discrete(breaks = function(x) as.character(pretty(as.numeric(x), n = 20)), name = "Genotype") +
  scale_x_discrete(breaks = function(x) as.character(pretty(as.numeric(x), n = 20)), name = "Trial") +
  scale_fill_viridis_d(na.value = "white", begin = 0.2, end = 0.8,
                       guide = guide_legend(override.aes = list(color = "black"))) +
  theme_presentation2(base_size = 16)


# Save
ggsave(filename = "genotype_trial_contingency_example.jpg", path = fig_dir, plot = g_pred_freq,
       height = 6, width = 8, dpi = 500)

## Plot tiles of years of genotype-location observations
g_pred_freq <- freq_df_plot %>%
  filter(trait == first(trait), obs == "phenotype") %>%
  ggplot(aes(x = trial_num, y = line_name_num, fill = obs)) +
  geom_tile() +
  # facet_grid(~ trait, drop = TRUE, scales = "free_x", space = "free_x") +
  scale_y_discrete(breaks = function(x) as.character(pretty(as.numeric(x), n = 20)), name = "Genotype") +
  scale_x_discrete(breaks = function(x) as.character(pretty(as.numeric(x), n = 20)), name = "Trial") +
  scale_fill_viridis_d(na.value = "white", begin = 0.2, end = 0.8,
                       guide = guide_legend(override.aes = list(color = "black"))) +
  theme_presentation2(base_size = 16)


# Save
ggsave(filename = "genotype_trial_contingency_example1.jpg", path = fig_dir, plot = g_pred_freq,
       height = 6, width = 8, dpi = 500)





## Table of phenotypic observations, predictions, and still missings
freq_df_plot %>%
  group_by(trait, obs) %>% 
  summarize(n = n()) %>% 
  spread(obs, n)





#######################
# Summarize error
#######################

## Summarize error variance by locations
met_var_comp <- met_mm_out_example %>% 
  unnest(var_comp) %>%
  filter(str_detect(component, "units")) %>%
  rename(environment = component) %>%
  mutate(environment = str_remove_all(environment, ":units"))

# Bind with trial metadata
met_var_comp1 <- met_var_comp %>%
  left_join(trial_metadata) %>%
  filter(variance != 0) %>%
  mutate(sdev = sqrt(variance))

## Number of observations per location/trait
met_var_comp_plot <- met_var_comp1 %>% 
  group_by(trait, location) %>% 
  summarize(nE = n_distinct(environment)) %>%
  left_join(met_var_comp1, .)

# Plot error variance by location and trait
g_met_varR <- met_var_comp1 %>%
  ggplot(aes(location, y = variance, fill = trait)) +
  geom_boxplot(position = "dodge") +
  facet_grid(trait ~ ., scales = "free_y") +
  theme_light()

# Save
ggsave(filename = "nursery_example_location_varR.jpg", path = fig_dir, plot = g_met_varR,
       height = 6, width = 8, dpi = 500)



#######################
# Model phenotypes
#######################




# Remove line names with any missing data
pheno_blup_tomodel <- pheno_blup %>%
  split(.$trait) %>%
  map_df(~group_by(.x, line_name) %>% filter(n() == n_distinct(.x$environment))) %>%
  ungroup() %>%
  left_join(trial_metadata)

## Convert to list of GE tables
pheno_blup_ge_list <- pheno_blup_tomodel %>%
  select(trait, environment, line_name, ppv) %>%
  as.data.frame() %>%
  split(.$trait) %>%
  map(~select(., -trait) %>% spread(environment, ppv) %>% column_to_rownames("line_name")  %>% as.matrix())


## Fit an AMMI model per trait ##



ammi_fit <- bilinear(x = subset(pheno_blup_tomodel, trait == "GrainYield"), G = "line_name", E = "environment",
                     y = "ppv", test = "bootstrap", B = 2)

# Plot
AMMIplot(bilinearObject = ammi_fit)

BBplot(bilinearObject = ammi_fit, Gnames = FALSE, nPC = 2, decorateGGE = T)



## Calculate the correlation between environments at a location
location_gen_cor <- pheno_blup_tomodel %>%
  group_by(trait, location) %>%
  filter(n_distinct(environment) > 1) %>%
  do(loc_cor_out = {
    df <- .
    select(df, line_name, environment, ppv) %>% 
      spread(environment, ppv) %>% 
      as.data.frame() %>% 
      column_to_rownames("line_name") %>% 
      cor() %>% 
      as.dist() %>%
      tidy() %>% 
      rename_all(~c("environment1", "environment2", "corr"))
  }) %>%
  ungroup() %>%
  mutate(mean_corr = map_dbl(loc_cor_out, ~mean(.$corr)))

## Calculate correlations based on phenotypic data
location_gen_cor_pheno <- inner_join(pheno_dat, distinct(pheno_blup_tomodel, trait, trial, line_name)) %>%
  group_by(trait, line_name, location, environment) %>%
  summarize_at(vars(value), mean) %>%
  group_by(trait, location) %>%
  filter(n_distinct(environment) > 1) %>%
  do(loc_cor_out = {
    df <- .
    select(df, line_name, environment, value) %>% 
      spread(environment, value) %>% 
      as.data.frame() %>% 
      column_to_rownames("line_name") %>% 
      cor(use = "pairwise.complete.obs") %>% 
      as.dist() %>%
      tidy() %>% 
      rename_all(~c("environment1", "environment2", "corr"))
  }) %>%
  ungroup() %>%
  mutate(mean_corr = map_dbl(loc_cor_out, ~mean(.$corr)))
  








