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


## Unnest the predicted phenotypic data
pheno_blup <- met_mm_out_example %>%
  unnest(pred_pheno) %>%
  filter_at(vars(contains("nonzero")), all_vars(.)) %>%
  left_join(., select(trial_metadata, environment, trial)) %>%
  mutate_at(vars(line_name, trial), as.factor) 


#######################
# Visualize phenotypic balance
#######################


## Visualize balance of lines in trials
ge_tab_pred <- pheno_blup %>% 
  # Complete trial-line_name cases
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
# Mega-environment analysis
#######################


# Remove line names with any missing data
pheno_blup_tomodel <- pheno_blup %>%
  select(trait, trial, environment, line_name, ppv) %>%
  complete(trial, line_name, trait, fill = list(ppv = NA)) %>%
  group_by(trait, line_name) %>%
  filter(all(!is.na(ppv))) %>%
  ungroup()

## Fit a GGE model per trait ##
bilinear_fit <- pheno_blup_tomodel %>%
  group_by(trait) %>%
  do(fit = {
    df <- .
    df1 <- droplevels(df)
    
    # Fit the model
    fit <- bilinear(x = df1, G = "line_name", E = "environment", y = "ppv", test = "bootstrap", B = 1, 
                    model = "AMMI")
    
    ## Return the fit
    fit[c("model", "mu", "Eeffect", "Geffect", "scores", "svdE", "ANOVA")]
    
  }) %>% ungroup()



## Extract genotype and environment scores
bilinear_scores <- bilinear_fit %>% 
  mutate(anova = map(fit, "ANOVA") %>% map(~rownames_to_column(., "term")),
         # Calculate variance explained by PCs using eigenvalues or sums of squares
         varprop = map(fit, "svdE") %>% map(~tibble(PC = paste0("PC", seq_along(.$d)), eigenvalue = .$d, 
                                                         propvar_eigen = eigenvalue / sum(eigenvalue))),
         varprop = map2(varprop, anova, left_join, by = c("PC" = "term")),
         varprop = map(varprop, ~mutate(.x, propvar_SS = SS / sum(SS, na.rm = TRUE), PC_num = parse_number(PC)) %>% 
                         select(-Df, -MS, -testStatistic)))


## Display proportion of variance explained using normalized eigenvalues
bilinear_scores %>%
  unnest(varprop) %>% 
  # Plot
  ggplot(aes(x = PC_num, y = propvar_SS)) +
  geom_point() + 
  facet_wrap(~ trait, scales = "free_x")

## Determine significant PCs using "elbow" method
# Tolerance for difference in variance explained
tol <- 0.03

bilinear_sig_PCs <- bilinear_scores %>%
  unnest(varprop) %>%
  #
  mutate(propvar = propvar_SS) %>%
  #
  arrange(trait, PC_num) %>%
  split(.$trait) %>% 
  ## Calculate the difference between steps of adding PCs. Find the first step when the difference is
  ## below the tolerance threshold
  map_df(~mutate(., propvar_diff = c(abs(diff(propvar)), 0), 
                 stop = which.min(propvar_diff >= tol), 
                 nPC = stop - 1))

## Summary df of number of sig PCs
bilinear_sig_PCs_summ <- bilinear_sig_PCs %>%
  group_by(trait) %>%
  filter(PC_num %in% seq(1, unique(nPC))) %>%
  summarize(total_propvar = sum(propvar), nPC = unique(nPC))


## Fit a model for that number of PCS
bilinear_fitN_fit <- bilinear_fit %>%
  left_join(., bilinear_sig_PCs_summ, by = "trait") %>%
  group_by(trait) %>%
  do({
    
    row <- .
    nPC <- row$nPC
    # Get the fitted ammi model
    fitted <- row$fit[[1]]
    
    # Get the environmental and genotypic effects
    g_effects <- fitted$Geffect %>%
      tibble(line_name = names(.), effect = .)
    e_effects <- fitted$Eeffect %>%
      tibble(environment = names(.), effect = .)
    
    # Get the environmental and genotypic scores
    g_scores <- fitted$scores$Gscores
    e_scores <- fitted$scores$Escores
    
    ## Sum the first nPC scores
    g_scores_sum <- rowSums(g_scores[,seq_len(nPC), drop = FALSE])
    e_scores_sum <- rowSums(e_scores[,seq_len(nPC), drop = FALSE])
    
    ## Combine into DF
    g_effects_scores <- cbind(g_effects, score = g_scores_sum, g_scores)
    e_effects_scores <- cbind(e_effects, score = e_scores_sum, e_scores)
    
    ## Predict y
    # First sum effects
    ge_effect_summ <- outer(
      X = g_effects_scores[,"effect"], 
      Y = e_effects_scores[,"effect"], 
      FUN = "+")
    
    ## Calculate phi, the fitted GxE effect vector using the principal components
    phi <- outer(X = g_scores_sum, Y = e_scores_sum)
    
    # Add the effects to phi, add intercept
    y_hat_mat <- c(fitted$mu) + ge_effect_summ + phi
    # Convert to df
    y_hat_df <- as.data.frame(y_hat_mat) %>%
      rownames_to_column("line_name") %>%
      gather(environment, y_hat, -line_name)
    
    
    ## Return tibble
    tibble(mu = fitted$mu, y_hat = list(y_hat_df), g_scores = list(g_effects_scores),
           e_scores = list(e_effects_scores), phi = list(phi))
    
  }) %>% ungroup()



## Determine mega-environments by looking at the winning genotypes per environment
# Set the direction of "best" traits
trait_dir_df <- tribble(
  ~trait, ~dir,
  "GrainYield", "high",
  "GrainProtein", "low"
)

bilinear_ranks <- bilinear_fitN_fit %>%
  unnest(y_hat) %>%
  left_join(., trait_dir_df) %>%
  split(.$trait) %>%
  map_df(~{
    
    # Rank genotypes by environment
    if (unique(.x$dir) == "high") {
      group_by(.x, trait, environment) %>%
        mutate(rank = desc(row_number(y_hat)),
               rank = rank - min(rank) + 1)
      
    } else {
      group_by(.x, trait, environment) %>%
        mutate(rank = row_number(y_hat))
      
    } 
  }) %>% ungroup()

## Assign mega-environments
bilinear_me <- bilinear_ranks %>% 
  split(.$trait) %>% 
  map_df(~group_by(.x, environment) %>% 
           filter(rank == min(rank)) %>% 
           ungroup() %>%
           mutate(me = as.numeric(as.factor(line_name)))) %>%
  select(-mu, -y_hat, -dir, -rank)




## Save the model output
save("bilinear_fit", "bilinear_fitN_fit", "bilinear_me", file = file.path(result_dir, "bilinear_model_fit.RData"))





#######################
# Model phenotypes
#######################


## Calculate correlations based on phenotypic data
location_gen_cor_pheno <- inner_join(pheno_dat, distinct(pheno_blup_tomodel, trait, trial, line_name)) %>%
  group_by(trait, line_name, location, environment) %>%
  summarize_at(vars(value), mean) %>%
  group_by(trait, location) %>%
  filter(n_distinct(environment) > 1) %>%
  do({
    df <- .
    pairwise_cor <- select(df, line_name, environment, value) %>% 
      spread(environment, value) %>% 
      as.data.frame() %>% 
      column_to_rownames("line_name") %>% 
      cor(use = "pairwise.complete.obs") %>% 
      as.dist() %>%
      tidy() %>% 
      rename_all(~c("environment1", "environment2", "corr")) %>%
      mutate_if(is.factor, as.character)
    
    # Calculate overlap in genotypes between environments
    pairwise_cor %>%
      mutate(overlap = map2_dbl(environment1, environment2, ~filter(df, environment %in% c(.x, .y)) %>% group_by(line_name) %>%
                                  filter(n_distinct(environment) == 2) %>% pull(line_name) %>% 
                                  n_distinct() ))
    
  }) %>%
  ungroup()

## Filter environmental correlations for those that are adjacent in years
location_adj_year_cor <- location_gen_cor_pheno %>%
  left_join(., select(trial_metadata, environment, year1 = year), by = c("environment1" = "environment")) %>%
  left_join(., select(trial_metadata, environment, year2 = year), by = c("environment2" = "environment")) %>%
  filter(abs(year1 - year2) == 1) %>%
  group_by(trait, location) %>%
  mutate(nCorr = n()) %>%
  ungroup()

## Summarize
location_adj_year_cor_summ <- location_adj_year_cor %>%
  group_by(trait, location) %>%
  summarize_at(vars(corr, nCorr), mean)


# Plot
location_adj_year_cor %>%
  ggplot(aes(x = location, y = corr)) +
  geom_jitter(width = 0.25) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(~ trait)








