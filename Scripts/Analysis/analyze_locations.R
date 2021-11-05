## Barley Nursery Analysis
## 
## This script uses the mixed model output to analyze the locations. 
## - Environments are clustered using the FA loadings estimated in the model
## - Calculate genetic correlations between environments using variance components
## - Calculate precision per environment
## - Calculate reliability, repeatability, and representativeness
## - Perform simple ranking of locations based on the above metrics
## - Run the optimization for selecting locations
## - Perform a sensitivity analysis  
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




# Color pallete for environments and line names
term_colors <- set_names(paletteer_d("wesanderson::Zissou1")[c(1,5)], "line_name", "environment")

# Function to rename terms
f_term_replace <- function(x) str_replace_all(x, c("line_name" = "Genotype", "environment" = "Environment"))

# A list of traits to analyze
traits <- unique(pheno_dat$trait)

# Minimum number of environments to include a location
min_env <- 3




# Load asreml model output ------------------------------------------------

load(file.path(result_dir, "phenotypic_analysis_asreml_fa.RData"))

# Remove WRN rainfed
analysis_models_fa <- analysis_models_fa %>%
  filter(!(nursery == "wrn" & management == "rainfed"),
         trait %in%  traits_keep)


## Extract variance components from the models and tidy
var_comp_tidy <- analysis_models_fa %>% 
  mutate(varcomp = map(summary, "varcomp"),
         varcomp = map(varcomp, ~rownames_to_column(., "term"))) %>%
  unnest(varcomp) %>%
  # Get the variance component type
  mutate(type = str_extract(term, "![A-Za-z0-9]*$") %>% str_remove(., "!"),
         type = case_when(type == "var" ~ "varG", type == "R" ~ "varR", TRUE ~ as.character(type)),
         environment = str_extract(term, "[A-Z]{3}[0-9]{2}[A-Z]"),
         # If environment is NA, that means that homogenous variance was used
         environment = ifelse(is.na(environment), "all", environment)) %>%
  # Add location and year data
  left_join(., distinct(trial_metadata, trial, location, year, environment, nursery, management)) %>%
  select(trait, nursery, management, trial, environment, location, year, varP, type, component, std.error) %>%
  arrange(nursery, management, trait, location, year)

## Create matrices of the loadings, residual genetic variance, and residual variance
var_comp_matrices <- var_comp_tidy %>%
  select(trait:management, environment, varP, type, component) %>%
  group_by(trait, nursery, management, varP) %>%
  nest() %>%
  ungroup() %>%
  ## Loadings
  mutate(Lambda = map(data, ~filter(.x, str_detect(type, "fa"))),
         # Residual genetic variance
         Psi = map(data, ~filter(.x, type == "varG")),
         # Residual
         varR = map(data, ~filter(.x, type == "varR"))) %>%
  # Convert all to matrices
  mutate_at(vars(Lambda, Psi, varR), ~map(., ~spread(.x, type, component) %>%  as.data.frame() %>% 
                                            column_to_rownames("environment") %>% as.matrix())) %>%
  mutate(Psi = map(Psi, ~`dimnames<-`(diag(.x[,1], nrow = nrow(.x)), replicate(2, row.names(.x), simplify = FALSE)))) %>%
  ## Calculate the rho matrix
  mutate(Rho = map2(Lambda, Psi, ~tcrossprod(.x) + .y),
         # D matrix
         D = map(Rho, ~`dimnames<-`(x = diag(diag(.x)^(-0.5)), value = dimnames(.x)))) %>%
  select(-data) %>%
  # Standardize the variance estimates
  mutate_at(vars(Psi, varR, Rho), list(scaled = ~map2(., varP, ~.x / .y)))
  


## Calculate the proportion of genetic variance explained by each factor
loadings_prop_var <- var_comp_matrices %>%
  mutate(prop_var = map2(Lambda, Rho, ~{
    # Variance explained by all factors
    LLT <- tcrossprod(.x)
    all_factor_prop_genvar <- sum(diag(LLT)) / sum(diag(LLT + .y))
    
    # Calculate the proportion of total factor variance explained by each factor
    factor_prop <- apply(X = .x, MARGIN = 2, FUN = function(x) {
      xxT <- tcrossprod(x)
      sum(diag(xxT)) / sum(diag(LLT))
    })
    
    # Multiply this proportion by the all_factor_prop_genvar to get the proportion 
    # of total genetic variance explained by each factor
    factor_prop_genvar <- factor_prop * all_factor_prop_genvar
    
    # Variance explained by residual genetic
    resid_prop_genvar <- 1 - all_factor_prop_genvar
    
    # Return a df
    as.data.frame(cbind(t(factor_prop_genvar), resid = resid_prop_genvar))
    
  })) %>% unnest(prop_var) %>%
  select(trait, nursery, management, varP, starts_with("fa"), resid)





# A dataframe of locations with at least min_env environments of observations
candidate_locations <- var_comp_tidy %>% 
  group_by(nursery, management, trait, location) %>% 
  filter(n_distinct(environment) > min_env) %>% 
  ungroup() %>% 
  distinct(nursery, management, trait, location)



## Calculate a correlation matrix among environments
env_cor_matrix <- var_comp_matrices %>%
  mutate(cor_mat = map2(Rho, D, ~.y %*% .x %*% .y)) %>%
  select(trait, nursery, management, cor_mat)


## Create a matrix of the genotypic scores from the blups
gen_score_matrix <- analysis_models_fa %>%
  mutate(geno_scores = map(geno_scores, ~spread(.x, component, estimate) %>% as.data.frame() %>%
                             column_to_rownames("line_name") %>% as.matrix()) ) %>%
  select(trait:management, geno_scores)


# Combine data frames of the loadings/scores
loadings_scores_df <- left_join(var_comp_matrices, gen_score_matrix) %>%
  mutate(loadings = map2(Lambda, geno_scores, rbind),
         loadings = map(loadings, ~as.data.frame(.x) %>% rownames_to_column("term") %>%
                          tbl_df() %>% mutate(type = ifelse(term %in% unique(var_comp_tidy$environment), "environment", "genotype"))) ) %>%
  select(trait, nursery, management, loadings)



## Comparison of pedigree/non-pedigree models using AIC
pedigree_model_comparison <- analysis_models_fa %>%
  select(nursery, trait, contains("aic"))









# Cluster environments using FA loadings ----------------------------------

nursery_mega_environment_output <- loadings_scores_df %>%
  group_by(trait, nursery, management) %>%
  do({
    
    row <- .
    tr <- row$trait
    nur <- row$nursery
    mgmt <- row$management
  
    
    ## Example loadings
    .x <- row$loadings[[1]] %>%
      # Adjust sign of genotype scores based on trait
      mutate_at(vars(fa1, fa2), ~ifelse(tr %in% neg_val_traits & type == "genotype", . * -1, .))
    
    # Adjust loadings using a scalar to make the length of the longest enviornmental
    # vector equal to the length of the longest genotypic vector
    
    # Calculate vector lengths for each (all vectors start at the origin)
    # Then calculate the ratio between the largest length and the smallest length
    # The square root of this value is the scaling parameter
    longest_vectors_df <- .x %>% 
      mutate(len = sqrt((fa1^2) + (fa2^2))) %>%
      arrange(desc(len)) %>%
      group_by(type) %>%
      slice(1) %>%
      ungroup()
    
    # Does genotype or environment have the longest vector?
    longest_vector_type <- subset(longest_vectors_df, len == max(len), type, drop = TRUE)
    
    # Calculate a scalar
    longest_vectors <- sort(longest_vectors_df$len, decreasing = TRUE)
    longest_vectors[2] <- 1 / longest_vectors[2]
    d_scalar <- sqrt(prod(longest_vectors))
    
    # Adjust the loadings
    .x1 <- mutate_at(.x, vars(contains("fa")), ~ ifelse(type == longest_vector_type, . / d_scalar, . * d_scalar))
    
    
    # Find a convex hull around the genotype points
    # These are the extreme genotypes
    extreme_genotypes <- filter(.x1, type == "genotype") %>% 
      slice(chull(select(., fa1, fa2)))
    
    ## Create mega-environments ##
    
    # First identify lines perpedicular with the polygon edges
    me_delineation <- extreme_genotypes %>%
      bind_rows(., slice(., 1)) %>% # Add the first line to complete the polygon
      mutate(rise = c(diff(fa2), 0), run = c(diff(fa1), 0),
             slope = rise / run, intercept = fa2 - (slope * fa1), slope_perp = -(1 / slope)) %>% # Calculate slope of perpendicular line 
      mutate(fa1_perp = intercept / (slope_perp - slope), fa1_temp = c(fa1[-1], fa1[1])) %>%
      # Determine if the perpendicular line intersects the segment
      mutate(left = pmin(fa1, fa1_temp), right = pmax(fa1, fa1_temp),
             intersect = fa1_perp>= left & fa1_perp <= right) %>%
      # Filter out the non-intersecting lines
      filter(intersect) %>%
      # Create coordinates for the segments
      mutate(fa2_end = fa1_perp * slope_perp, fa1_end = fa1_perp, fa1 = 0, fa2 = 0, 
             mega_environment = LETTERS[seq(nrow(.))]) %>%
      # Select relevant columns
      select(mega_environment, term, slope = slope_perp, fa1, fa2, fa1_end, fa2_end)
    
    
    
    ## Separate environments into mega-environments
    
    # Find the beginning and end angles of each ME
    me_delineation1 <- me_delineation %>%
      mutate(angle = map2_dbl(fa1_end, fa2_end, ~angle(v1 = c(.x, .y), v2 = c(0, 1000))),
             angle = ifelse(sign(fa1_end) == -1, 360 - angle, angle)) %>%
      arrange(angle) %>%
      bind_rows(., slice(., 1)) %>%
      mutate(angle_end = c(angle[-1], last(angle))) %>%
      head(-1) %>%
      select(mega_environment, angle, angle_end) %>%
      mutate(check = angle_end < angle)
    
    # Add additional range to capture angles close to 360
    if (any(me_delineation1$check)) {
      me_delineation1 <- add_row(me_delineation1, angle = 0, angle_end = subset(me_delineation1, check, angle_end, drop = TRUE),
                                 mega_environment = subset(me_delineation1, check, mega_environment, drop = TRUE),
                                 check = FALSE) %>%
        mutate(angle_end = ifelse(check, 360, angle_end))
      
    }
    
    

    # Find the angle of each environment vector
    environment_me_assign <- filter(.x1, type == "environment") %>%
      mutate(angle = map2_dbl(fa1, fa2, ~angle(v1 = c(.x, .y), v2 = c(0, 1000))),
             angle = ifelse(sign(fa1) == -1, 360 - angle, angle),
             angle = ifelse(angle < min(me_delineation1$angle), angle + 360, angle)) %>%
      # For each angle, determine the megaenvironment
      mutate(mega_environment = map_chr(angle, ~subset(me_delineation1, .x >= me_delineation1$angle & .x <= me_delineation1$angle_end, 
                                                       mega_environment, drop = TRUE) ))
    
    # Calulate the number of environments per ME
    me_delineation2 <- me_delineation1 %>%
      left_join(., summarize(group_by(environment_me_assign, mega_environment), nEnv = n()), by = "mega_environment") %>%
      mutate(nEnv = ifelse(is.na(nEnv), 0, nEnv),
             nEnvSort = ifelse(nEnv == 0, 0, 1),
             mega_environment1 = fct_reorder2(mega_environment, angle, nEnvSort, .desc = TRUE),
             mega_environment1 = paste0("ME", as.numeric(mega_environment1)),
             me_annotation = paste0(mega_environment1, " (n = ", nEnv, ")"),
             me_annotation = as.factor(me_annotation))
    
    
    ## Data to plot
    dat_toplot <- left_join(.x1, select(environment_me_assign, term, mega_environment), by = "term") %>%
      left_join(., me_delineation2, by = "mega_environment")
    
    ## Plot
    g_me <- ggplot(data = dat_toplot, aes(x = fa1, y = fa2)) +
      # geom_hline(yintercept = 0, lwd = 0.5) +
      # geom_vline(xintercept = 0, lwd = 0.5) +
      geom_point(data = filter(dat_toplot, type == "genotype"), size = 0.5, color = "grey85", shape = 3) +
      geom_point(data = filter(dat_toplot, type == "environment"), aes(color = me_annotation), size = 1) +
      geom_polygon(data = extreme_genotypes, fill = alpha("white", 0), color = "grey85") +
      geom_text(data = extreme_genotypes, aes(label = term), size = 2.5, hjust = "inward", vjust = "inward") +
      geom_segment(data = me_delineation, aes(x = fa1, y = fa2, xend = fa1_end, yend = fa2_end), 
                   lty = 2, inherit.aes = FALSE) +
      scale_color_discrete(name = NULL) +
      scale_x_continuous(name = "Loading 1", breaks = pretty) +
      scale_y_continuous(name = "Loading 2", breaks = pretty) +
      labs(subtitle = paste0(c(str_add_space(tr), toupper(nur), str_to_title(mgmt)), collapse = ", ")) +
      theme_acs(base_size = 8) +
      theme(legend.justification = c("left", "top"), legend.position = c(0.001, 0.999),
            legend.key.height = unit(0.5, "line"), legend.background = element_rect(fill = alpha("white", 0)))
    
    ## Edit the plot
    if (nlevels(droplevels(dat_toplot$me_annotation)) > 4) {
      suppressMessages(g_me <- g_me + scale_color_discrete(name = NULL, guide = guide_legend(ncol = 2)))
      
    }
    
    # Return the data and the plot
    tibble(me_data = list(dat_toplot), plot = list(g_me))
    
  }) %>% ungroup()





## Calculate average ME loadings and average location loadings
environmental_loadings <- nursery_mega_environment_output %>% 
  unnest(me_data) %>% 
  filter(type == "environment") %>% 
  rename(environment = term) %>% 
  left_join(distinct(trial_metadata, environment, location)) %>% 
  select(trait, nursery, management, environment, mega_environment = mega_environment1, me_annotation, location, fa1, fa2)

average_me_loadings <- environmental_loadings %>%
  group_by(nursery, management, trait, mega_environment, me_annotation) %>%
  summarize_at(vars(contains("fa")), mean) %>%
  ungroup()

average_location_loadings <- environmental_loadings %>%
  inner_join(., candidate_locations)  %>%
  group_by(nursery, management, trait, location) %>%
  summarize(fa1 = mean(fa1), fa2 = mean(fa2), n = n()) %>%
  ungroup()


## For each location, calculate the angle between the average location axis and
## the ME axis; then calculate the angle between each composite environment
## and the location axis. This will measure representativeness of the location
## and the stability of that representativeness

# First combine data
nursery_mega_environment_output4 <- environmental_loadings %>% 
  inner_join(., rename_at(average_me_loadings, vars(contains("fa")), ~paste0("me_", .))) %>% 
  inner_join(., rename_at(average_location_loadings, vars(contains("fa")), ~paste0("loc_", .)))

# Calculate angles
nursery_mega_environment_angles <- nursery_mega_environment_output4 %>%
  mutate(environment_loc_angles = pmap_dbl(list(fa1, fa2, loc_fa1, loc_fa2), ~angle(v1 = c(..1, ..2), v2 = c(..3, ..4)))) %>%
  # NA angles are 0
  mutate_at(vars(contains("angles")), ~ifelse(is.na(.), 0, .))

# For each location, calculate the mean absolute angle (maa) about the mean
nursery_location_angles_variance <- nursery_mega_environment_angles %>% 
  group_by(nursery, management, trait, location) %>%
  summarize_at(vars(environment_loc_angles), list(maa = mean)) %>%
  ungroup()

# For each ME, calculate the angle between the average location axis and the ME axis
nursery_mega_environment_location_angles <- nursery_mega_environment_angles %>%
  select(nursery, management, trait, mega_environment, me_fa1, me_fa2, location, loc_fa1, loc_fa2) %>%
  distinct() %>% 
  {left_join(x = distinct(select(., -location, -loc_fa1, -loc_fa2)), 
             y = distinct(select(., -mega_environment, -me_fa1, -me_fa2)))} %>%
  mutate(location_me_angles = pmap_dbl(list(me_fa1, me_fa2, loc_fa1, loc_fa2), ~angle(v1 = c(..1, ..2), v2 = c(..3, ..4))))
  

## Produce a final summary that includes the angle between each location and each ME and the variance
## of environmental angles about the location mean
nursery_environment_angles_summary <- nursery_mega_environment_location_angles %>%
  left_join(., nursery_location_angles_variance)

# For each location, find the closest mega-environment
location_me_representativeness <- nursery_environment_angles_summary %>%
  select(nursery, management, trait, mega_environment, location, location_me_angles, maa) %>%
  # Assign locations to MEs
  arrange(nursery, management, trait, location, location_me_angles) %>%
  group_by(nursery, management, trait, location) %>%
  slice(1) %>%
  group_by(nursery, management, trait) %>%
  nest(.key = "representativeness") %>%
  ungroup() %>%
  mutate(representativeness = map(representativeness, ~arrange(.x, mega_environment, location_me_angles, maa)),
         representativeness = map(representativeness, ~mutate(.x, mega_environment = paste0("ME", as.numeric(as.factor(mega_environment))))))




# Summarize reliability to determine discriminatory locations -----------------


# Reliability will be measured using the estimated residual variance,




# ## Example differences in LSD based on within-trial variance
# var_comp_tidy %>% 
#   inner_join(., candidate_locations) %>%
#   filter(type == "varR") %>%
#   left_join(., summarize(group_by(pheno_dat, trait, trial), n = n())) %>%
#   mutate(LSD = qt(p = 0.05 / 2, df = n - 1, lower.tail = FALSE) * (sqrt(2 * component))) %>%
#   group_by(trait, location) %>%
#   # group_by(nursery, management, trait, location) %>%
#   summarize_at(vars(component, LSD), mean) %>%
#   summarize_at(vars(component, LSD), list(~min, ~max, ~mean)) %>%
#   as.data.frame()



 

 
## Calculate the average reliabilty per location
location_avg_reliability <- var_comp_tidy %>%
  inner_join(., candidate_locations) %>%
  filter(type == "varR") %>% 
  group_by(trait, nursery, management, location) %>%
  summarize(varR_mean = mean(component)) %>%
  nest(.key = "reliability") %>%
  ungroup()



# Calculate correlations of locations over years ---------------------

env_cor_df <- env_cor_matrix %>%
  mutate(cor_df = map(cor_mat, ~as.data.frame(.x) %>% rownames_to_column("environment.x") %>%
                        gather(environment.y, correlation, -environment.x))) %>%
  select(-cor_mat)


# Tidy the correlations within location and between locations
location_correlations <- env_cor_df %>%
  mutate(cor_df = map(cor_df, ~filter(., environment.x != environment.y))) %>%
  unnest(cor_df) %>%
  # Add locations
  left_join(., distinct(trial_metadata, environment, location), by = c("environment.x" = "environment")) %>%
  left_join(., distinct(trial_metadata, environment, location), by = c("environment.y" = "environment")) %>%
  inner_join(., candidate_locations, by = c("trait", "nursery", "management", "location.x" = "location")) %>%
  inner_join(., candidate_locations, by = c("trait", "nursery", "management", "location.y" = "location")) %>%
  group_by(nursery, management, trait) %>% 
  nest(.key = "correlations") %>%
  ungroup()
  



location_repeatability <- location_correlations %>%
  rename(repeatability = correlations)


## Alternative ME representativeness - calculate the average correlation between environments
## in a location with those in a mega-environment
## 
# 
  
# Df of environments in mega-environments
envs_in_mes <- select(environmental_loadings, nursery, trait, environment, mega_environment) %>%
  group_by(nursery, trait, mega_environment) %>%
  nest(envs_in_me = c(environment)) %>%
  ungroup()

# df of environments in locations
envs_in_locs <- unnest(env_cor_df, cols = cor_df) %>% 
  select(nursery, trait, environment = environment.x) %>% 
  distinct() %>% 
  left_join(., distinct(trial_metadata, nursery, environment, location)) %>%
  group_by(nursery, trait, location) %>%
  nest(envs_in_loc = c(environment)) %>%
  ungroup()

# Cross within traits
location_me_representativeness_correlations <- full_join(envs_in_mes, envs_in_locs) %>%
  # add correlations
  left_join(., env_cor_df) %>%
  # Calculate average correlations between envs in the me and envs in the loc
  mutate(
    me_loc_correlation = pmap(select(., envs_in_me, envs_in_loc, cor_df), 
                              ~filter(..3, environment.x %in% ..1$environment, environment.y %in% ..2$environment) ),
    correlation = map_dbl(me_loc_correlation, ~mean(pull(.x, correlation)))
  ) %>%
  select(nursery, management, trait, mega_environment, location, me_loc_correlation, correlation) %>%
  ## Subset for location-me combinations of assigned locations within mes
  group_by(nursery, management, trait) %>%
  nest(cor_repre = c(location, mega_environment, correlation, me_loc_correlation )) %>%
  ungroup() %>%
  left_join(., location_me_representativeness) %>%
  mutate(cor_repre = map2(cor_repre, representativeness, ~right_join(.x, select(.y, mega_environment, location), by = c("location", "mega_environment"))),
         cor_repre = map(cor_repre, ~arrange(., desc(correlation)))) %>%
  select(-representativeness)
  
  



# Calculate precision (variance of genotype mean) per environment ---------

environmental_precision <- var_comp_matrices %>%
  group_by(trait, nursery, management) %>%
  do({
    
    row <- .
    
    tr <- row$trait
    G <- row$Rho_scaled[[1]]
    R <- `dimnames<-`(diag(as.numeric(row$varR_scaled[[1]])), 
                      replicate(2, row.names(row$varR_scaled[[1]]), simplify = FALSE))
    G_us <- row$Rho[[1]]
    R_us <- `dimnames<-`(diag(as.numeric(row$varR[[1]])), 
                         replicate(2, row.names(row$varR[[1]]), simplify = FALSE))
    
    # Calculate varY
    env_varY <- tibble(environment = row.names(G_us), varY = map_dbl(seq_len(nrow(G_us)), ~varYbar(i = .x, G = G_us, R = R_us)))
    
    
    # Locations to use
    location_candidates_tr <- left_join(select(row, trait, nursery, management), candidate_locations, 
                                        by = c("trait", "nursery", "management"))$location
    
    
    # Create a list of environments per locations per year
    env_varY1 <- pheno_dat %>%
      inner_join(., distinct(inner_join(trial_metadata, row, by = c("nursery", "management")), trait, environment), 
                 by = c("environment", "trait")) %>%
      filter(location %in% location_candidates_tr) %>%
      distinct(location, year, environment) %>%
      droplevels() %>%
      # Add varY
      left_join(., env_varY, by = "environment")
      
    # Return the df
    tibble(varY = list(env_varY1), G = list(G_us), R = list(R_us))
    
  }) %>% ungroup()


# ## Add environmental mean to compare precision with mean
# 
# # Read in the phenotype data used for modeling
# load(file.path(data_dir, "data_for_asreml_modeling.RData"))
# 
# environmental_precision_mean <- pheno_to_model %>% 
#   unnest(data) %>% 
#   group_by(trait, nursery, management, location, year, environment) %>%
#   summarize(mean = mean(value)) %>%
#   ungroup() %>%
#   left_join(unnest(environmental_precision), .)
# 
# # Calculate average precision and mean per location
# location_precision_mean <- environmental_precision_mean %>% 
#   group_by(trait, nursery, location, management) %>%
#   summarize_at(vars(varY, mean), mean) %>%
#   ungroup()
# 
# 
# 
# # Correlate mean and variance per trait/nursery/management
# location_precision_mean %>%
#   group_by(trait, nursery, management) %>%
#   do(cor_test = cor.test(.$varY, .$mean)) %>%
#   ungroup() %>%
#   mutate(cor_varY_mean = map_dbl(cor_test, "estimate"),
#          p_value = map_dbl(cor_test, "p.value")) %>%
#   select(-cor_test) %>%
#   as.data.frame()
# 
# # Plot
# location_precision_mean %>%
#   ggplot(aes(x = mean, y = varY)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ trait, scales = "free")
# 
# # Filter outliers and plot again
# location_precision_mean %>%
#   filter(
#     !(trait == "BetaGlucan" & varY > 2e06),
#     !(trait == "DiastaticPower" & varY > 7500),
#     !(trait == "HeadingDate" & varY > 50),
#     !(trait == "PlumpGrain" & varY > 2),
#     !(trait == "PlantHeight" & mean > 150)
#   ) %>%
#   ggplot(aes(x = mean, y = varY, color = nursery, group = nursery)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ trait, scales = "free")
# 
# ## Rank 



## No consistent relationship between mean and variance; 
## most concerning traits are DP and GP (wrn)





# Calculate trait heritabilities ------------------------------------------

# Calculate the heritability within each environment using the estimates of 
# within-environment genetic and residual variances
environmental_precision %>%
  mutate(varG_bar = map_dbl(G, ~mean(diag(.))),
         varGE = map_dbl(G, calc_varGE),
         varRbar = map_dbl(R, ~mean(diag(.))),
         varYbar = map2_dbl(G, R, ~varYbar(i = seq_len(nrow(.x)), G = .x, R = .y)),
         nE = map_dbl(G, nrow),
         h2 = varG_bar / (varG_bar + varGE + varRbar)) %>%
  select(-varY, -G, -R) %>%
  # summarize
  mutate(varGE2varG = varGE / varG_bar) %>%
  as.data.frame() %>%
  arrange(nursery, varGE2varG)

  
# Combine reliability, repeatability, and representativeness -------------------


# Combine data
location_performance <- list(location_repeatability, location_avg_reliability, location_me_representativeness, location_me_representativeness_correlations) %>%
  reduce(full_join)

## Create a df for plotting
location_performance_toplot <- location_performance %>%
  mutate(repeatability = map(repeatability, ~ filter(.x, location.x == location.y) %>% 
                               mutate(delim = map2_chr(environment.x, environment.y, ~paste0(sort(c(.x, .y)), collapse = ":"))) %>%
                               filter(!duplicated(delim)) %>% select(-delim) %>% rename(location = location.x) %>%
                               group_by(location) %>% summarize(repeatability = mean(correlation), .groups = "drop")),
         # Combine data
         performance = pmap(list(repeatability, reliability, representativeness, cor_repre), ~reduce(list(..1, ..2, ..3, ..4), left_join, by = "location"))) %>%
  select(-repeatability:-cor_repre) %>%
  unnest(performance) %>%
  rename(precision = varR_mean, representativeness1 = location_me_angles, representativeness2 = correlation)




# Combine information for choosing locations ------------------------------

## Calculate genotype-location interaction matrices
genotype_location_matrices <- var_comp_matrices %>%
  group_by(trait, nursery, management) %>%
  do({
    
    row <- .
    
    tr <- row$trait
    G <- row$Rho_scaled[[1]]
    R <- `dimnames<-`(diag(as.numeric(row$varR_scaled[[1]])), 
                      replicate(2, row.names(row$varR_scaled[[1]]), simplify = FALSE))
    G_us <- row$Rho[[1]]
    R_us <- `dimnames<-`(diag(as.numeric(row$varR[[1]])), 
                         replicate(2, row.names(row$varR[[1]]), simplify = FALSE))
    
    # Locations to use
    location_candidates_tr <- left_join(select(row, trait, nursery, management), candidate_locations, 
                                        by = c("trait", "nursery", "management"))$location
    
    
    # Create a list of environments per locations per year
    pheno_dat1 <- pheno_dat %>%
      inner_join(., distinct(inner_join(trial_metadata, row, by = c("nursery", "management")), trait, environment), 
                 by = c("environment", "trait")) %>%
      filter(location %in% location_candidates_tr) %>%
      distinct(location, year, environment) %>%
      droplevels()
    
    env_list <- pheno_dat1 %>%
      split(.$year) %>%
      map(droplevels) %>%
      map(~split(., .$location) %>% map_chr("environment"))
    
    # Transpose to get environments per location
    env_by_loc <- pheno_dat1 %>%
      split(.$location) %>%
      map("environment")
    
    
    ####
    
    ## Calculate GL and RL based on unscaled covariance matrices
    # For each year, extract the variance-covariance of environments
    env_vcov_list <- map(env_list, ~G_us[row.names(G_us) %in% .x, row.names(G_us) %in% .x, drop = FALSE])
    
    # Calculate the average covariance between pairs of locations
    loc_vcov <- crossing(location1 = sort(unique(unlist(map(env_list, names)))), location2 = location1) %>%
      mutate(covariance = map2_dbl(location1, location2, ~{
        mean(unlist(sapply(env_vcov_list, function(V) V[ row.names(V) %in% env_by_loc[[.x]], colnames(V) %in% env_by_loc[[.y]] ])))
      }))
    
    # Convert to matrix
    GL_us <- loc_vcov %>% 
      # Convert NA to 0
      mutate(covariance = ifelse(is.na(covariance), 0, covariance)) %>%
      spread(location2, covariance) %>% 
      as.data.frame() %>% 
      column_to_rownames("location1") %>%
      as.matrix()
    
    # Calculate the average residual variance per location
    RL_us <- map(env_by_loc, ~R_us[row.names(R_us) %in% .x, row.names(R_us) %in% .x, drop = FALSE]) %>%
      subset(., !sapply(., is_empty)) %>%
      map(~as.matrix(mean(diag(.x)))) %>%
      imap_dbl(~`dimnames<-`(.x, list(.y, .y))) %>%
      {`dimnames<-`(diag(x = .), list(names(.), names(.))) }
    
    
    # Return the matrices
    tibble(GL = list(GL_us), RL = list(RL_us))
    
  }) %>% ungroup()


# Summarize location repeatability - convert to a matrix
location_repeatability1 <- location_repeatability %>%
  mutate(repeatability = map(repeatability, ~ filter(.x, location.x == location.y) %>% 
                               mutate(delim = map2_chr(environment.x, environment.y, ~paste0(sort(c(.x, .y)), collapse = ":"))) %>%
                               filter(!duplicated(delim)) %>% select(-delim) %>% rename(location = location.x) %>%
                               group_by(location) %>% summarize(repeatability = mean(correlation), .groups = "drop"))) %>%
  mutate(repeatability = map(repeatability, ~as.data.frame(.x) %>% column_to_rownames("location") %>% as.matrix())) %>%
  ungroup()

# Representativeness
location_me_representativeness1 <- location_performance_toplot %>%
  group_by(nursery, management, trait) %>%
  nest() %>%
  ungroup() %>%
  mutate(data = map(data, ~select(.x, location, contains("repre")) %>% as.data.frame() %>% 
                      column_to_rownames("location") %>% as.matrix() ),
         representativeness1 = map(data, ~.[,"representativeness1", drop = FALSE]),
         representativeness2 = map(data, ~.[,"representativeness2", drop = FALSE]) ) %>%
  select(-data)
  


# Combine data.frames
location_opimization_input <- reduce(list(genotype_location_matrices, location_repeatability1, 
                                          location_me_representativeness1), full_join) %>%
  mutate(locations = map(GL, ~sort(row.names(.)))) %>%
  mutate_at(vars(GL, RL), ~map2(.x = ., .y = locations, ~.x[.y, .y])) %>%
  # mutate_at(vars(repeatability, representativeness), ~map2(.x = ., .y = locations, ~mutate(.x, location = factor(location, levels = .y)) %>%
  #                                                            arrange(location) ))
  mutate_at(vars(repeatability, contains("repres")), ~map2(.x = ., .y = locations, ~.x[.y,,drop = FALSE])) %>%
  # Calculate a scaling factor for varY based on the average varY when using a single location 
  # (this should presumably be the highest varY)
  mutate(varY = map2(GL, RL, ~matrix(data = sapply(X = seq_along(diag(.x)), FUN = varYbar, G = .x, R = .y), 
                                     dimnames = list(row.names(.x), "varY"))),
         varY_scaling = map_dbl(varY, mean))





# Rank the best locations for each trait ----------------------------------


# Rank each location for each metric within each trait
# Sum the ranks to get an overall score within each trait
# Find the locations with the lowest score per trait
# Summarize across traits
#

# Reorganize the data; calculate varY for each location
location_opimization_input_torank <- location_opimization_input %>%
  mutate_at(vars(varY, repeatability, contains("repr")), 
            ~map(., ~as.data.frame(.) %>% rownames_to_column(., "location"))) %>%
  mutate(data = pmap(select(., varY, repeatability, contains("repr")), 
                     ~reduce(list(..1, ..2, ..3, ..4), full_join, by = "location"))) %>%
  select(trait, nursery, management, data) %>%
  unnest(data) %>%
  # Add ME information
  left_join(., select(unnest(location_me_representativeness, representativeness), nursery:location))

# Calculate rank
location_performance_rank <- location_opimization_input_torank %>%
  group_by(nursery, management, trait, mega_environment) %>%
  mutate_at(vars(representativeness1, varY), list(rank = rank)) %>%
  mutate_at(vars(repeatability, representativeness2), list(rank = ~rank(-.))) %>%
  ungroup()

# Sum the ranks
location_performance_score <- location_performance_rank %>%
  select(-representativeness1_rank) %>%
  mutate(score = rowSums(select(., contains("rank")))) %>%
  group_by(nursery, management, trait, mega_environment) %>%
  mutate(score_rank = rank(score, ties.method = "min")) %>%
  select(nursery, management, trait, mega_environment, location, score, score_rank) %>%
  ungroup()


##

# Get the location with the lowest score for each trait
location_performance_score1 <- location_performance_score %>%
  group_by(nursery, management, trait, mega_environment) %>%
  top_n(x = ., n = 1, wt = -score) %>%
  ungroup()

# For each nursery, list the subset of optimal locations
# Count traits per best location
nursery_best_rank_locations <- location_performance_score1 %>%
  group_by(nursery, management, location) %>%
  nest(traits = c(trait)) %>%
  mutate(traits = map(traits, 1))

# Combine the ranks 
location_performance_score1 <- location_performance_score %>%
  full_join(., nursery_best_rank_locations) %>%
  mutate(nTraits = map_dbl(traits, length)) %>%
  arrange(nursery, management, trait, score, desc(nTraits)) %>%
  group_by(nursery, management, trait, mega_environment) %>%
  slice(1) %>%
  ungroup()

# For each nursery, list the subset of optimal locations
# Count traits per best location
nursery_best_rank_locations <- location_performance_score1 %>%
  select(-traits) %>%
  group_by(nursery, management, location) %>%
  nest(traits = c(trait)) %>%
  mutate(traits = map(traits, 1))





# Resample environments within locations and calculate precision -----------------------

# Number of reps of sampling
nRep <- 10

# Resample by trait, nursery, management
varY_resampling <- environmental_precision %>%
  group_by(trait, nursery, management) %>%
  do(out = {
    row <- .
    
    G <- row$G[[1]]
    R <- row$R[[1]]
    envnames <- row.names(G)
    
    # Get the environments within locations
    envs_within_locs <- row$varY[[1]] %>%
      distinct(location, environment) %>%
      split(.$location)
    
    # Cross locations by number of environments to resample
    location_resampling_plan <- envs_within_locs %>%
      imap_dfr(~{
        df <- .x
        crossing(location = .y, nEnvSample = seq(2, nrow(df) - 1), rep = seq_len(nRep)) %>%
          mutate(environments = map(nEnvSample, ~sample(df$environment, size = .)))
      })
    
    # Calculate varY for each resample
    varY_location_resamples <- location_resampling_plan %>%
      mutate(varY = map_dbl(environments, ~varYbar(i = which(envnames %in% .x), G = G, R = R)))
    
    select(varY_location_resamples, -environments)
    
  }) %>% ungroup()









# Save data for figure construction ---------------------------------------

save_list <- c("environmental_precision", "loadings_prop_var", "location_correlations", "location_performance_toplot", 
               "candidate_locations", "nursery_mega_environment_output", "pedigree_model_comparison", 
               "location_opimization_input", "varY_resampling", "nursery_best_rank_locations")

save(list = save_list, file = file.path(result_dir, "location_analysis_results.RData"))


