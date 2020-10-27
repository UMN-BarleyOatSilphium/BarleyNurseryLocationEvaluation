## Barley Nursery Analysis
## 
## Analyze MET from barley nursery
## 
## Use output from mixed model
## 
## Author: Jeff Neyhart
## Last modified: 26 Feb. 2020
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
term_colors <- set_names(paletteer_d(package = "wesanderson", palette = "Zissou1")[c(1,5)], "line_name", "environment")

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
  


## Calculate the proportion of genetic variance explained by loadings vs residual genetic variance
loadings_prop_var <- var_comp_matrices %>%
  mutate(prop_var = map2(Psi, Rho, ~{
    # Variance explained by residual genetic variance
    tibble(environment = row.names(.x), resid_var_prop = diag(.x) / diag(.y)) %>%
      mutate(loadings_var_prop = 1 - resid_var_prop)
  })) %>% unnest(prop_var)





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
      theme_acs(base_size = 10) +
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
  summarize_at(vars(environment_loc_angles), list(maa = ~mean)) %>%
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
  nest(.key = representativeness) %>%
  mutate(representativeness = map(representativeness, ~arrange(.x, mega_environment, location_me_angles, maa)),
         representativeness = map(representativeness, ~mutate(.x, mega_environment = paste0("ME", as.numeric(as.factor(mega_environment))))))




# Summarize reliability to determine discriminatory locations -----------------


# Reliability will be measured using the estimated residual variance,




## Example differences in LSD based on within-trial variance
var_comp_tidy %>% 
  inner_join(., candidate_locations) %>%
  filter(type == "varR") %>%
  left_join(., summarize(group_by(pheno_dat, trait, trial), n = n())) %>%
  mutate(LSD = qt(p = 0.05 / 2, df = n - 1, lower.tail = FALSE) * (sqrt(2 * component))) %>%
  group_by(trait, location) %>%
  # group_by(nursery, management, trait, location) %>%
  summarize_at(vars(component, LSD), mean) %>%
  summarize_at(vars(component, LSD), list(~min, ~max, ~mean)) %>%
  as.data.frame()



 

 
## Calculate the average reliabilty per location
location_avg_reliability <- var_comp_tidy %>%
  inner_join(., candidate_locations) %>%
  filter(type == "varR") %>% 
  group_by(trait, nursery, management, location) %>%
  summarize(varR_mean = mean(component)) %>%
  nest(.key = "reliability")



# Calculate correlations of locations over years ---------------------


# Tidy the correlations within location and between locations
location_correlations <- env_cor_matrix %>%
  mutate(cor_df = map(cor_mat, ~as.data.frame(.x) %>% rownames_to_column("environment.x") %>%
                        gather(environment.y, correlation, -environment.x) %>% filter(environment.x != environment.y))) %>%
  unnest(cor_df) %>%
  # Add locations
  left_join(., distinct(trial_metadata, environment, location), by = c("environment.x" = "environment")) %>%
  left_join(., distinct(trial_metadata, environment, location), by = c("environment.y" = "environment")) %>%
  inner_join(., candidate_locations, by = c("trait", "nursery", "management", "location.x" = "location")) %>%
  inner_join(., candidate_locations, by = c("trait", "nursery", "management", "location.y" = "location")) %>%
  group_by(nursery, management, trait) %>% nest(.key = "correlations")
  



location_repeatability <- location_correlations %>%
  rename(repeatability = correlations)




# Calculate precision (variance of genotype mean) per environment ---------

environmental_precision <- var_comp_matrices %>%
  group_by(trait, nursery, management) %>%
  do(varY = {
    
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
    env_varY1
    
  }) %>% ungroup()





# Combine reliability, repeatability, and representativeness -------------------


# Combine data
location_performance <- reduce(list(location_repeatability, location_avg_reliability, location_me_representativeness), full_join)

## Create a df for plotting
location_performance_toplot <- location_performance %>%
  mutate(repeatability = map(repeatability, ~ filter(.x, location.x == location.y) %>% 
                               mutate(delim = map2_chr(environment.x, environment.y, ~paste0(sort(c(.x, .y)), collapse = ":"))) %>%
                               filter(!duplicated(delim)) %>% select(-delim) %>% rename(location = location.x) %>%
                               group_by(location) %>% summarize(repeatability = mean(correlation))),
         # Combine data
         performance = pmap(list(repeatability, reliability, representativeness), ~reduce(list(..1, ..2, ..3), left_join, by = "location"))) %>%
  unnest(performance) %>%
  rename(precision = varR_mean, representativeness = location_me_angles)




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
                               group_by(location) %>% summarize(repeatability = mean(correlation)))) %>%
  mutate(repeatability = map(repeatability, ~as.data.frame(.x) %>% column_to_rownames("location") %>% as.matrix())) %>%
  ungroup()

# Representativeness
location_me_representativeness1 <- location_me_representativeness %>%
  mutate(representativeness = map(representativeness, ~select(., location, mega_environment, location_me_angles))) %>%
  mutate(representativeness = map(representativeness, ~select(., location, repre = location_me_angles) %>% as.data.frame(.x) %>% 
                                    column_to_rownames("location") %>% as.matrix())) %>%
  ungroup()

# Combine data.frames
location_opimization_input <- reduce(list(genotype_location_matrices, location_repeatability1, 
                                          location_me_representativeness1), full_join) %>%
  mutate(locations = map(GL, ~sort(row.names(.)))) %>%
  mutate_at(vars(GL, RL), ~map2(.x = ., .y = locations, ~.x[.y, .y])) %>%
  # mutate_at(vars(repeatability, representativeness), ~map2(.x = ., .y = locations, ~mutate(.x, location = factor(location, levels = .y)) %>%
  #                                                            arrange(location) ))
  mutate_at(vars(repeatability, representativeness), ~map2(.x = ., .y = locations, ~.x[.y,,drop = FALSE])) %>%
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
  mutate_at(vars(varY, repeatability, representativeness), ~map(., ~as.data.frame(.) %>% rownames_to_column(., "location"))) %>%
  mutate(data = pmap(select(., varY, repeatability, representativeness), .f = ~reduce(list(..1, ..2, ..3), full_join, by = "location"))) %>%
  unnest(data) %>%
  # Add ME information
  left_join(., select(unnest(location_me_representativeness, representativeness), nursery:location))
  
# Calculate rank
location_performance_rank <- location_opimization_input_torank %>%
  group_by(nursery, management, trait, mega_environment) %>%
  mutate_at(vars(repre, varY), list(rank = ~rank)) %>%
  mutate(repeatability_rank = rank(-repeatability)) %>%
  ungroup()

# Sum the ranks
location_performance_score <- location_performance_rank %>%
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
  nest(trait, .key = "traits") %>%
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
  group_by(nursery, management, location) %>%
  nest(trait, .key = "traits") %>%
  mutate(traits = map(traits, 1))



# Optimize selection of locations -----------------------------------------

# Load optimization package
library(gramEvol)


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


# Run the optimization procedure 
optimized_locations_all_traits <- location_opimization_input %>%
  crossing(loc_penalty = loc_penalties, tibble(trait_weight_group = names(trait_weights), trait_weights)) %>%
  arrange(nursery, management, loc_penalty, trait_weight_group) %>%
  group_by(nursery, management, loc_penalty, trait_weight_group) %>%
  do({
    df <- .
    
    # Extract lists of matrices
    GL <- df$GL
    RL <- df$RL
    corL <- df$repeatability
    repL <- df$representativeness
    # Location penalty
    pen <- df$loc_penalty[[1]]
    # Fitness component weights
    comp_weight <- c(varY = 1, repAvg = 0.75, reprAvg = 0.25)
    varYscaling <- df$varY_scaling
    # List of traits
    tr_list <- df$trait
    # Trait weights
    tr_weights <- df$trait_weights[[1]]
    
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
                                evalFunc = function(x) -fitness_int(x = x, locs = all_locations, traits = tr_list, G.list = GL,
                                                                    R.list = RL, P.list = corL, M.list = repL, varY.scaling = varYscaling,
                                                                    component.weights = comp_weight, trait.weights = tr_weights,
                                                                    non.zero.env.list = non_zero_loc_list, env.penalty = pen))
    
 
    # Return the fitness components
    optim_x <- optim_out$best$genome
    agro_optim_loc <- all_locations[optim_x != 0]
    maltq_optim_locs <- all_locations[optim_x == 2]
    fit_out <- fitness_int(x = optim_x, locs = all_locations, traits = tr_list, G.list = GL, R.list = RL, P.list = corL, M.list = repL, 
                           varY.scaling = varYscaling, return.all = TRUE) %>%
      as.data.frame() %>%
      rownames_to_column("trait")
    # Scale and average the fitness components
    fit_out_scaled_trait <- fit_out %>%
      mutate(varY = 1 - (varY / varYscaling), reprAvg = 1 - (reprAvg / 90))
    # Summarize 
    fit_out_scaled <- summarize_at(fit_out_scaled_trait, vars(-trait), weighted.mean, w = tr_weights, na.rm = TRUE)
    
    # Create a tibble that summarizes the optimization results
    optim_results <- tibble(method = "optimization", fitness = list(fit_out), fitness_scaled = list(fit_out_scaled), 
                            fitness_scaled_trait = list(fit_out_scaled_trait), 
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
                       varY.scaling = varYscaling, return.all = TRUE) %>% as.data.frame() %>% rownames_to_column("trait") )
    
    # Rescale fitness
    random_loc_fitness_scaled_trait <- random_loc_fitness %>%
      map(~mutate(., varY = 1 - (varY / varYscaling), reprAvg = 1 - (reprAvg / 90)))
    # Summarize
    random_loc_fitness_scaled <- random_loc_fitness_scaled_trait %>%
      map_df(., ~summarize_at(.x, vars(-trait), weighted.mean, w = tr_weights, na.rm = TRUE))
                            
    # Save random fitness results
    random_results <- tibble(method = "random", fitness = list(random_loc_fitness), fitness_scaled = list(random_loc_fitness_scaled),
                             fitness_scaled_trait = list(random_loc_fitness_scaled_trait))
    
    
    ## Calculate the fitness of the ranked choice environments ##
    best_rank_locations <- subset(nursery_best_rank_locations, nursery %in% df$nursery, location, drop = TRUE)
    # Convert this to 0, 1, 2
    best_rank_x <- loc_values
    best_rank_x[!names(best_rank_x) %in% best_rank_locations] <- 0
    
    best_rank_fitness <- fitness_int(x = best_rank_x, locs = all_locations, traits = tr_list, G.list = GL, R.list = RL, P.list = corL, 
                                     M.list = repL,  varY.scaling = varYscaling, return.all = TRUE) %>% 
      as.data.frame() %>% 
      rownames_to_column("trait")
    
    # Rescale fitness
    best_rank_fitness_scaled_trait <- best_rank_fitness %>%
      mutate(., varY = 1 - (varY / varYscaling), reprAvg = 1 - (reprAvg / 90)) 
    # Summarize 
    best_rank_fitness_scaled <- summarize_at(best_rank_fitness_scaled_trait, vars(-trait), weighted.mean, w = tr_weights, na.rm = TRUE)
    
    # Save ranked locations fitness results
    best_rank_results <- tibble(method = "best_rank", fitness = list(best_rank_fitness), fitness_scaled = list(best_rank_fitness_scaled),
                                fitness_scaled_trait = list(best_rank_fitness_scaled_trait),
                                optim_loc = list(list(agro_loc = all_locations[best_rank_x != 0], 
                                                      maltq_loc = all_locations[best_rank_x == 2])))
                                
    # Combine and return results
    bind_rows(optim_results, random_results, best_rank_results)
    
  }) %>% ungroup()



  



# Save data for figure construction ---------------------------------------

save("environmental_precision", "loadings_prop_var", "location_correlations", "location_performance_toplot", 
     "candidate_locations", "nursery_mega_environment_output", "pedigree_model_comparison", 
     "optimized_locations_all_traits", "location_opimization_input", "trait_weights", 
     file = file.path(result_dir, "figure_data.RData"))


