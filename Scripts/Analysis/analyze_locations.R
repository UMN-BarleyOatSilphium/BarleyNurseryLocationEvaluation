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
  select(trait, nursery, management, trial, environment, location, year, type, component, std.error) %>%
  arrange(nursery, management, trait, location, year) %>%
  # Remove any locations with less than x environments
  group_by(nursery, management, trait, location) %>% filter(n_distinct(environment) >= min_env) %>% ungroup()

## Create matrices of the loadings, residual genetic variance, and residual variance
var_comp_matrices <- var_comp_tidy %>%
  select(trait:management, environment, type, component) %>%
  group_by(trait, nursery, management) %>%
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
  select(-data)


## Calculate the proportion of genetic variance explained by loadings vs residual genetic variance
loadings_prop_var <- var_comp_matrices %>%
  mutate(prop_var = map2(Psi, Rho, ~{
    # Variance explained by residual genetic variance
    var_prop_resid <- mean(diag(.x) / diag(.y))
    # Return tibble
    tibble(loadings_var_prop = 1 - var_prop_resid, resid_var_prop = var_prop_resid)
  })) %>% unnest(prop_var) %>%
  select(trait, nursery, management, contains("var_prop"))
  



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
             mega_environment = paste0("ME", seq(nrow(.)))) %>%
      # Select relevant columns
      select(mega_environment, term, slope = slope_perp, fa1, fa2, fa1_end, fa2_end)
    
    
    
    ## Separate environments into mega-environments
    
    # Find the beginning and end angles of each ME
    me_delineation1 <- me_delineation %>%
      mutate(angle = map2_dbl(fa1_end, fa2_end, ~angle(v1 = c(.x, .y), v2 = c(0, 1000))),
             angle = ifelse(sign(fa1_end) == -1, 360 - angle, angle)) %>%
      bind_rows(., slice(., 1)) %>%
      mutate(angle_end = c(angle[-1], last(angle))) %>%
      head(-1) %>%
      select(mega_environment, angle, angle_end) %>%
      mutate(angle_end = ifelse(angle_end < angle, 360 + angle_end, angle_end))
    
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
             me_annotation = paste0(mega_environment, " (n = ", nEnv, ")"),
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
      scale_color_discrete(name = NULL, drop = FALSE) +
      scale_x_continuous(name = "Loading 1", breaks = pretty) +
      scale_y_continuous(name = "Loading 2", breaks = pretty) +
      labs(subtitle = paste0(c(str_add_space(tr), toupper(nur), str_to_title(mgmt)), collapse = ", ")) +
      theme_acs(base_size = 10) +
      theme(legend.justification = c("left", "top"), legend.position = c(0.001, 0.999),
            legend.key.height = unit(0.5, "line"), legend.background = element_rect(fill = alpha("white", 0)))
    
    ## Edit the plot
    if (nlevels(dat_toplot$me_annotation) > 4) {
      suppressMessages(g_me <- g_me + scale_color_discrete(name = NULL, drop = FALSE, guide = guide_legend(ncol = 2)))
      
    }
    
    # Return the data and the plot
    tibble(me_data = list(dat_toplot), plot = list(g_me))
    
  }) %>% ungroup()


## Save plots
nursery_mega_environment_output_toplot <- nursery_mega_environment_output %>%
  mutate(delim = paste(nursery, management, sep = "_")) %>%
  split(.$delim)

## Iterate and save plots
for (nursery_output in nursery_mega_environment_output_toplot) {
  
  # Create a filename
  filename <- paste0("mega_environments_", unique(nursery_output$nursery), "_",
                     unique(nursery_output$management), ".jpg")
  # Create a plot grid
  g_plot_grid <- plot_grid(plotlist = nursery_output$plot, ncol = 4)
  # Number of rows
  plot_rows <- ceiling(nrow(nursery_output) / 4)
  
  # Save
  ggsave(filename = filename, plot = g_plot_grid, path = fig_dir, width = 16, height = plot_rows * 4,
         dpi = 1000)
  
  
}




## Calculate average ME loadings and average location loadings
environmental_loadings <- nursery_mega_environment_output %>% 
  unnest(me_data) %>% 
  filter(type == "environment") %>% 
  rename(environment = term) %>% 
  left_join(distinct(trial_metadata, environment, location)) %>% 
  select(trait, nursery, management, environment, mega_environment, me_annotation, location, fa1, fa2)

average_me_loadings <- environmental_loadings %>%
  group_by(nursery, management, trait, mega_environment, me_annotation) %>%
  summarize_at(vars(contains("fa")), mean) %>%
  ungroup()

average_location_loadings <- environmental_loadings %>%
  group_by(nursery, management, trait, location) %>%
  summarize(fa1 = mean(fa1), fa2 = mean(fa2), n = n()) %>%
  ungroup()

## Edit each plot to include the average ME axis
nursery_mega_environment_output2 <- nursery_mega_environment_output %>%
  left_join(., nest(group_by(average_me_loadings, nursery, management, trait), .key = "average_me_loadings")) %>%
  left_join(., nest(group_by(average_location_loadings, nursery, management, trait), .key = "average_location_loadings"))


## Plot the average ME axis
nursery_mega_environment_output3 <- nursery_mega_environment_output2 %>%
  mutate(plot1 = pmap(list(plot, average_me_loadings, average_location_loadings), ~{
    g <- ..1
    me_df <- ..2
    loc_df <- ..3
    
    # Create an annotation for locations
    loc_df1 <- mutate(loc_df, annotation = paste0(abbreviate(location), " (", n, ")"))
    
    ## Remove the polygon
    g1 <- g
    g1$layers[3:5] <- NULL
    
    # Add average ME and location axes; return the plot
    g1 + 
      geom_segment(data = me_df, aes(x = 0, y = 0, xend = fa1, yend = fa2, color = me_annotation), lwd = 1) +
      geom_segment(data = loc_df1, aes(x = 0, y = 0, xend = fa1, yend = fa2), lwd = 0.25) +
      geom_label(data = loc_df1, aes(x = fa1, y = fa2, label = annotation), size = 2, hjust = "outward", vjust = "outward",
                 fill = alpha("white", 0.7), label.size = NA)
    
  }))
    

## Save plots
nursery_mega_environment_output_toplot <- nursery_mega_environment_output3 %>%
  mutate(delim = paste(nursery, management, sep = "_")) %>%
  split(.$delim)
   
## Iterate and save plots
for (nursery_output in nursery_mega_environment_output_toplot) {
  
  # Create a filename
  filename <- paste0("mega_environment_location_axes_", unique(nursery_output$nursery), "_",
                     unique(nursery_output$management), ".jpg")
  # Create a plot grid
  g_plot_grid <- plot_grid(plotlist = nursery_output$plot1, ncol = 4)
  # Number of rows
  plot_rows <- ceiling(nrow(nursery_output) / 4)
  
  # Save
  ggsave(filename = filename, plot = g_plot_grid, path = fig_dir, width = 16, height = plot_rows * 4,
         dpi = 1000)
  
  
}


## For each location, calculate the angle between the average location axis and
## the ME axis; then calculate the angle between each composite environment
## and the location axis. This will measure representativeness of the location
## and the stability of that representativeness

# First combine data
nursery_mega_environment_output4 <- environmental_loadings %>% 
  left_join(., rename_at(average_me_loadings, vars(contains("fa")), ~paste0("me_", .))) %>% 
  left_join(., rename_at(average_location_loadings, vars(contains("fa")), ~paste0("loc_", .)))

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
  nest(.key = representitiveness) %>%
  mutate(representitiveness = map(representitiveness, ~arrange(.x, mega_environment, location_me_angles, maa)))




# Summarize reliability to determine discriminatory locations ----------


# Reliability will be measured using the estimated residual variance,


## Plot residuals variance by location and trait
varR_plot_list <- var_comp_tidy %>% 
  mutate(location_abbr = abbreviate(location)) %>%
  filter(type == "varR") %>% 
  group_by(trait, nursery, management) %>%
  do(plot = {
    df <- .
    # Refactor
    df %>% 
      mutate_at(vars(location, location_abbr), ~fct_reorder(., component, median)) %>%
      # plot
      ggplot(., aes(x = location_abbr, y = log10(component))) + 
      geom_jitter(width = 0.25, color = "grey85", size = 0.5) +
      geom_boxplot(alpha = 0, outlier.shape = NA) +
      scale_y_continuous(name = expression(log10(italic(sigma[R]^2)))) +
      scale_x_discrete(name = "Location") +
      # facet_grid(trait ~ nursery, scales = "free_x", space = "free_x", switch = "y",
      #            labeller = labeller(nursery = toupper, trait = str_add_space)) +
      facet_grid(trait ~ ., scales = "free_x", space = "free_x", switch = "y",
                 labeller = labeller(nursery = toupper, trait = str_add_space)) +
      labs(subtitle = paste0(c(str_add_space(unique(df$trait)), toupper(unique(df$nursery)), str_to_title(unique(df$management))), 
                             collapse = ", ")) +
      theme_genetics(10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside")
  }) %>% ungroup()

## Save plots
varR_plot_list %>%
  split(.$nursery) %>%
  imap(~{

  ## Modify plots
  plot_list1 <- map(.x$plot, ~. + theme(axis.title.x = element_blank()))
  
  # Combine plots
  g_plot <- add_sub(plot = plot_grid(plotlist = plot_list1, ncol = 4), label = "Location")
  ggsave(filename = paste0("varR_boxplot_by_location_", .y, ".jpg"), plot = g_plot, path = fig_dir,
         height = 3 * ceiling(length(plot_list1) / 4), width = 12, dpi = 1000)
  
})
 


## Calculate the average reliabilty per location
location_avg_reliability <- var_comp_tidy %>%
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
  group_by(nursery, management, trait) %>% nest(.key = "correlations")
  

# Plot
location_cor_plot_list <- location_correlations %>%
  group_by(nursery, management, trait) %>%
  do(plot = {
    row <- .
    
    nur <- row$nursery
    mgmt <- row$management
    tr <- row$trait
    
    df <- row$correlations[[1]] %>%
      # Filter for matching xy locations
      filter(location.x == location.y) %>%
      mutate(delim = map2_chr(environment.x, environment.y, ~paste0(sort(c(.x, .y)), collapse = ":"))) %>%
      filter(!duplicated(delim)) %>%
      select(-delim)
    
    
    # Reorder location factors - create annotation
    df %>%
      mutate(location_abbr = paste0(abbreviate(location.x, 4), "."),
             location_abbr = fct_reorder(location_abbr, .x = correlation, .fun = median, .desc = TRUE)) %>%
      ggplot(aes(x = location_abbr, y = correlation)) +
      geom_boxplot() +
      scale_y_continuous(name = "Within-location correlation", breaks = pretty) +
      scale_x_discrete(name = "Location") +
      labs(subtitle = paste0(c(str_add_space(tr), toupper(nur), str_to_title(mgmt)), collapse = ", ")) +
      theme_genetics(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside")
  }) %>% 
  ungroup()

## Save plots
location_cor_plot_list %>%
  split(.$nursery) %>%
  imap(~{
  
  ## Modify plots
  plot_list1 <- map(.x$plot, ~. + theme(axis.title = element_blank()))
  
  # Combine plots
  g_plot <- add_sub(plot = plot_grid(plotlist = plot_list1, ncol = 4), label = "Location")
  ggsave(filename = paste0("location_correlations_", .y, ".jpg"), plot = g_plot, path = fig_dir,
         height = 3 * ceiling(length(plot_list1) / 4), width = 12, dpi = 1000)
  
})


location_repeatability <- location_correlations %>%
  rename(repeatability = correlations)




# Combine reliability, repeatability, and representitiveness -------------------


# Combine data
location_performance <- reduce(list(location_repeatability, location_avg_reliability, location_me_representativeness), full_join)

## Create a df for plotting
location_performance_toplot <- location_performance %>%
  mutate(repeatability = map(repeatability, ~ filter(.x, location.x == location.y) %>% 
                               mutate(delim = map2_chr(environment.x, environment.y, ~paste0(sort(c(.x, .y)), collapse = ":"))) %>%
                               filter(!duplicated(delim)) %>% select(-delim) %>% rename(location = location.x) %>%
                               group_by(location) %>% summarize(repeatability = mean(correlation))),
         # Combine data
         performance = pmap(list(repeatability, reliability, representitiveness), ~reduce(list(..1, ..2, ..3), left_join, by = "location"))) %>%
  unnest(performance) %>%
  rename(precision = varR_mean, representitiveness = location_me_angles)


# Function to scale from 0 to 1
scale_01 <- function(x, high.favorable = TRUE) {
  if (high.favorable) {
    ( x - min(x) ) / ( max(x) - min(x) )
  } else {
    x1 <- -(x - max(x))
    x1 / max(x1 - min(x1))
  }}

## Plot
location_performance_plots <- location_performance_toplot %>%
  group_by(nursery, management, trait) %>%
  do(plot = {
    df <- .
    df1 <- mutate(df, representitiveness = scale_01(representitiveness, high.favorable = FALSE))
    
    # Color scale breaks
    color_breaks <- seq(0, 1, by = 0.25)
    names(color_breaks) <- c("Low", rep("", 3), "High")
    
    ggplot(data = df1, aes(x = -log10(precision), y = repeatability)) +
      geom_point(aes(color = representitiveness, shape = mega_environment), size = 2.5) +
      geom_text_repel(aes(label = location), size = 3) +
      scale_x_continuous(name = expression("Precision ["*-log10(sigma[R]^2)*"]"), breaks = pretty) +
      scale_y_continuous(name = expression("Repeatability (within-location "~rho[G]*")"), breaks = pretty) +
      scale_shape_discrete(name = "ME") +
      scale_color_gradient(low = "grey75", high = "blue", breaks = color_breaks, name = "ME Rep") +
      labs(subtitle = paste0(c(str_add_space(unique(df$trait)), toupper(unique(df$nursery)), str_to_title(unique(df$management))), 
                             collapse = ", ")) +
      theme_genetics(10)
    
  }) %>% ungroup()


# Split
location_performance_plots1 <- location_performance_plots %>%
  mutate(delim = paste(nursery, management, sep = "_")) %>%
  split(.$delim)

## Iterate and save plots
for (df in location_performance_plots1) {
  
  # Create a filename
  filename <- paste0("location_overall_performance_", unique(df$nursery), "_",
                     unique(df$management), ".jpg")
  # Create a plot grid
  g_plot_grid <- plot_grid(plotlist = df$plot, ncol = 4)
  # Number of rows
  plot_rows <- ceiling(nrow(df) / 4)
  
  # Save
  ggsave(filename = filename, plot = g_plot_grid, path = fig_dir, width = 20, height = plot_rows * 4,
         dpi = 1000)
  
  
}



## Rank the best locations for each trait
# 
# Rank each location for each quality within each trait
# Sum the ranks to get an overall score within each trait
# Find the locations with the lowest score per trait
# Summarize across traits
# 

# Calculate rank
location_performance_rank <- location_performance_toplot %>%
  group_by(nursery, management, trait, mega_environment) %>%
  mutate_at(vars(representitiveness, precision), list(rank = ~rank)) %>%
  mutate(repeatability_rank = rank(-repeatability)) %>%
  ungroup()

# Sum the ranks
location_performance_score <- location_performance_rank %>%
  mutate(score = rowSums(select(., contains("rank")))) %>%
  select(nursery, management, trait, mega_environment, location, score)


# Get the location with the lowest score for each trait
location_performance_score1 <- location_performance_score %>%
  group_by(nursery, management, trait, mega_environment) %>%
  top_n(x = ., n = 1, wt = -score) %>%
  ungroup()


# For each nursery, list the subset of optimal locations
nursery_optimal_locations <- location_performance_score1 %>%
  distinct(nursery, management, location)

## How many of these were included in trials within the last 5 years?
nursery_optimal_locations %>%
  inner_join(., distinct(filter(trial_metadata, year %in% 2014:2019), nursery, management, location)) %>%
  group_by(nursery, management) %>%
  summarize(nLoc = n_distinct(location))


nursery_optimal_locations %>%
  group_by(nursery, management) %>%
  summarize(nLoc = n_distinct(location))

# nursery management  nLoc
# 1 mvn     rainfed       11
# 2 wrn     irrigated     13
# 3 wrn     rainfed        7

# How many originally?
location_performance_toplot %>% 
  group_by(nursery, management) %>% 
  summarize(nLoc = n_distinct(location))

# nursery management  nLoc
# 1 mvn     rainfed       15
# 2 wrn     irrigated     21
# 3 wrn     rainfed       12







