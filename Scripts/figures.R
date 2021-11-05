## Barley Nursery Analysis
## 
## Script to produce figures
## 
## Author: Jeff Neyhart
## Last modified: 26 Feb. 2020
## 


# Base script
proj_dir <- getwd()
source(file.path(proj_dir, "startup.R"))


library(wesanderson)
library(ggrepel)
library(cowplot)
library(patchwork)
library(grid)

# Color pallete for environments and line names
term_colors <- set_names(wes_palette("Zissou1")[c(1,5)], "line_name", "environment")

# Traits to highlight in the paper
figure_traits <- c("GrainYield", "GrainProtein")

# Load results from the analyses
load(file.path(result_dir, "location_analysis_results.RData"))
load(file.path(result_dir, "location_optimization_results.RData"))
load(file.path(result_dir, "sensitivity_test_results.RData"))


# Theme to use
theme_acs2 <- function(base_size = 8) {
  
  theme_genetics(base_size = base_size) %+replace%
    theme(panel.border = element_rect(fill = alpha("white", 0)))
  
}

# Function for adding padding to formatted numbers
format_numbers2 <- function(x, signif.digits = 2, width = 4) {
  str_pad(string = format_numbers(x = x, signif.digits = signif.digits), width = width, side = "left")
}


# Function for pretty discrete axis text
# Will skip some levels if there are too many
pretty_labels <- function(x, max.levels = 12) {
  if (length(x) > max.levels) {
    by <- ceiling(length(x) / max.levels)
    x[-seq(1, length(x), by = by)] <- ""
  }
  return(x)
}


# Data.frame of location abbreviations and a key
location_abbr_key <- pheno_dat %>%
  mutate(location_abbr = str_sub(string = environment, start = 1, end = 3)) %>%
  select(location, location_abbr) %>%
  distinct() %>%
  mutate_all(as.character)






# Map of locations --------------------------------------------------------


## Plot locations

# Get the map data for canada
canada <- rnaturalearth::ne_states(country = "canada") %>%
  tidy(x = ., region = "name_en") %>%
  mutate(group = as.numeric(as.factor(group)))

# Download map data for US by county
usa_county <- map_data(map = "county")
# Download state data
usa_state <- map_data(map = "state")

# Adjust the groups in the states
usa_state <- usa_state %>%
  mutate(group = group + max(canada$group))

# Adjust the groups in the counties
usa_county <- usa_county %>%
  mutate(group = group + max(usa_state$group))

# Tidy and combine
north_america <- bind_rows(usa_state, usa_county, canada)

# Edit the trial information
trial_info_toplot <- trial_metadata %>%
  filter(!(nursery == "wrn" & management == "rainfed")) %>%
  group_by(nursery, location, lat, long) %>% 
  summarize(nYear = n_distinct(year))

# Annotate common locations
trial_info_toplot1 <- trial_info_toplot %>% 
  group_by(location, lat, long) %>% 
  summarize(nursery = ifelse(n_distinct(nursery) == 1, nursery, "Both")) %>%
  ungroup() %>%
  mutate(nursery_abbr = ifelse(nursery == "Both", nursery, toupper(nursery)),
         nursery = ifelse(nursery != "Both", f_nursery_expand(nursery), nursery),
         nursery = str_remove_all(nursery, " Nursery")) %>%
  mutate_at(vars(contains("nursery")), ~fct_relevel(., "Both", after = Inf))


# Create a box for x and y axes
xlim <- range(pretty(trial_info_toplot$long))
ylim <- range(pretty(trial_info_toplot$lat))

# Map
g_map <- ggplot(data = north_america, aes(x = long, y = lat)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = canada, aes(group = group), fill = "grey95", color = "grey50", lwd = 0.75) + # Add canada
  geom_polygon(data = usa_state, aes(group = group), fill = "grey95", color = "grey50", lwd = 0.75) +
  geom_point(data = trial_info_toplot1, aes(color = nursery), size = 2.5) +
  scale_color_viridis_d(begin = 0.2, end = 0.8, name = NULL) +
  coord_map(projection = "bonne", lat0 = mean(ylim), xlim = xlim, ylim = ylim) +
  theme_void(base_size = 14) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_rect(fill = "white"),
        legend.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "line"))


# Save the figure
ggsave(filename = "trial_location_map.jpg", plot = g_map, path = fig_dir,
       width = 8, height = 5, dpi = 1000)









# Figure 1: Variance of genotype means and correlation plot --------------------

# Edits to plot
environmental_precision_toplot <- environmental_precision %>%
  unnest(varY) %>%
  filter(!is.na(varY)) %>%
  # Calculate median varY per location
  group_by(nursery, trait, location) %>%
  mutate(varY_median = median(varY)) %>%
  ungroup() %>%
  # Add location abbreviations
  left_join(., location_abbr_key) %>%
  # Arrange by ascending varY within nursery-trait combinations
  arrange(nursery, trait, varY_median) %>%
  mutate(x = fct_inorder(paste(nursery, trait, location_abbr, sep = ".")))



# Plot genotypic variance for the figure traits
g_env_varY_figure_traits <- environmental_precision_toplot %>%
  filter(trait %in% figure_traits) %>%
  ggplot(., aes(x = x, y = log10(varY))) + 
  geom_jitter(width = 0.1, color = "grey85", size = 0.2) +
  geom_boxplot(alpha = 0, outlier.shape = NA, lwd = 0.25) +
  scale_y_continuous(name = expression('Precision ['*log[10](italic(V[bar(Y)]))*']'), breaks = pretty, 
                     labels = function(x) ifelse(x > 1e4, formatC(x, format = "e", digits = 1), x), trans = "reverse") +
  scale_x_discrete(name = NULL, labels = function(x) map_chr(str_split(x, "\\."), 3)) +
  facet_wrap(~ nursery + trait, scales = "free", shrink = TRUE,
             labeller = labeller(nursery = toupper, trait = str_add_space, .multi_line = FALSE)) +
  theme_genetics(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), axis.text.y = element_text(size = 5), strip.placement = "outside")
  

# Similarly format the within-location correlation data
location_cor_tidy <- location_correlations %>%
  unnest(correlations) %>%
  inner_join(., candidate_locations, by = c("nursery", "management", "trait", "location.x" = "location")) %>%
  inner_join(., candidate_locations, by = c("nursery", "management", "trait", "location.y" = "location")) %>%
  filter(location.x == location.y) %>%
  mutate(location = location.x, 
         delim = map2_chr(environment.x, environment.y, ~paste0(sort(c(.x, .y)), collapse = ":"))) %>%
  # Add location abbreviations
  left_join(., location_abbr_key) %>%
  group_by(trait, nursery) %>%
  filter(!duplicated(delim)) %>%
  ungroup() %>%
  select(-delim) %>%
  # Arrange by nursery, trait, median location varR
  group_by(nursery, trait, location) %>%
  mutate(loc_median_corR = median(correlation)) %>%
  ungroup() %>%
  arrange(nursery, trait, desc(loc_median_corR))
  

# Plot within-location correlation for the figure traits
g_correlation_figure_traits <- location_cor_tidy %>%
  filter(trait %in% figure_traits) %>%
  # Merge nursery, trait, and location for sorting
  mutate(x = fct_inorder(paste(nursery, trait, location_abbr, sep = "."))) %>%
  ggplot(., aes(x = x, y = correlation)) + 
  geom_jitter(width = 0.1, color = "grey85", size = 0.2) +
  geom_boxplot(alpha = 0, outlier.shape = NA, lwd = 0.25) +
  scale_y_continuous(name = expression("Within-location"~italic(rho[G])), breaks = pretty) +
  scale_x_discrete(name = NULL, labels = function(x) map_chr(str_split(x, "\\."), 3)) +
  facet_wrap(~ nursery + trait, scales = "free", shrink = TRUE,
             labeller = labeller(nursery = toupper, trait = str_add_space, .multi_line = FALSE)) +
  theme_genetics(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), axis.text.y = element_text(size = 5), strip.placement = "outside")

# # Combine plots
# g_precision_correlation <- plot_grid(g_env_varY_figure_traits, g_correlation_figure_traits, nrow = 1,
#                                      labels = letters, label_size = 8, align = "hv")
# # Save
# ggsave(filename = "figure1_reliability_repeatability_example.jpg", path = fig_dir, plot = g_precision_correlation,
#        height = 8, width = 17.5, units = "cm", dpi = 1000)
# 


# Format the representativeness
representativeness_tidy <- location_performance_toplot %>%
  select(trait, nursery, management, location, me = mega_environment.x, me_loc_correlation) %>%
  unnest(me_loc_correlation) %>%
  left_join(., location_abbr_key) %>%
  # Arrange by nursery, trait, median location varR
  group_by(nursery, trait, location, me) %>%
  mutate(loc_median_corR = median(correlation)) %>%
  ungroup() %>%
  arrange(nursery, trait, desc(loc_median_corR))


## Representativeness
g_representativeness <- representativeness_tidy %>%
  filter(trait %in% figure_traits) %>%
  mutate(x = fct_inorder(paste(nursery, trait, location_abbr, sep = "."))) %>%
  ggplot(., aes(x = x, y = correlation)) + 
  geom_jitter(width = 0.1, color = "grey85", size = 0.2) +
  geom_boxplot(alpha = 0, outlier.shape = NA, lwd = 0.25) +
  scale_y_continuous(name = expression("ME representativeness ("*italic(rho[G])*")"), breaks = pretty) +
  scale_x_discrete(name = NULL, labels = function(x) map_chr(str_split(x, "\\."), 3)) +
  facet_wrap(~ nursery + trait, scales = "free", shrink = TRUE,
             labeller = labeller(nursery = toupper, trait = str_add_space, .multi_line = FALSE)) +
  theme_genetics(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), axis.text.y = element_text(size = 5), strip.placement = "outside")


# Combine plots
g_metrics <- plot_grid(g_env_varY_figure_traits, g_correlation_figure_traits, g_representativeness, 
                       ncol = 1, labels = letters, label_size = 8, align = "hv")
# Save
ggsave(filename = "figure1_metrics_example_alt.jpg", path = fig_dir, plot = g_metrics,
       height = 18, width = 8.5, units = "cm", dpi = 1000)






# ## Range in representativeness per trait
# location_opimization_input_toplot %>% 
#   distinct(nursery, management, trait, location, representativeness) %>% 
#   group_by(trait, nursery, management) %>% 
#   summarize_at(vars(representativeness), list(min = min, max = max, mean = mean))


# Combine plots
g_precision_correlation <- plot_grid(g_env_varY_figure_traits, g_correlation_figure_traits, nrow = 1,
                                     labels = letters, label_size = 8, align = "hv")
# Save
ggsave(filename = "figure1_reliability_repeatability_example.jpg", path = fig_dir, plot = g_precision_correlation,
       height = 8, width = 17.5, units = "cm", dpi = 1000)



# # Rank representativeness per trait/nursery
# location_opimization_input_toplot %>% 
#   distinct(nursery, management, trait, location, representativeness) %>% 
#   split(list(.$trait, .$nursery, .$management)) %>%
#   map_df(~arrange(., desc(representativeness))) %>%
#   as.data.frame()










# Figure 2: Location performance plots ------------------------------------

figure_traits1 <- c("GrainYield", "HeadingDate", "GrainProtein", "DiastaticPower", "MaltExtract")

location_opimization_input_toplot <- location_opimization_input %>%
  # Scale representativeness to the maximum angle within a trait/nursery
  # mutate(representativeness = map(representativeness, ~ 1 - (. / max(.)))) %>%
  # Combine data
  mutate_at(vars(repeatability, representativeness2, varY), ~map(., ~as.data.frame(.) %>% rownames_to_column("location"))) %>%
  mutate(data_toplot = pmap(select(., repeatability, representativeness2, varY), list) %>% map(~reduce(., full_join, by = "location"))) %>%
  select(-repeatability, -representativeness2, -varY) %>%
  unnest(data_toplot) %>%
  rename(representativeness = representativeness2) %>%
  # Add location abbreviations
  left_join(., location_abbr_key)


## Plot
location_performance_plots <- location_opimization_input_toplot %>%
  group_by(nursery, management, trait) %>%
  do(plot = {
    df <- .
    df1 <- df %>%
      # mutate(representativeness1 = scale_01(repre, high.favorable = FALSE),
      #        representativeness = 1 - (repre / 90) ) %>%
      mutate(representativeness_bin = cut(x = representativeness, breaks = 5)) %>%
      # Add me
      left_join(., distinct_at(location_performance_toplot, vars(-precision, -repeatability, -contains("repre"), -maa)),
                by = c("trait", "nursery", "management", "location")) %>%
      select(-mega_environment.y) %>% rename(mega_environment = mega_environment.x)
    
    # Color scale breaks
    # color_breaks <- pretty(df1$representativeness, min.n = 5, n = 5)
    # color_breaks <- seq_len(nlevels(df1$representativeness_bin))
    # names(color_breaks) <- c("", "Low", rep("", length(color_breaks) - 4),  "High", "")
    # names(color_breaks) <- c("Low", rep("", length(color_breaks) - 2), "High")
    
    # Use common size breaks
    size_breaks <- seq_len(nlevels(droplevels(df1$representativeness_bin))) / 2
    names(size_breaks) <- c("Low", rep("", length(size_breaks) - 2),  "High")
    
    # Determine the shapes to use (closed for MVN; open for WRN)
    shapes <- if (unique(df1$nursery) == "wrn") c(1, 2, 0, 5, 3, 4) else c(16, 17, 15, 18, 3, 4)
    
    ggplot(data = df1, aes(x = log10(varY), y = repeatability)) +
      geom_text_repel(aes(label = location_abbr), size = 3) +
      geom_point(aes(size = representativeness_bin, shape = mega_environment)) +
      scale_x_reverse(name = expression('Precision ['*log[10](italic(V[bar(Y)]))*']'), breaks = pretty,
                      labels = function(x) ifelse(x > 1e4, formatC(x, format = "e", digits = 1), x)) +
      scale_y_continuous(name = expression("Repeatability (within-location "*italic(rho[G])*")"), breaks = pretty,
                         expand = expansion(mult = 0.1)) +
      scale_shape_manual(name = "Mega environment (ME)", guide = guide_legend(order = 1), values = shapes) +
      scale_size_manual(values = size_breaks, labels = names(size_breaks),  
                        name = expression("ME representativeness ("*italic(rho[G])*")")) +
      labs(subtitle = paste0(c(str_add_space(unique(df1$trait)), toupper(unique(df1$nursery)), 
                               str_to_title(unique(df1$management))), collapse = ", ")) +
      theme_acs2(base_size = 8)
    
  }) %>% ungroup()



# Subset relevant traits
performance_toplot <- filter(location_performance_plots, trait %in% figure_traits1, 
                             !(nursery == "wrn" & management == "rainfed"))

# Generate a common size break
size_breaks <- seq_len(5) / 2
names(size_breaks) <- c("Low", rep("", length(size_breaks) - 2),  "High")

  

# Create merged plots by nursery
merged_plots_by_nursery <- performance_toplot %>%
  mutate(n = as.numeric(as.factor(nursery))) %>%
  group_by(nursery) %>%
  do({
    df <- .
    
    # Determine the shapes to use (closed for MVN; open for WRN)
    shapes <- if (unique(df$nursery) == "wrn") c(1, 2, 0, 5) else c(16, 17, 15, 18)
    
    
    # Plot gy
    performance_partA <- subset(df, trait == "GrainYield", plot, drop = TRUE)[[1]] +
      labs(subtitle = NULL) +
      facet_grid(nursery ~ trait, labeller = labeller(nursery = f_nursery_expand, trait = str_add_space), 
                 switch = "y") +
      scale_shape_manual(name = "Mega environment (ME)", labels = NULL, values = shapes, 
                         guide = guide_legend(order = 1, keywidth = unit(0.5, "pt"))) +
      guides(size = guide_legend(label.position = "top")) +
      theme(strip.placement = "outside", panel.border = element_rect(fill = alpha("white", 0)), 
            plot.title = element_text(hjust = 0.5),
            legend.position  = "bottom", legend.background = element_rect(fill = alpha("white", 0), linetype = 0), 
            legend.box.background = element_rect(fill = alpha("white", 0), linetype = 0), legend.box.just = "bottom")
    
    # # Adjust point and text size
    performance_partA$layers[[1]]$aes_params$size <- 2 # Text
    # performance_partA$layers[[2]]$aes_params$size <- 2
    
    # Plot all other traits
    partB_list <- subset(df, trait != "GrainYield", plot, drop = TRUE)
    performance_partB_plot1 <- partB_list[[1]]
    performance_partB_plot1$data <- bind_rows(map(partB_list, "data")) %>%
      mutate(group = paste0(nursery, "_", trait),
             representativeness_bin = cut(x = representativeness, breaks = 5))
    
    # Use common size breaks
    size_breaks <- seq_len(nlevels(droplevels(performance_partB_plot1$data$representativeness_bin))) / 2
    names(size_breaks) <- c("Low", rep("", length(size_breaks) - 2),  "High")
      
    
    # Adjust text and point size
    performance_partB_plot1$layers[[1]]$aes_params$size <- 1.5 # Text
    # performance_partB_plot1$layers[[2]]$aes_params$size <- 1.5
    
    performance_partB <- performance_partB_plot1 +
      labs(subtitle = NULL) +
      facet_wrap(~ group, ncol = 2, scales = "free", 
                 labeller = labeller(group = function(x) str_add_space(str_remove(x, "[a-z]*_")))) +
      scale_size_manual(values = size_breaks, labels = names(size_breaks), name = "ME representativeness") +
      theme(strip.placement = "outside",
            panel.border = element_rect(fill = alpha("white", 0)), legend.position = "none",
            axis.title = element_blank(), axis.text = element_text(size = 5))
    
    ## Combine plot A and plot B
    if (unique(df$n) == 1) {
      main_plot <- plot_grid(performance_partA + theme(legend.position = "none", axis.title = element_blank()), 
                             performance_partB, nrow = 1, align = "hv", axis = "tbrl",
                             labels = letters[1:2], label_size = 9)
    } else {
      main_plot <- plot_grid(performance_partA + theme(legend.position = "none", axis.title = element_blank()), 
                             performance_partB, nrow = 1, align = "hv", axis = "tbrl")
    }
    
    # Return the main plot and partA to extract titles and the legend
    tibble(main = list(main_plot), partA = list(performance_partA))
    
  }) %>% ungroup()

## Merge
g_performance_merge <- plot_grid(plotlist = merged_plots_by_nursery$main, ncol = 1, 
                                 align = "h", axis = "tb", label_size = 9)
## Add axis titles and the legend
y_axis <- get_plot_component(merged_plots_by_nursery$partA[[1]], "ylab-l")
x_axis <- get_plot_component(merged_plots_by_nursery$partA[[1]], "xlab-b")
legend <- get_legend(merged_plots_by_nursery$partA[[1]])

g_performance_merge1 <- plot_grid(y_axis, g_performance_merge, nrow = 1, rel_widths = c(0.03, 1))
g_performance_merge2 <- plot_grid(g_performance_merge1, x_axis, legend, ncol = 1, rel_heights = c(1, 0.05, 0.09))


# Save
ggsave(filename = "figure2_location_performance.jpg", path = fig_dir, plot = g_performance_merge2,
       height = 12, width = 17.5, units = "cm", dpi = 1000)





# Figure 3: Visualize optimization results --------------------------------

## Part A - composite fitness


# Plot fitness for optimized or random methods of location selection
# Use the equal-weight trait results
optimized_fitness_scaled <- optimized_locations_all_traits %>%
  # Count the number of locations
  mutate(nOptimLoc = map(optim_loc, "agro_loc") %>% map_dbl(length)) %>%
  select(-fitness, -fitness_scaled, -fitness_summarize) %>%
  unnest(fitness_scaled_summarize) %>%
  # Convert varY such that higher values are favorable
  mutate(varY = 1 - varY) %>%
  gather(fitness_component, value, varY, repAvg, reprAvg) %>%
  arrange(nursery, management, loc_penalty, trait_weight_group, comp_weight_group, scale_components, 
          method, fitness_component) %>% # Order
  mutate(loc_penalty = fct_inseq(as.character(loc_penalty)),
         fitness_component = fct_inorder(fitness_component) %>% fct_relevel("varY"),
         method = str_replace_all(method, "_", " ") %>% str_to_title() %>% str_remove_all(" "),
         method = ifelse(method == "BestRank", "Best ranked locations", method),
         method = fct_inorder(method) %>% fct_relevel("Optimization"))

# Subset to plot
optimized_fitness_scaled1 <- optimized_fitness_scaled %>%
  filter(trait_weight_group == "equal", comp_weight_group == "equal", scale_components) %>%
  # Merge method with number of optimal locations
  mutate(nLoc = ifelse(method == "Optimization", nOptimLoc, NA)) %>%
  group_by(nursery, management, loc_penalty, trait_weight_group, comp_weight_group) %>%
  mutate(nLoc = mean(nLoc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(loc_penalty1 = paste0(loc_penalty, "\n(n=", nLoc, ")"),
         method1 = paste0(method, " (n=", nOptimLoc, ")"))

## Add an annotation df
optimized_fitness_scaled_ann <- optimized_fitness_scaled1 %>%
  filter(loc_penalty == last(loc_penalty), str_detect(method, "Best"), fitness_component == "varY") %>%
  distinct_at(vars(nursery, management, contains("loc_penalty"), contains("method"), fitness_component, value))


g_optimized_fitness_comparison <- ggplot(data = NULL, aes(x = loc_penalty1, y = value)) +
  geom_hline(data = distinct(filter(optimized_fitness_scaled1, str_detect(method, "Best")), nursery, method, fitness_component, value),
             aes(yintercept = value, lty = method)) +
  geom_boxplot(data =  filter(optimized_fitness_scaled1, method == "Random"), aes(fill = "Random"), 
               lwd = 0.25, width = 0.5, outlier.shape = NA) +
  geom_point(data = filter(optimized_fitness_scaled1, method == "Optimization"), aes(shape = method), size = 1.25) +
  geom_text(data = optimized_fitness_scaled_ann, aes(label = method1, y = 1.1), size = 2, hjust = 1, vjust = 1, nudge_x = 0.5) +
  facet_grid(fitness_component ~ nursery, scales = "free", switch = "y",
             labeller = labeller(nursery = toupper, fitness_component = f_fitness_comp_replace)) +
  scale_y_continuous(name = NULL, breaks = pretty) +
  scale_x_discrete(name = "Location cost penalty\n(number of selected locations)") +
  scale_shape_discrete(name = NULL, guide = guide_legend(order = 1)) +
  scale_fill_manual(name = NULL, values = alpha("grey85", 0.5)) +
  scale_linetype_manual(name = NULL, values = 2, guide = FALSE) +
  theme_acs2() +
  theme(legend.background = element_rect(fill = alpha("white", 0)), legend.box.background = element_rect(fill = alpha("white", 0)),
        strip.placement = "outside", legend.box = "vertical", legend.spacing.y = unit(-0.1, "lines"), 
        legend.position = c(0.01, 0.69), legend.justification = c(0, 0),
        legend.key.height = unit(0.5, "line"), legend.key.width = unit(0.5, "line"))

# Save
ggsave(filename = "figure3_optimized_fitness_comparison_partA.jpg", path = fig_dir, plot = g_optimized_fitness_comparison,
       height = 9, width = 8, units = "cm", dpi = 1000)



## Part B - individual fitness for select traits

# Plot the individual fitness components for traits
optimized_fitness_trait <- optimized_locations_all_traits %>%
  # Count the number of locations
  mutate(nOptimLoc = map(optim_loc, "agro_loc") %>% map_dbl(length)) %>%
  # Bind rows
  mutate_at(vars(contains("fitness")), ~map(., bind_rows)) %>%
  # Convert varY such that higher values are favorable; only for the scaled version
  mutate_at(vars(contains("fitness_scaled")), ~map(., ~mutate(.x, varY = 1 - varY))) %>%
  # Unnest the desired column
  unnest(fitness) %>% mutate(varY = log10(varY)) %>%
  gather(fitness_component, value, varY, repAvg, reprAvg) %>%
  arrange(nursery, management, loc_penalty, trait_weight_group, method, fitness_component) %>% # Order
  mutate(fitness_component = fct_inorder(fitness_component) %>% fct_relevel("varY"),
         method = str_replace_all(method, "_", " ") %>% str_to_title() %>% str_remove_all(" "),
         method = ifelse(method == "BestRank", "Best ranked locations", method),
         method = fct_inorder(method) %>% fct_relevel("Optimization"))

# Subset to plot
optimized_fitness_trait1 <- optimized_fitness_trait %>%
  filter(trait_weight_group == "equal", comp_weight_group == "equal") %>%
  filter(trait %in% figure_traits1)



## Only show the zero penalty and the selected penalty
select_loc_penalty <- 0.01
optimized_fitness_trait2 <- optimized_fitness_trait1 %>%
  filter(loc_penalty %in% c(0, select_loc_penalty))


# Create a factor that combines the method and location penalty level; assign
# a color to each level of that factor;
# name levels with the method and the number of selected locations
group_factor_colors <- optimized_fitness_trait2 %>%
  distinct(nursery, method, loc_penalty, nOptimLoc, scale_components) %>%
  # Random nOptimLoc is the same as optimization at the same location penalty
  split(list(.$nursery, .$loc_penalty)) %>%
  map_dfr(~{
    mutate(.x, nOptimLoc = ifelse(method == "Random", subset(.x, method == "Optimization", nOptimLoc, drop = TRUE), nOptimLoc))
  }) %>%
  # Boxplot fill colors
  mutate(loc_penalty1 = ifelse(str_detect(method, "Best"), "none", loc_penalty),
         boxplot_fill = as.numeric(fct_inorder(loc_penalty1)),
         boxplot_fill = c(as.character(NA), neyhart_palette("umn2")[-1:-2])[boxplot_fill],
         boxplot_color = ifelse(is.na(boxplot_fill), "white", "black")) %>%
  # Point color and shape
  mutate(color_group = case_when(
    method == "Optimization" ~ paste0(method, " (penalty = ", loc_penalty,"; n = ", nOptimLoc, ")"),
    TRUE ~ paste0(method, " (n = ", nOptimLoc, ")")),
    color_group = fct_inorder(color_group),
    point_color = ifelse(is.na(boxplot_fill), "black", boxplot_fill),
    point_shape = color_group
  ) %>%
  select(-loc_penalty1, -nOptimLoc)
  
 
# Quantiles of random output
rand_quantiles <- c(0.10, 0.90)


## Split by nursery and create plots
optimized_fitness_comparison_trait_plotlist <- optimized_fitness_trait2 %>%
  group_by(scale_components, nursery, management) %>%
  nest() %>%
  ungroup() %>%
  arrange(scale_components, nursery, management) %>%
  mutate(n = seq_len(nrow(.)),
         plot = list(NULL))

# Iterate over rows
for (i in seq_len(nrow(optimized_fitness_comparison_trait_plotlist))) {
  row <- optimized_fitness_comparison_trait_plotlist[i,]
  
  df <- row$data[[1]] %>%
    crossing(., select(row, nursery, management, scale_components))
  n <- row$n
  
  # Create a custom grouping variable
  df_split <- df %>%
    mutate(method = as.character(method),
           # method = ifelse(method == "BestRank", method1, method),
           method = fct_relevel(method, "Best ranked locations", "Optimization", "Random")) %>%
    arrange(trait, method, loc_penalty) %>%
    mutate(group = ifelse(str_detect(method, "Best"), as.character(method), as.character(loc_penalty)),
           group = paste(trait, group, sep = "_"),
           group = fct_inorder(group),
           fitness_component = fct_inorder(f_fitness_comp_replace2(fitness_component))) %>%
    # Split the fitness component
    split(.$fitness_component)
  
  # Empty list of plots
  subplot_list <- list()
  
  # Iterate over fitness components for the subplots
  for (j in seq_along(df_split)) {
    df1 <- df_split[[j]]
    
    # Boxplots
    df1_boxplot <- filter(df1, method != "Optimization", !(str_detect(method, "Best") & loc_penalty == last(loc_penalty)))  %>%
      left_join(., group_factor_colors)
    # Vector of boxplot colors and fill
    boxplot_colors <- df1_boxplot %>%
      distinct(fitness_component, trait, group, boxplot_color) %>%
      arrange(fitness_component, trait, group) %>%
      pull()
    boxplot_fills <- df1_boxplot %>%
      distinct(fitness_component, trait, group, boxplot_fill) %>%
      arrange(fitness_component, trait, group) %>%
      pull()
    
    
    # Points
    df1_points <- filter(df1, method != "Random") %>%
      left_join(., group_factor_colors) %>%
      distinct(trait, value, group, point_shape, point_color, nursery, fitness_component)
    point_colors <- df1_points %>%
      distinct(point_shape, point_color) %>%
      {setNames(object = .$point_color, nm = .$point_shape)}
    
    # Plot
    g_sub <- ggplot(data = NULL, aes(x = trait, y = value, group = group)) +
      # Place an invisible point at 1
      geom_point(data = mutate(distinct(df1, nursery, trait, fitness_component, group), value = 1),
                 shape = NA, key_glyph = "polygon", position = position_dodge(0.75)) +
      geom_boxplot(data = df1_boxplot, position = position_dodge(0.90), alpha = 0.25, width = 0.75,
                   outlier.shape = NA, lwd = 0.25, color = boxplot_colors, fill = boxplot_fills) +
      geom_point(data = df1_points, aes(shape = point_shape, color = point_shape),
                 size = 1.5, position = position_dodge(0.90)) +
      facet_grid(fitness_component ~ nursery, scales = "free", switch = "y",
                 labeller = labeller(nursery = toupper, fitness_component = label_parsed)) +
      scale_y_continuous(name = NULL, breaks = pretty) +
      scale_x_discrete(name = NULL, labels = str_add_space, guide = guide_axis(check.overlap = TRUE, n.dodge = 2)) +
      scale_color_manual(values = point_colors, name = NULL) +
      scale_shape_discrete(name = NULL) +
      # scale_shape_discrete(guide = FALSE, labels = str_add_space) +
      theme_acs2() +
      # theme(legend.background = element_rect(fill = alpha("white", 0)), legend.box.background = element_rect(fill = alpha("white", 0))) +
      theme(strip.placement = "outside", legend.box = "horizontal", legend.spacing.y = unit(-0.1, "lines"),
            legend.direction = "vertical", legend.position = "top", legend.justification = c(0, 0),
            legend.box.just = "right", legend.title = element_text(size = 8), legend.key.height = unit(0.5, "line"),
            legend.key.width = unit(0.5, "line"))
    
    # If the fitness component is VarY, reverse the y axis scale
    if (str_detect(names(df_split)[j], "Precision")) g_sub <- g_sub + scale_y_continuous(name = NULL, breaks = pretty, trans = "reverse") 

    # If j is 1, keep the legend and the top facet strip
    if (j == 1) {
      g_sub <- g_sub +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      
    } else if (j == 2) {
      # If j is 2, remove everything
      g_sub <- g_sub +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              legend.position = "none", strip.text.x = element_blank())
      
    } else {
      # If j is 3, keep only the x axis
      g_sub <- g_sub +
        theme(legend.position = "none", strip.text.x = element_blank())
      
    }
    
    # Add to the list
    subplot_list[[j]] <- g_sub
    
  } # Close the subplot loop
  
  # If n is odd, remove the y axis strip
  # For odd values of n
  if (n %% 2 == 0) {
    subplot_list <- map(subplot_list, ~. + theme(strip.text.y = element_blank()))
    
  }
  
  # Combine the subplots
  g <- patchwork::wrap_plots(... = subplot_list, ncol = 1)
  
  optimized_fitness_comparison_trait_plotlist$plot[[i]] <- g
  
}

# Merge plots together
g_combine_list <- optimized_fitness_comparison_trait_plotlist %>%
  mutate(scale_components = paste0("scale_components", scale_components)) %>%
  split(.$scale_components) %>%
  map(~patchwork::wrap_plots(.x$plot, nrow = 1))

# Save
for (i in seq_along(g_combine_list)) {
  ggsave(filename = paste0("figure3_optimized_fitness_comparison_partB_", names(g_combine_list)[i], ".jpg"), 
         path = fig_dir, plot = g_combine_list[[i]], height = 12, width = 12, units = "cm", dpi = 1000)
}



## Merge plots
g_combine_all <- plot_grid(g_optimized_fitness_comparison, g_combine_list$scale_componentsTRUE, 
                           nrow = 1, align = "h", axis = "tblr", rel_widths = c(0.45, 0.68), labels = letters[1:2], 
                           label_size = 10, label_y = 0.90)


ggsave(filename = "figure3_optimized_fitness_comparison.jpg", path = fig_dir, plot = g_combine_all,
       height = 11, width = 17.5, units = "cm", dpi = 1000)




# Figure 4: Map of selected locations --------------------------------



## Map the selected locations for the penalty of 0.01
selected_locations <- optimized_locations_all_traits %>%
  filter(method == "optimization", loc_penalty == 0.01, comp_weight_group == "equal", 
         trait_weight_group == "equal", scale_components) %>% 
  mutate(optim_loc_df = map(optim_loc, ~as.data.frame(table(unlist(.x))))) %>%
  unnest(optim_loc_df) %>%
  select(-fitness:-optim_loc) %>%
  rename(location = Var1, freq = Freq) %>%
  mutate(trait_group = ifelse(freq == 2, "Malt quality + agronomic traits", "Agronomic traits"),
         nursery_abbr = ifelse(nursery == "Both", nursery, toupper(nursery)),
         nursery = f_nursery_expand(nursery)) %>%
  left_join(., select(trial_info_toplot1, -nursery, -nursery_abbr)) %>%
  mutate(nursery = str_remove_all(nursery, " Nursery"))

# Add remaining locations
selected_locations_toplot <- selected_locations %>%
  bind_rows(., filter(trial_info_toplot1, ! location %in% selected_locations$location)) %>%
  mutate(trait_group = ifelse(is.na(trait_group), "unselected", trait_group),
         selected = trait_group != "unselected",
         nursery_abbr = fct_inorder(nursery_abbr)) %>%
  select(nursery_abbr, location, lat, long, trait_group, selected)


# Reset x and y lims
xlim <- c(-122, -85)
ylim <- c(35, 51)

# First build the base map plot
g_base_map <- ggplot(data = north_america, aes(x = long, y = lat)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = canada, aes(group = group), fill = "grey95", color = "grey50", lwd = 0.5) + # Add canada
  geom_polygon(data = usa_state, aes(group = group), fill = "grey95", color = "grey50", lwd = 0.5) +
  coord_map(projection = "bonne", lat0 = mean(ylim), xlim = xlim, ylim = ylim) +
  theme_void(base_size = 10)
  
# Add points
g_map1 <- g_base_map + 
  geom_point(data = selected_locations_toplot, aes(color = nursery_abbr, shape = trait_group, size = selected)) +
  geom_label_repel(data = subset(selected_locations_toplot, selected), aes(label = location, color = nursery_abbr),
                   fill = alpha("white", 0.75), point.padding = unit(0.5, "line"),
                   min.segment.length = 0, force = 10, key_glyph = "blank", size = 2.5) +
  scale_color_manual(values = neyhart_palette("umn1")[c(3,4,2)], name = NULL) +
  # scale_color_viridis_d(begin = 0.2, end = 0.8, name = NULL) +
  scale_shape_manual(name = NULL, values = c(15, 17, 16), breaks = setdiff(unique(selected_locations_toplot$trait_group), "unselected")) +
  scale_size_manual(values = c(1, 2), guide = FALSE) +
  # Adjust the orientation of guide
  guides(color = guide_legend(order = 1, direction = "horizontal", keywidth = unit(0.25, "line")),
         shape = guide_legend(order = 2, direction = "vertical", keyheight = unit(0.75, "line"), keywidth = unit(1, "line"))) +
  theme(legend.box.background = element_rect(fill = "white", color = "white"),
        legend.position = c(1, 0), legend.box = "vertical", legend.box.just = "left",
        legend.justification = c(1, 0), legend.box.margin = margin(), legend.spacing = unit(0, "line"), 
        legend.text = element_text(size = 8))


# Save the figure
ggsave(filename = "figure4_trial_location_map_optimized.jpg", plot = g_map1, path = fig_dir,
       width = 8.5, height = 6, units = "cm", dpi = 1000)









# Table 1: data summary ---------------------------------------------------

pheno_dat %>%
  left_join(., select(trial_metadata, trial, nursery, management)) %>%
  filter(!(nursery == "wrn" & management == "rainfed")) %>%
  summarize_at(vars(location, trial, year, line_name, trait), n_distinct)

# location trial  year line_name trait
#       46   426    25       722    20

# Summarize data available for each nursery (exclude rainfed WRN)
(pheno_dat_summary <- pheno_dat %>%
  left_join(., select(trial_metadata, trial, nursery, management)) %>%
  filter(!(nursery == "wrn" & management == "rainfed")) %>%
  group_by(nursery) %>%
  summarize_at(vars(location, environment, year, line_name, trait), n_distinct))

# nursery location environment  year line_name trait
# 1 mvn           20         175    25       401    19
# 2 wrn           32         251    23       393    18
    

# Save
pheno_dat_summary %>%
  mutate(nursery = f_nursery_expand(nursery)) %>%
  rename(genotypes = line_name) %>%
  rename_all(str_to_title) %>%
  write_csv(x = ., file = file.path(fig_dir, "table1_data_summary.csv"))

## Data used in analysis
pheno_dat %>%
  left_join(., select(trial_metadata, trial, nursery, management)) %>%
  filter(!(nursery == "wrn" & management == "rainfed")) %>%
  group_by(trait, nursery, management) %>%
  summarize_at(vars(trial, location, year), n_distinct) %>%
  arrange(location) %>% as.data.frame()





# Supplemental Figure 1: location precision (varR) -------------------

# Plot genotypic variance for the figure traits
precision_plot_list <- environmental_precision_toplot %>%
  group_by(trait, nursery, management) %>%
  do({
    df <- .
    
    # Find outliers
    df1 <- df %>%
      filter(varY <= mean(varY) + (sd(varY)*3))
      
    nur <- unique(df1$nursery)
    mgmt <- unique(df1$management)
    tr <- unique(df1$trait) %>% ifelse(. == "SolubleProteinTotalProtein", f_trait_abbr(.), .)
    
    plot <- ggplot(df1, aes(x = x, y = log10(varY))) + 
      geom_jitter(width = 0.1, color = "grey85", size = 0.75) +
      geom_boxplot(alpha = 0, outlier.shape = NA, lwd = 0.25) +
      scale_y_continuous(name = expression('Variance of a genotype mean ['*log[10](italic(V[bar(Y)]))*']'), breaks = pretty, 
                         labels = function(x) ifelse(x > 1e4, formatC(x, format = "e", digits = 1), x)) +
      scale_x_discrete(name = NULL, labels = function(x) map_chr(str_split(x, "\\."), 3)) +
      labs(subtitle = paste0(c(f_trait_rename(tr), toupper(nur)), collapse = ", ")) +
      theme_genetics(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), axis.text.y = element_text(size = 6), strip.placement = "outside")
    
    tibble(plot = list(plot), nRemoved = nrow(df) - nrow(df1))
    
  }) %>% ungroup()
    

# Combine plots
plot_combined_list <- precision_plot_list %>%
  split(.$nursery) %>%
  map(~{
    
    nursery <- unique(.x$nursery)
    
    ## Modify plots
    plot_list1 <- map(.x$plot, ~. + theme(axis.title = element_blank()))
    
    # Add y axis grob
    y_axis <- grid::textGrob(label = expression('Variance of a genotype mean ['*log[10](italic(V[bar(Y)]))*']'), rot = 90, 
                             gp = gpar(fontsize = 10))

    # Merge plots
    plot_grid(y_axis, add_sub(plot = plot_grid(plotlist = plot_list1, ncol = 4), label = "Location", size = 10), rel_widths = c(0.05, 1))
    
  })

# Merge A and B together
g_combined <- plot_grid(plotlist = plot_combined_list, ncol = 1, labels = letters[seq_along(plot_combined_list)],
                        label_size = 9)

# Save
ggsave(filename = "figureS1_precision_boxplot_by_location.jpg", plot = g_combined, path = fig_dir,
       height = 10, width = 8, dpi = 1000)











# Supplemental Figure 2: location repeatability ----------------------


# Plot
location_cor_plot_list <- location_correlations %>%
  group_by(nursery, management, trait) %>%
  do(plot = {
    row <- .
    
    nur <- row$nursery
    mgmt <- row$management
    tr <- row$trait %>% ifelse(. == "SolubleProteinTotalProtein", f_trait_abbr(.), .)
    
    df <- row$correlations[[1]] %>%
      # Filter for matching xy locations
      filter(location.x == location.y) %>%
      mutate(delim = map2_chr(environment.x, environment.y, ~paste0(sort(c(.x, .y)), collapse = ":")),
             trait = tr, nursery = nur, management = mgmt) %>%
      filter(!duplicated(delim)) %>%
      select(-delim)
    
    
    # Reorder location factors - create annotation
    df %>%
      left_join(., location_abbr_key, by = c("location.x" = "location")) %>%
      mutate(location_abbr = fct_reorder(location_abbr, .x = correlation, .fun = median, .desc = TRUE)) %>%
      ggplot(aes(x = location_abbr, y = correlation)) +
      geom_jitter(width = 0.25, color = "grey85", size = 0.5) +
      geom_boxplot(alpha = 0, outlier.shape = NA, lwd = 0.25) +
      scale_y_continuous(name = expression("Repeatability (within-location"~italic(rho[G])*")"), breaks = pretty) +
      scale_x_discrete(name = "Location") +
      labs(subtitle = paste0(c(f_trait_rename(tr), toupper(nur)), collapse = ", ")) +
      theme_genetics(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), axis.text.y = element_text(size = 6), strip.placement = "outside")
    
  }) %>% 
  ungroup()

# Combine plots
plot_combined_list <- location_cor_plot_list %>%
  split(.$nursery) %>%
  map(~{
    
    nursery <- unique(.x$nursery)
    
    ## Modify plots
    plot_list1 <- map(.x$plot, ~. + theme(axis.title = element_blank()))
    
    # Add y axis grob
    y_axis <- grid::textGrob(label = expression("Repeatability (within-location"~italic(rho[G])*")"), rot = 90, 
                             gp = gpar(fontsize = 10))
    
    # Merge plots
    plot_grid(y_axis, add_sub(plot = plot_grid(plotlist = plot_list1, ncol = 4), label = "Location", size = 10), rel_widths = c(0.05, 1))
    
  })

# Merge A and B together
g_combined <- plot_grid(plotlist = plot_combined_list, ncol = 1, labels = letters[seq_along(plot_combined_list)],
                        label_size = 9)

# Save
ggsave(filename = "figureS2_location_correlations.jpg", plot = g_combined, path = fig_dir,
       height = 10, width = 8, dpi = 1000)




# Supplemental Figure 3: mega-environment assignment -----------------

## Save mega-environment plots
plot_combined_list <- nursery_mega_environment_output %>%
  # Add proportions of variance
  left_join(., loadings_prop_var) %>%
  mutate(delim = paste(nursery, management, sep = "_"),
         plot = pmap(select(., plot, fa1, fa2), ~{
           .x <- ..1; .y <- ..2; .z <- ..3
           .x$labels$subtitle <- f_trait_rename(.x$labels$subtitle) %>% str_remove(", Rainfed|, Irrigated")
           .x + 
             scale_x_continuous(name = paste0("Loading 1 (", round(.y * 100, 1), "%)"), breaks = pretty) +
             scale_y_continuous(name = paste0("Loading 2 (", round(.z * 100, 1), "%)"), breaks = pretty)
           })
  ) %>%
  split(.$delim) %>%
  # Create a plot grid
  map(~plot_grid(plotlist = .x$plot, ncol = 4))

# Merge A and B together
g_combined <- plot_grid(plotlist = plot_combined_list, ncol = 1, labels = letters[seq_along(plot_combined_list)],
                        label_size = 9)

# Save
ggsave(filename = "figureS3_mega_environments.jpg", plot = g_combined, path = fig_dir,
       height = 18, width = 14, dpi = 1000)
  

## summarize varprop of loadings per trait/ nursery
loadings_prop_var_total <- loadings_prop_var %>%
  group_by(trait, nursery, management) %>%
  summarize(prop_var_loadings = fa1 + fa2) %>%
  ungroup() %>%
  arrange(desc(prop_var_loadings))




# Supplemental Figure 4: representativeness ----------------



# Plot
location_repr_plot_list <- representativeness_tidy %>%
  group_by(nursery, management, trait) %>%
  nest() %>%
  rename(correlations = data) %>%
  do(plot = {
    row <- .
    
    nur <- row$nursery
    mgmt <- row$management
    tr <- row$trait %>%
      ifelse(. == "SolubleProteinTotalProtein", f_trait_abbr(.), .)
    
    df <- row$correlations[[1]] %>%
      mutate(delim = map2_chr(environment.x, environment.y, ~paste0(sort(c(.x, .y)), collapse = ":")),
             trait = tr, nursery = nur, management = mgmt) %>%
      filter(!duplicated(delim)) %>%
      select(-delim)
    
    # Create a title
    title_text <- parse(text = paste0("'", f_trait_rename(tr), ", ", toupper(nur), ", '*", 
                                      paste0("N[ME]==", n_distinct(df$me))))
      
      
    
    # Reorder location factors - create annotation
    df %>%
      mutate(location_abbr = fct_reorder(location_abbr, .x = correlation, .fun = median, .desc = TRUE)) %>%
      ggplot(aes(x = location_abbr, y = correlation)) +
      geom_jitter(width = 0.25, color = "grey85", size = 0.5) +
      geom_boxplot(aes(color = me), alpha = 0, outlier.shape = NA, lwd = 0.25) +
      scale_y_continuous(name = "Mega-environment representativeness", breaks = pretty) +
      scale_x_discrete(name = "Location") +
      scale_color_discrete(guide = FALSE) +
      labs(subtitle = title_text) +
      theme_genetics(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6), strip.placement = "outside")
    
  }) %>% 
  ungroup()

# Combine plots
plot_combined_list <- location_repr_plot_list %>%
  split(.$nursery) %>%
  map(~{
    
    nursery <- unique(.x$nursery)
    
    ## Modify plots
    plot_list1 <- map(.x$plot, ~. + theme(axis.title = element_blank()))
    
    # Add y axis grob
    y_axis <- grid::textGrob(label = "Mega-environment representativeness", rot = 90, 
                             gp = gpar(fontsize = 10))
    
    # Merge plots
    plot_grid(y_axis, add_sub(plot = plot_grid(plotlist = plot_list1, ncol = 4), label = "Location", size = 10), rel_widths = c(0.05, 1))
    
  })

# Merge A and B together
g_combined <- plot_grid(plotlist = plot_combined_list, ncol = 1, labels = letters[seq_along(plot_combined_list)],
                        label_size = 9)

# Save
ggsave(filename = "figureS4_location_me_representativeness.jpg", plot = g_combined, path = fig_dir,
       height = 10, width = 8, dpi = 1000)



# Supplemental Figure 5: map of mega-environments -------------------------


# Create a box for x and y axes
xlim <- range(pretty(trial_info_toplot$long))
ylim <- range(pretty(trial_info_toplot$lat))

# Distinct df
me_location_assignment <- representativeness_tidy %>%
  distinct(trait, nursery, management, location, me)

# Number of MEs per trait/nursery
aggregate(formula = me ~ trait + nursery, data = me_location_assignment, FUN = n_distinct) %>% arrange(me)


# Alternative supplemental figure
# 
# Plot a map and then locations by mega-environment
set.seed(1220)
g_map_megaEnvironments_list <- me_location_assignment %>%
  split(.$trait) %>%
  map(~{
    df1 <- .
    # # unnest
    # df1 <- unnest(df1, cols = c(me_data)) %>%
    #   select(-plot)
    
    df2 <- left_join(df1, trial_info_toplot1, by = "location") %>%
      unite(me_factor, c(nursery.x, me), remove = FALSE) %>%
      rename(nursery = nursery.x) %>%
      mutate(nursery_abbr = toupper(nursery))
    
    # Separate df for points only in one nursery versus those that are in two
    df2_single_nursery <- df2 %>% 
      group_by(location) %>% 
      filter(n() == 1) %>% 
      ungroup()
    df2_multiple_nursery <- setdiff(df2, df2_single_nursery) %>%
      mutate(nursery_abbr = toupper(nursery)) %>%
      mutate_at(vars(lat, long), list(end = ~jitter(., factor = 1)))
    
    ggplot(data = north_america, aes(x = long, y = lat)) +
      geom_polygon(fill = "white") +
      geom_polygon(data = canada, aes(group = group), fill = "grey95", color = "grey50", lwd = 0.75) + # Add canada
      geom_polygon(data = usa_state, aes(group = group), fill = "grey95", color = "grey50", lwd = 0.75) +
      geom_point(data = df2_single_nursery, aes(shape = nursery_abbr, color = me_factor), size = 3.5) +
      # Add segments connecting the jittered points to the original point
      geom_segment(data = df2_multiple_nursery, aes(xend = long_end, yend = lat_end)) +
      geom_point(data = df2_multiple_nursery, aes(x = long_end, y = lat_end, shape = nursery_abbr, color = me_factor), 
                 size = 3.5) +
      scale_shape_discrete(name = NULL) +
      scale_color_discrete(guide = FALSE) +
      labs(subtitle = str_add_space(unique(df2$trait))) +
      coord_map(projection = "bonne", lat0 = mean(ylim), xlim = xlim, ylim = ylim) +
      theme_void(base_size = 14) +
      theme(legend.position = c(0.20, 0.10), legend.background = element_rect(fill = "white"),
            legend.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "line"))
    
  })


# Create a grid of plots
g_combined <- plot_grid(plotlist = g_map_megaEnvironments_list, ncol = 3)

# Save
ggsave(filename = "figureS5_mega_environments_location_maps.jpg", plot = g_combined, path = fig_dir,
       height = 20, width = 18, dpi = 1000)




# Supplemental Figure 6: overall location performance ----------------

# Combine plots within nurseries; omit x and y axes and legends
location_performance_plots1 <- location_performance_plots %>%
  mutate(delim = paste(nursery, management, sep = "_"),
         plot = map(plot, ~{
           .x$labels$subtitle <- f_trait_rename(.x$labels$subtitle) %>% str_remove(", Rainfed|, Irrigated")
           # Change positions of legends to bottom
           .x + 
             theme(legend.position = "bottom", legend.key.width = unit(0.5, "lines")) +
             guides(size = guide_legend(label.position = "top"),
                    shape = guide_legend(label = FALSE))
           
         }) )

# Extract the plots, remove axes and legends
location_performance_plots_list <- location_performance_plots1 %>%
  split(.$delim) %>%
  map("plot") %>%
  map(~ map(.x, ~. + theme(legend.position = "none", axis.title = element_blank())) )

# Create the x and y axes
y_axis <- get_plot_component(location_performance_plots1$plot[[1]], "ylab-l")
x_axis <- get_plot_component(location_performance_plots1$plot[[1]], "xlab-b")
# Get the legend from the plot with the most MEs
legend_list <- location_performance_plots1 %>% 
  mutate(nME = map(plot, "data") %>% map_dbl(~n_distinct(.x$mega_environment))) %>% 
  split(.$delim) %>%
  map(~{filter(.x, nME == max(nME)) %>% 
      head(1) %>% 
      pull(plot) %>% 
      first() %>%
      get_legend()
  })

## combine plots by nursery
plot_combined_list <- location_performance_plots_list %>%
  map(~plot_grid(plotlist = .x, ncol = 4)) %>%
  # Add the axes
  map(~plot_grid(y_axis, .x, nrow = 1, rel_widths = c(0.05, 1))) %>%
  map2(., legend_list, ~plot_grid(.x, x_axis, .y, ncol = 1, rel_heights = c(1, 0.05, 0.07)))

# Merge A and B together
g_combined <- plot_grid(plotlist = plot_combined_list, ncol = 1, labels = letters[seq_along(plot_combined_list)],
                        label_size = 9)

# Save
ggsave(filename = "figureS6_location_overall_performance.jpg", plot = g_combined, path = fig_dir,
       height = 12, width = 10, dpi = 1000)


# Supplemental Figure 7: optimization sensitivity analysis --------------

# Plot the results of the summarized fitness values
g_sensitivity_model_effects_list <- sensitivity_test_fitness_summary_fitted_models_effects1 %>%
  mutate(metricWeight = fct_rev(str_remove(metricWeight, "_weight"))) %>%
  split(.$nursery) %>%
  map(~{
    ggplot(data = .x, aes(x = weight, y = fit)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey85", color = "black",
                  alpha = 0.5, lty = 2) +
      geom_line(color = "blue") +
      geom_text(data = distinct(.x, metricWeight, metric, annotation),
                aes(x = Inf, y = -Inf, label = annotation), hjust = 1.1, vjust = -0.5, size = 3, parse = TRUE) +
      facet_grid(metric ~ metricWeight, scales = "free_y", switch = "both",
                 labeller = labeller(metric = f_fitness_comp_replace, metricWeight = f_fitness_comp_replace)) +
      scale_y_continuous(breaks = pretty, name = "Scaled value of each metric") +
      scale_x_continuous(name = "Weight assigned to each metric in the objective function") +
      labs(subtitle = f_nursery_expand(x = unique(.x$nursery))) +
      theme_acs() +
      theme(strip.placement = "outside", strip.background = element_blank())
  })

# Combine
g_combine1 <- plot_grid(plotlist = g_sensitivity_model_effects_list, nrow = 1)




# Plot
g_sensitivity_model_effects_list2 <- sensitivity_test_fitness_fitted_models_effects1 %>%
  distinct(metric, nursery, trait, metricWeight, estimate, signif_asterick) %>%
  mutate(metricWeight = str_remove(metricWeight, "_weight")) %>%
  split(.$nursery) %>%
  map(~{
    ggplot(data = .x, aes(x = metric, y = estimate, fill = metricWeight)) +
      geom_col(position = position_dodge(0.90)) +
      geom_text(aes(label = signif_asterick), position = position_dodge(0.90), size = 3) +
      facet_wrap(~ trait, labeller = labeller(trait = str_add_space)) +
      scale_y_continuous(name = "Effect size", breaks = pretty) +
      scale_x_discrete(labels = f_fitness_comp_replace, name = "Metric responding to varying weight") +
      scale_fill_discrete(labels = f_fitness_comp_replace, name = "Metric varying\nin weight",
                          guide = guide_legend(title.position = "left")) +
      labs(subtitle = f_nursery_expand(x = unique(.x$nursery))) +
      theme_acs2() +
      theme(legend.position = c(0.90, 0.10), axis.text.x = element_text(angle = 30, hjust = 1))
  })

# Combine
g_combine2 <- plot_grid(plotlist = g_sensitivity_model_effects_list2, nrow = 1)


# Combine and save

# Save
ggsave(filename = "figureS7_optimization_sensitivity_analysis.jpg", 
       plot = plot_grid(g_combine1, g_combine2, ncol = 1, labels = letters[1:2], label_size = 9), 
       path = fig_dir, height = 12, width = 16, dpi = 1000)


# Supplemental Figure 8 and 9: optimization results for all traits --------------

# Subset to plot
optimized_fitness_trait2 <- optimized_fitness_trait %>%
  filter(comp_weight_group == "equal")

# Create a factor that combines the method and location penalty level; assign
# a color to each level of that factor;
# name levels with the method and the number of selected locations
group_factor_colors <- optimized_fitness_trait2 %>%
  distinct(nursery, method, loc_penalty, nOptimLoc, scale_components, trait_weight_group) %>%
  # Random nOptimLoc is the same as optimization at the same location penalty
  split(list(.$nursery, .$loc_penalty)) %>%
  map_dfr(~{
    mutate(.x, nOptimLoc = ifelse(method == "Random", subset(.x, method == "Optimization", nOptimLoc, drop = TRUE), nOptimLoc))
  }) %>%
  # Boxplot fill colors
  mutate(loc_penalty1 = ifelse(str_detect(method, "Best"), "none", loc_penalty),
         boxplot_fill = as.numeric(fct_inorder(loc_penalty1)),
         boxplot_fill = c(as.character(NA), neyhart_palette("umn2")[-1:-2])[boxplot_fill],
         boxplot_color = ifelse(is.na(boxplot_fill), "white", "black")) %>%
  # Point color and shape
  mutate(color_group = case_when(
    method == "Optimization" ~ paste0(method, " (penalty = ", loc_penalty,"; n = ", nOptimLoc, ")"),
    TRUE ~ paste0(method, " (n = ", nOptimLoc, ")"))) %>%
  arrange(desc(method), loc_penalty) %>%
  mutate(color_group = fct_inorder(color_group),
         point_color = ifelse(is.na(boxplot_fill), "black", boxplot_fill),
         point_shape = color_group) %>%
  select(-loc_penalty1, -nOptimLoc)


# Quantiles of random output
rand_quantiles <- c(0.10, 0.90)


## Split by nursery and create plots
optimized_fitness_comparison_trait_plotlist <- optimized_fitness_trait2 %>%
  group_by(scale_components, trait_weight_group, nursery, management) %>%
  nest() %>%
  ungroup() %>%
  arrange(scale_components, trait_weight_group, nursery, management) %>%
  mutate(n = seq_len(nrow(.)),
         plot = list(NULL))

# Iterate over rows
for (i in seq_len(nrow(optimized_fitness_comparison_trait_plotlist))) {
  row <- optimized_fitness_comparison_trait_plotlist[i,]
  
  df <- row$data[[1]] %>%
    crossing(., select(row, nursery, management, scale_components, trait_weight_group))
  n <- row$n
  
  # Create a custom grouping variable
  df_split <- df %>%
    mutate(method = as.character(method),
           # method = ifelse(method == "BestRank", method1, method),
           method = fct_relevel(method, "Best ranked locations", "Optimization", "Random")) %>%
    arrange(trait, method, loc_penalty) %>%
    mutate(group = ifelse(str_detect(method, "Best"), as.character(method), as.character(loc_penalty)),
           group = paste(trait, group, sep = "_"),
           group = fct_inorder(group),
           fitness_component = fct_inorder(f_fitness_comp_replace2(fitness_component))) %>%
    # Split the fitness component
    split(.$fitness_component)
  
  # Empty list of plots
  subplot_list <- list()
  
  # Iterate over fitness components for the subplots
  for (j in seq_along(df_split)) {
    df1 <- df_split[[j]]
    
    # Boxplots
    df1_boxplot <- filter(df1, method != "Optimization", !(str_detect(method, "Best") & loc_penalty == last(loc_penalty)))  %>%
      left_join(., group_factor_colors)
    # Vector of boxplot colors and fill
    boxplot_colors <- df1_boxplot %>%
      distinct(fitness_component, trait, group, boxplot_color) %>%
      arrange(fitness_component, trait, group) %>%
      pull()
    boxplot_fills <- df1_boxplot %>%
      distinct(fitness_component, trait, group, boxplot_fill) %>%
      arrange(fitness_component, trait, group) %>%
      pull()
    
    
    # Points
    df1_points <- filter(df1, method != "Random") %>%
      left_join(., group_factor_colors) %>%
      distinct(trait, value, group, point_shape, point_color, nursery, fitness_component)
    point_colors <- df1_points %>%
      distinct(point_shape, point_color) %>%
      {setNames(object = .$point_color, nm = .$point_shape)}
    
    # Plot
    g_sub <- ggplot(data = NULL, aes(x = trait, y = value, group = group)) +
      # Place an invisible point at 1
      geom_point(data = mutate(distinct(df1, nursery, trait, fitness_component, group), value = 1),
                 shape = NA, key_glyph = "polygon", position = position_dodge(0.75)) +
      geom_boxplot(data = df1_boxplot, position = position_dodge(0.90), alpha = 0.25, width = 0.75,
                   outlier.shape = NA, lwd = 0.25, color = boxplot_colors, fill = boxplot_fills) +
      geom_point(data = df1_points, aes(shape = point_shape, color = point_shape),
                 size = 1.5, position = position_dodge(0.90)) +
      facet_grid(fitness_component ~ nursery, scales = "free", switch = "y",
                 labeller = labeller(nursery = toupper, fitness_component = label_parsed)) +
      scale_y_continuous(name = NULL, breaks = pretty) +
      scale_x_discrete(name = NULL, labels = str_add_space, guide = guide_axis(check.overlap = TRUE, n.dodge = 2)) +
      scale_color_manual(values = point_colors, name = NULL, guide = guide_legend(ncol = 2)) +
      scale_shape_discrete(name = NULL, guide = guide_legend(ncol = 2)) +
      # scale_shape_discrete(guide = FALSE, labels = str_add_space) +
      theme_acs2() +
      # theme(legend.background = element_rect(fill = alpha("white", 0)), legend.box.background = element_rect(fill = alpha("white", 0))) +
      theme(strip.placement = "outside", legend.box = "horizontal", legend.spacing.y = unit(-0.1, "lines"),
            legend.direction = "vertical", legend.position = "top", legend.justification = c(0, 0),
            legend.box.just = "right", legend.title = element_text(size = 8), legend.key.height = unit(0.5, "line"),
            legend.key.width = unit(0.5, "line"))
    
    # If the fitness component is VarY, reverse the y axis scale
    if (str_detect(names(df_split)[j], "Precision")) g_sub <- g_sub + scale_y_continuous(name = NULL, breaks = pretty, trans = "reverse") 
    
    # If j is 1, keep the legend and the top facet strip
    if (j == 1) {
      g_sub <- g_sub +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      
    } else if (j == 2) {
      # If j is 2, remove everything
      g_sub <- g_sub +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              legend.position = "none", strip.text.x = element_blank())
      
    } else {
      # If j is 3, keep only the x axis
      g_sub <- g_sub +
        theme(legend.position = "none", strip.text.x = element_blank())
      
    }
    
    # Add to the list
    subplot_list[[j]] <- g_sub
    
  } # Close the subplot loop
  
  # Combine the subplots
  g <- patchwork::wrap_plots(... = subplot_list, ncol = 1)
  
  optimized_fitness_comparison_trait_plotlist$plot[[i]] <- g
  
}

# Merge plots together
g_combine_list_df <- optimized_fitness_comparison_trait_plotlist %>%
  mutate(scale_components = paste0("scale_components", scale_components)) %>%
  group_by(scale_components, trait_weight_group) %>%
  do(plot = patchwork::wrap_plots(.$plot, ncol = 1)) %>%
  ungroup()

g_combine_list <- subset(g_combine_list_df, trait_weight_group == "equal", plot, drop = TRUE)
  
# Save
for (i in seq_along(g_combine_list)) {
  ggsave(filename = paste0("figureS8_optimized_fitness_comparison_partB_", names(g_combine_list)[i], ".jpg"), 
         path = fig_dir, plot = g_combine_list[[i]], height = 24, width = 18, units = "cm", dpi = 1000)
}


# Supplemental figure 7
ggsave(filename = "figureS8_optimized_fitness_comparison_all_traits.jpg", path = fig_dir, 
       plot = subset(g_combine_list_df, trait_weight_group == "equal" & scale_components == "scale_componentsTRUE", plot, drop = TRUE)[[1]],
       height = 24, width = 18, units = "cm", dpi = 1000)

# Supplemental figure 8
ggsave(filename = "figureS9_optimized_fitness_comparison_all_traits_weighted.jpg", path = fig_dir, 
       plot = subset(g_combine_list_df, trait_weight_group == "weighted" & scale_components == "scale_componentsTRUE", plot, drop = TRUE)[[1]],
       height = 24, width = 18, units = "cm", dpi = 1000)



# Supplemental Figure 10: precision by number of years of sampling --------------


# Summarize
varY_resampling_summary <- varY_resampling %>%
  unnest(out) %>%
  group_by(trait, nursery, management, location, nEnvSample) %>%
  summarize(varY = mean(varY), .groups = "drop") %>%
  # Add a boolean for locations selected via optimization
  mutate(nursery = toupper(nursery)) %>%
  left_join(., select(selected_locations, nursery = nursery_abbr, location, method)) %>%
  mutate(method = ifelse(is.na(method), "not_selected", method))


## Plot results
# Filter traits to highlight
varY_resampling_summary1 <- varY_resampling_summary %>%
  filter(trait %in% c("GrainProtein", "GrainYield")) %>%
  mutate(nursery = tolower(nursery)) 

# Display locations selected via optimization in bold
# 
g_vary_resampling <- varY_resampling_summary %>%
  filter(trait %in% c("GrainProtein", "GrainYield")) %>%
  mutate(nursery = tolower(nursery)) %>%
  ggplot(aes(x = nEnvSample, y = log10(varY), color = location)) +
  geom_line(aes(size = method)) +
  # geom_label(data = top_n(x = group_by(filter(varY_resampling_summary1, method == "optimization"), trait, nursery, location), n = 1, wt = nEnvSample),
  #            aes(label = location), force = 10) +
  geom_label_repel(data = top_n(x = group_by(filter(varY_resampling_summary1, method == "optimization"), trait, nursery, location), n = 1, wt = nEnvSample),
                   aes(x = Inf, y = log10(varY) + 0.5, label = location), min.segment.length = 100, direction = "y",
                   size = 3, fill = alpha("white", 0.75)) +
  scale_color_discrete(guide = FALSE) +
  scale_size_manual(values = c(0.25, 0.75), guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = expression('Precision ['*log[10](italic(V[bar(Y)]))*']'), trans = "reverse") + 
  scale_x_continuous(breaks = pretty, name = "Number of sampled years at a location") +
  facet_grid(trait ~ nursery, scales = "free_y", switch = "y", 
             labeller = labeller(trait = str_add_space, nursery = f_nursery_expand)) +
  theme_acs2(base_size = 8) +
  theme(strip.placement = "outside")


# Save
ggsave(filename = "figureS10_precision_location_resampling.jpg", plot = g_vary_resampling, path = fig_dir,
       height = 4, width = 6, dpi = 1000)



# # Supplemental Figure XX: compare optimization with or without scaling the metrics --------------
# 
# # List of plots
# g_optimized_fitness_comparison_list <- optimized_fitness_scaled_ann %>%
#   mutate(plot = pmap(.l = ., .f = ~{
#     subtitle <- paste0("Components scaled: ", ..1)
#     
#     ggplot(data = NULL, aes(x = loc_penalty1, y = value)) +
#       geom_hline(data = distinct(filter(..3, method == "BestRank"), nursery, method, fitness_component, value),
#                  aes(yintercept = value, lty = method)) +
#       geom_boxplot(data =  filter(..3, method == "Random"), aes(fill = "Random"), 
#                    lwd = 0.25, width = 0.5, outlier.shape = NA) +
#       geom_point(data = filter(..3, method == "Optimization"), aes(shape = method), size = 1.25) +
#       geom_text(data = ..2, aes(label = method1, y = 1.1), size = 2, hjust = 1, vjust = 1, nudge_x = 0.5) +
#       facet_grid(fitness_component ~ nursery, scales = "free", switch = "y",
#                  labeller = labeller(nursery = toupper, fitness_component = f_fitness_comp_replace)) +
#       scale_y_continuous(name = NULL, breaks = pretty) +
#       scale_x_discrete(name = "Location cost penalty\n(number of selected locations)") +
#       scale_shape_discrete(name = NULL, guide = guide_legend(order = 1)) +
#       scale_fill_manual(name = NULL, values = alpha("grey85", 0.5)) +
#       scale_linetype_manual(name = NULL, values = 2, guide = FALSE) +
#       labs(subtitle = subtitle) +
#       theme_acs() +
#       theme(legend.background = element_rect(fill = alpha("white", 0)), legend.box.background = element_rect(fill = alpha("white", 0)),
#             strip.placement = "outside", legend.box = "vertical", legend.spacing.y = unit(-0.1, "lines"), 
#             legend.position = c(0.01, 0.69), legend.justification = c(0, 0),
#             legend.key.height = unit(0.5, "line"), legend.key.width = unit(0.5, "line"))
#     
#   }))
# 
# plot_grid(plotlist = g_optimized_fitness_comparison_list$plot, ncol = 1)






# Supplemental Table 1: AIC comparison when using pedigree information --------

pedigree_model_comparison %>%
  rename(Unrelated = modelI_aic, Pedigree = modelA_aic, Trait = trait) %>%
  gather(model, aic, -nursery, -Trait) %>%
  arrange(nursery, desc(model)) %>%
  unite(model, nursery, model, sep = ".") %>%
  mutate(Trait = str_add_space(Trait), model = fct_inorder(model)) %>%
  spread(model, aic) %>%
  mutate_at(vars(-Trait), ~formatC(x = ., digits = 2, format = "f")) %>%
  write_csv(x = ., path = file.path(fig_dir, "tableS1_model_aic_comparison.csv"))




# Supplemental Table 2: number of mega-environments per trait -------------


## Summarize the number of ME per trait
me_summary_table <- nursery_mega_environment_output %>%
  unnest(me_data) %>%
  filter(type == "environment") %>%
  distinct(nursery, management, trait, mega_environment, nEnv) %>%
  group_by(nursery, management, trait) %>%
  summarize(nME = n_distinct(mega_environment), MinNe = min(nEnv), MaxNe = max(nEnv), MeanNe = mean(nEnv),
            TotalE = sum(nEnv)) %>%
  mutate(nMEenv = paste0(round(MeanNe, 2), " (", MinNe, ", ", MaxNe, ")")) %>%
  ungroup() %>%
  select(-contains("Ne")) %>%
  arrange(nursery, nME) %>%
  as.data.frame()

# Clean and save
me_summary_table %>%
  select(-management) %>%
  mutate(nursery = f_nursery_expand(nursery), trait = str_add_space(trait)) %>%
  rename_at(vars(nursery, trait), str_to_title) %>%
  write_csv(x = ., path = file.path(fig_dir, "tableS2_mega_environment_summary.csv"))




# Supplemental Table 3: location metadata -------------------------

trial_metadata %>% 
  filter(!(nursery == "wrn" & management == "rainfed")) %>% 
  mutate(abbreviation = str_sub(environment, 1, 3)) %>% 
  group_by(nursery, management, location) %>% 
  summarize(nYears = n_distinct(year), latitude = unique(lat), longitude = unique(long), abbreviation = unique(abbreviation)) %>% 
  ungroup() %>% 
  select(-management) %>% 
  spread(nursery, nYears) %>% 
  mutate_at(vars(mvn, wrn), ~ifelse(is.na(.), 0, .)) %>%
  select(location, abbreviation, names(.)) %>%
  rename_all(str_to_title) %>%
  rename(MVN = Mvn, WRN = Wrn) %>%
  write_csv(x = ., path = file.path(fig_dir, "tableS3_location_metadata.csv"))








# Other figures -----------------------------------------------------------






# Number of locations selected in each optimization procedure




# Plot the number of selected locations by penalty
optimized_locations_all_traits %>%
  filter(method == "optimization") %>%
  mutate(nSelectedLoc = map(optim_loc, ~as_tibble(map(.x, length)))) %>%
  unnest(nSelectedLoc) %>%
  mutate(agro_loc = agro_loc - maltq_loc) %>%
  gather(trait_type, nLoc, agro_loc, maltq_loc) %>%
  mutate(trait_type = ifelse(trait_type == "agro_loc", "justAgronomicTraits", "allTraits"),
         loc_penalty = fct_inseq(as.character(loc_penalty))) %>%
  ggplot(aes(x = loc_penalty, y = nLoc, fill = trait_type)) +
  geom_col() +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  facet_grid(~ nursery)








# Nursery genotype stats --------------------------------------------------

# Number of shared entries between nurseries
pheno_dat %>%
  left_join(., select(trial_metadata, trial, nursery, management)) %>%
  filter(!(nursery == "wrn" & management == "rainfed")) %>%
  distinct(line_name, nursery, year) %>%
  group_by(line_name) %>%
  summarize_at(vars(year, nursery), n_distinct) %>%
  arrange(desc(nursery), desc(year)) %>%
  as.data.frame()


# Row type per nursery
pheno_dat %>%
  left_join(., select(trial_metadata, trial, nursery, management)) %>%
  filter(!(nursery == "wrn" & management == "rainfed")) %>%
  distinct(line_name, nursery) %>%
  left_join(., line_metadata) %>%
  group_by(nursery, row_type) %>%
  summarize(n = n())

# Number of years of testing per entry
entry_obs <- pheno_dat %>%
  left_join(., select(trial_metadata, trial, nursery, management)) %>%
  filter(!(nursery == "wrn" & management == "rainfed")) %>%
  distinct(nursery, line_name, year) %>%
  group_by(nursery, line_name) %>%
  summarize(nYear = n_distinct(year)) 

summarize_at(entry_obs, vars(nYear), list(~min, ~max, ~median))

# How many tested in 3+ years?
filter(entry_obs, nYear >= 3) %>% 
  summarize(n = n()) %>%
  left_join(., summarize(entry_obs, nTotal = n_distinct(line_name))) %>%
  mutate(per = n / nTotal)

# Number of entry per trial
pheno_dat %>%
  left_join(., select(trial_metadata, trial, year1 = year, nursery, management)) %>%
  filter(!(nursery == "wrn" & management == "rainfed")) %>%
  group_by(nursery, trial) %>%
  summarize(nLine = n_distinct(line_name)) %>%
  summarize_at(vars(nLine), list(~min, ~max, ~median))


# Percent completeness (number of total trials/year) per trait
pheno_dat %>%
  left_join(., select(trial_metadata, trial, year1 = year, nursery, management)) %>%
  filter(!(nursery == "wrn" & management == "rainfed"))  %>% 
  distinct(nursery, trial, year, trait) %>% 
  group_by(nursery, year) %>% 
  mutate(nTotalTrial = n_distinct(trial)) %>% 
  group_by(nursery, trait, year) %>% 
  summarize(nTrialTrait = n_distinct(trial), nTotalTrial = mean(nTotalTrial)) %>% 
  mutate(percComplete = nTrialTrait / nTotalTrial) %>%
  summarize_at(vars(nTrialTrait, percComplete), list(~min, ~max, ~median)) %>%
  as.data.frame()



# What locations were selected by optimization?
optimized_locations_all_traits %>%
  filter(method != "random", loc_penalty == 0.01, trait_weight_group == "equal") %>%
  mutate(locations = map(optim_loc, "agro_loc")) %>%
  unnest(locations) %>%
  as.data.frame()

# How many years of data are available for the selected locations
selected_locations %>%
  mutate(nursery = tolower(nursery_abbr)) %>%
  select(nursery, location) %>%
  left_join(., aggregate(year ~ location + nursery + management, data = trial_metadata, FUN = n_distinct)) %>% 
  group_by(nursery, management) %>% 
  summarize_at(vars(year), list(min = min, max = max, mean = mean))

# Average number of years of data for any one location
aggregate(year ~ location + nursery + management, data = trial_metadata, FUN = n_distinct) %>%
  aggregate(year ~ nursery, data = ., mean)


# How many locations were tested in each of the past 5 years?
trial_metadata %>% 
  filter(year %in% tail(unique(.$year), 5)) %>% 
  group_by(nursery, management, year) %>% 
  summarize(nLoc = n_distinct(location))



## Why was Bottineau not selected in the MVN? Was it due to excessive GxE?
location_opimization_input %>%
  filter(nursery == "mvn") %>%
  mutate(avg_GL = map(GL, ~as.data.frame(.x) %>% rownames_to_column("location1") %>% 
                        gather(location2, GL, -location1) %>% filter(location1 != location2) %>% 
                        group_by(location1) %>% summarize(GL = mean(GL)) %>% mutate(rank_GL = rank(GL)) )) %>%
  unnest(avg_GL) %>% 
  group_by(location1) %>%
  summarize(rank_GL = mean(rank_GL))




## Compare best rank locations with optimization
optimization_bestrank_compare <- optimized_locations_all_traits %>%
  filter(loc_penalty == 0.01, method != "random", trait_weight_group == "equal", comp_weight_group == "equal",
         scale_components) %>%
  select(nursery, management, method, optim_loc) %>%
  mutate(optim_loc = map(optim_loc, unlist) %>% map(unique)) %>%
  spread(method, optim_loc) %>%
  mutate(intersect = map2(optimization, best_rank, intersect))

# Locations in both
optimization_bestrank_compare$intersect
  
optimized_locations <- optimized_locations_all_traits %>%
  filter(method == "optimization", trait_weight_group == "equal") %>%
  select(nursery, management, method, loc_penalty, optim_loc)
  
