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

# Load figure data
figure_results_objects <- load(file.path(result_dir, "figure_data.RData"))


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
  mutate(nursery = ifelse(nursery != "Both", f_nursery_expand(nursery), nursery),
         nursery = fct_relevel(nursery, "Both", after = Inf))


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
  scale_y_continuous(name = expression('Variance of a genotype mean ['*log[10](italic(V[bar(Y)]))*']'), breaks = pretty, 
                     labels = function(x) ifelse(x > 1e4, formatC(x, format = "e", digits = 1), x)) +
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
  
# Combine plots
g_precision_correlation <- plot_grid(g_env_varY_figure_traits, g_correlation_figure_traits, nrow = 1, 
                                     labels = letters, label_size = 8, align = "hv")
# Save
ggsave(filename = "figure1_reliability_repeatability_example.jpg", path = fig_dir, plot = g_precision_correlation,
       height = 8, width = 17.5, units = "cm", dpi = 1000)



# Figure 2: Location performance plots ------------------------------------

figure_traits1 <- c("GrainYield", "HeadingDate", "GrainProtein", "DiastaticPower", "MaltExtract")

location_opimization_input_toplot <- location_opimization_input %>%
  # Combine data
  mutate_at(vars(repeatability, representativeness, varY), ~map(., ~as.data.frame(.) %>% rownames_to_column("location"))) %>%
  mutate(data_toplot = pmap(select(., repeatability, representativeness, varY), list) %>% map(~reduce(., full_join, by = "location"))) %>%
  unnest(data_toplot) %>%
  mutate(representativeness = 1 - (repre / 90)) %>%
  # Add location abbreviations
  left_join(., location_abbr_key)


## Plot
location_performance_plots <- location_opimization_input_toplot %>%
  group_by(nursery, management, trait) %>%
  do(plot = {
    df <- .
    df1 <- df %>%
      mutate(representativeness1 = scale_01(repre, high.favorable = FALSE),
             representativeness = 1 - (repre / 90)) %>%
      # Add me
      left_join(., distinct_at(location_performance_toplot, vars(-precision, -repeatability, -representativeness, -maa)),
                by = c("trait", "nursery", "management", "location"))
    
    # Color scale breaks
    color_breaks <- pretty(df1$representativeness, min.n = 5, n = 5)
    names(color_breaks) <- c("", "Low", rep("", length(color_breaks) - 4),  "High", "")
    # names(color_breaks) <- c("Low", rep("", length(color_breaks) - 2), "High")
    
    ggplot(data = df1, aes(x = log10(varY), y = repeatability)) +
      geom_text_repel(aes(label = location_abbr), size = 3) +
      geom_point(aes(color = representativeness, shape = mega_environment), size = 2.5) +
      scale_x_reverse(name = expression('Precision ['*log[10](italic(V[bar(Y)]))*']'), breaks = pretty,
                      labels = function(x) ifelse(x > 1e4, formatC(x, format = "e", digits = 1), x)) +
      scale_y_continuous(name = expression("Repeatability (within-location "~rho[G]*")"), breaks = pretty) +
      scale_shape_discrete(name = "Mega environment (ME)", guide = guide_legend(order = 1)) +
      scale_color_gradient(low = "grey75", high = "blue", breaks = color_breaks, labels = names(color_breaks), 
                           name = "ME representativeness") +
      labs(subtitle = paste0(c(str_add_space(unique(df$trait)), toupper(unique(df$nursery)), str_to_title(unique(df$management))), 
                             collapse = ", ")) +
      theme_acs2(base_size = 8)
    
  }) %>% ungroup()



# Subset relevant traits
performance_toplot <- filter(location_performance_plots, trait %in% figure_traits1, 
                             !(nursery == "wrn" & management == "rainfed"))

# Color scale breaks
color_breaks <- performance_toplot$plot %>% 
  map("data") %>%
  map("representativeness") %>%
  unlist() %>%
  pretty(., min.n = 5, n = 5)
names(color_breaks) <- c("", "Low", rep("", length(color_breaks) - 4),  "High", "")


# Create merged plots by nursery
merged_plots_by_nursery <- performance_toplot %>%
  mutate(n = as.numeric(as.factor(nursery))) %>%
  group_by(nursery) %>%
  do({
    df <- .
    
    # Plot gy
    performance_partA <- subset(df, trait == "GrainYield", plot, drop = TRUE)[[1]] +
      labs(subtitle = NULL) +
      facet_grid(nursery ~ trait, labeller = labeller(nursery = f_nursery_expand, trait = str_add_space), 
                 switch = "y") +
      scale_shape_discrete(name = "Mega environment (ME)", labels = NULL, 
                           guide = guide_legend(order = 1, keywidth = unit(0.5, "pt"))) +
      scale_color_gradient(low = "grey90", high = "blue", breaks = color_breaks, name = "ME representativeness", 
                           guide = guide_colorbar(label.position = "top")) +
      theme(strip.placement = "outside", panel.border = element_rect(fill = alpha("white", 0)), plot.title = element_text(hjust = 0.5),
            legend.position  = "bottom", legend.background = element_rect(fill = alpha("white", 0), linetype = 0), 
            legend.box.background = element_rect(fill = alpha("white", 0), linetype = 0), legend.box.just = "bottom")
    
    # Adjust point and text size
    performance_partA$layers[[1]]$aes_params$size <- 2
    performance_partA$layers[[2]]$aes_params$size <- 2
    
    # Plot all other traits
    partB_list <- subset(df, trait != "GrainYield", plot, drop = TRUE)
    performance_partB_plot1 <- partB_list[[1]]
    performance_partB_plot1$data <- bind_rows(map(partB_list, "data")) %>%
      mutate(group = paste0(nursery, "_", trait))
    
    # Adjust text and point size
    performance_partB_plot1$layers[[2]]$aes_params$size <- 1.5
    performance_partB_plot1$layers[[1]]$aes_params$size <- 1.5
    
    performance_partB <- performance_partB_plot1 +
      labs(subtitle = NULL) +
      facet_wrap(~ group, ncol = 2, scales = "free", labeller = labeller(group = function(x) str_add_space(str_remove(x, "[a-z]*_")))) +
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
g_performance_merge2 <- plot_grid(legend, g_performance_merge1, x_axis, ncol = 1, rel_heights = c(0.09, 1, 0.05))


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
  unnest(fitness_scaled) %>%
  gather(fitness_component, value, varY, repAvg, reprAvg) %>%
  arrange(nursery, management, loc_penalty, trait_weight_group, method, fitness_component) %>% # Order
  mutate(loc_penalty = fct_inseq(as.character(loc_penalty)),
         fitness_component = fct_inorder(fitness_component) %>% fct_relevel("varY"),
         method = str_replace_all(method, "_", " ") %>% str_to_title() %>% str_remove_all(" ") %>%
           fct_inorder() %>% fct_relevel("Optimization"))

# Subset to plot
optimized_fitness_scaled1 <- optimized_fitness_scaled %>%
  filter(trait_weight_group == "equal") %>%
  # Merge method with number of optimal locations
  mutate(nLoc = ifelse(method == "Optimization", nOptimLoc, NA)) %>%
  group_by(nursery, management, loc_penalty, trait_weight_group) %>%
  mutate(nLoc = mean(nLoc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(loc_penalty1 = paste0(loc_penalty, "\n(n=", nLoc, ")"),
         method1 = paste0(method, " (n=", nOptimLoc, ")"))

## Add an annotation df
optimized_fitness_scaled_ann <- optimized_fitness_scaled1 %>%
  filter(loc_penalty == last(loc_penalty), method == "BestRank", fitness_component == "varY") %>%
  distinct_at(vars(nursery, management, contains("loc_penalty"), contains("method"), fitness_component, value))


g_optimized_fitness_comparison <- ggplot(data = NULL, aes(x = loc_penalty1, y = value)) +
  geom_hline(data = distinct(filter(optimized_fitness_scaled1, method == "BestRank"), nursery, method, fitness_component, value),
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
optimized_fitness_trait_scaled <- optimized_locations_all_traits %>%
  # Count the number of locations
  mutate(nOptimLoc = map(optim_loc, "agro_loc") %>% map_dbl(length),
         fitness_scaled_trait = map(fitness_scaled_trait, bind_rows)) %>%
  unnest(fitness_scaled_trait) %>%
  gather(fitness_component, value, varY, repAvg, reprAvg) %>%
  arrange(nursery, management, loc_penalty, trait_weight_group, method, fitness_component) %>% # Order
  mutate(fitness_component = fct_inorder(fitness_component) %>% fct_relevel("varY"),
         method = str_replace_all(method, "_", " ") %>% str_to_title() %>% str_remove_all(" ") %>%
           fct_inorder() %>% fct_relevel("Optimization"))


# Subset to plot
optimized_fitness_trait_scaled1 <- optimized_fitness_trait_scaled %>%
  filter(trait_weight_group == "equal") %>%
  filter(trait %in% figure_traits1) %>% # Subset traits
  # Merge method with number of optimal locations
  mutate(nLoc = ifelse(method == "Optimization", nOptimLoc, NA)) %>%
  group_by(nursery, management, loc_penalty, trait_weight_group) %>%
  mutate(nLoc = mean(nLoc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(loc_penalty1 = paste0(loc_penalty, " (n=", nLoc, ")"),
         method1 = paste0(method, " (n=", nOptimLoc, ")"))



## Only show the zero penalty and the selected penalty
select_loc_penalty <- 0.01

optimized_fitness_trait_scaled2 <- optimized_fitness_trait_scaled1 %>%
  filter(loc_penalty %in% c(0, select_loc_penalty))
  
# Create a list of colors
penalty_colors <- optimized_fitness_trait_scaled2 %>% distinct(loc_penalty, loc_penalty1) %>%
  mutate(i = as.numeric(as.factor(loc_penalty)), 
         color = map_chr(i, ~neyhart_palette("umn2")[-1:-2][.x])) %>%
  {setNames(object = .$color, nm = .$loc_penalty1)}

# Quantiles of random output
rand_quantiles <- c(0.10, 0.90)


## Split by nursery and create plots
optimized_fitness_comparison_trait_plotlist <- optimized_fitness_trait_scaled2 %>%
  group_by(nursery, management) %>%
  nest() %>%
  mutate(n = seq_len(nrow(.))) %>%
  group_by(nursery, management) %>%
  do(plot = {
    row <- .
    
    df <- row$data[[1]] %>%
      crossing(., select(row, nursery, management))
    n <- row$n
  
    # Create a custom grouping variable
    df1 <- df %>%
      mutate(method = as.character(method),
             method = ifelse(method == "BestRank", method1, method),
             method = fct_inorder(method)) %>%
      arrange(trait, method, loc_penalty) %>%
      mutate(group = ifelse(str_detect(method, "BestRank"), as.character(method), as.character(loc_penalty)),
             group = paste(trait, group, sep = "_"),
             group = fct_inorder(group))
    
    # Boxplots
    df1_boxplot <- filter(df1, method != "Optimization", !(method == "BestRank" & loc_penalty == last(loc_penalty)))  %>% 
      mutate(fill = ifelse(method == "Random", loc_penalty1, "bestRank"),
             fill = fct_inorder(fill))
    # Boxplot colors
    boxplot_colors <- ifelse(str_detect(distinct(arrange(df1_boxplot, fitness_component, trait, group), 
                                                 group, fitness_component, method)$method, "BestRank"), "white", "black")
    
    
    # Linerange
    df1_linerange <- filter(df1, method != "Optimization", !(method == "BestRank" & loc_penalty == last(loc_penalty)))  %>% 
      group_by(group, nursery, trait, loc_penalty1, method, fitness_component) %>%
      summarize_at(vars(value), list(lower = ~quantile(., rand_quantiles[1]), upper = ~quantile(., rand_quantiles[2]), 
                                     value = ~mean(.)), na.rm = TRUE) %>%
      ungroup() %>%
      mutate(fill = ifelse(method == "Random", loc_penalty1, "bestRank"),
             fill = fct_inorder(fill))
      
    # Points
    df1_points <- filter(df1, method != "Random") %>%
      mutate(color = ifelse(method == "Optimization", loc_penalty1, "bestRank"),
             color = fct_inorder(color))
    
    # Text annotation
    df1_text <- df1_points %>%
      filter(fitness_component == "varY", trait == "DiastaticPower", loc_penalty == first(loc_penalty)) %>%
      distinct(group, trait, value, method, fitness_component)
      
    df1_text1 <- df1_text %>%
      filter(method == "Optimization") %>%
      mutate(arrow_end_y = value + 0.02, arrow_end_x = 1.1, arrow_begin_y = 0.98, arrow_begin_x = 2)
    
    
    
    # Plot
    g <- ggplot(data = NULL, aes(x = trait, y = value, group = group)) +
      # Place an invisible point at 1
      geom_point(data = mutate(distinct(df1, nursery, trait, fitness_component, group), value = 1), 
                 shape = NA, key_glyph = "polygon", position = position_dodge(0.75)) +
      geom_boxplot(data = df1_boxplot, aes(fill = fill), position = position_dodge(0.75), alpha = 0.25, width = 0.5,
                   outlier.shape = NA, lwd = 0.25, color = boxplot_colors) +
      # geom_linerange(data = df1_linerange, aes(ymin = lower, ymax = upper, color = fill, alpha = "Random"), key_glyph = "rect", 
      #                lwd = 1.5, position = position_dodge(0.75)) +
      geom_point(data = df1_points, aes(color = color, shape = method), size = 1.5, position = position_dodge(0.75)) +
      facet_grid(fitness_component ~ nursery, scales = "free", switch = "y",
                 labeller = labeller(nursery = toupper, fitness_component = f_fitness_comp_replace)) +
      scale_y_continuous(name = NULL, breaks = pretty) +
      scale_x_discrete(name = NULL, labels = str_add_space, guide = guide_axis(check.overlap = TRUE, n.dodge = 2)) +
      scale_fill_manual(values = penalty_colors, guide = FALSE) +
      scale_color_manual(name = NULL, values = c("black", penalty_colors), breaks = names(penalty_colors), guide = guide_legend(ncol = 1)) +
      # scale_alpha_manual(name = NULL, values = 0.75, guide = guide_legend(order = 2, override.aes = list(fill = alpha("grey", 0.5)))) +
      scale_alpha_manual(guide = FALSE, values = 0.25) +
      scale_shape_discrete(name = NULL, labels = str_add_space, guide = guide_legend(order = 1, ncol = 1)) +
      # scale_shape_discrete(guide = FALSE, labels = str_add_space) +
      theme_acs2() +
      # theme(legend.background = element_rect(fill = alpha("white", 0)), legend.box.background = element_rect(fill = alpha("white", 0))) + 
      theme(strip.placement = "outside", legend.box = "horizontal", legend.spacing.y = unit(-0.1, "lines"), legend.direction = "horizontal",
            legend.position = "top", legend.justification = c(0, 0), legend.box.just = "right", legend.title = element_text(size = 8),
            legend.key.height = unit(0.5, "line"), legend.key.width = unit(0.5, "line"))
    
    if (n == 1) {
      # g <- g + 
      #   geom_text(data = df1_text1, aes(x = arrow_begin_x, y = arrow_begin_y, label = method), size = 2, hjust = 0, nudge_x = 0.05) +
      #   geom_curve(data = df1_text1, aes(x = arrow_begin_x, y = arrow_begin_y, xend = arrow_end_x, yend = arrow_end_y),
      #              arrow = arrow(angle = 20, ends = "last", length = unit(0.25, "line")), curvature = 0.25, size = 0.25)
        
    } else {
      g <- g + theme(strip.text.y = element_blank())
      
    }
    
    g
    
  }) %>% ungroup()

# Merge plots together
g_combine <- plot_grid(plotlist = optimized_fitness_comparison_trait_plotlist$plot, nrow = 1)



# Save
ggsave(filename = "figure3_optimized_fitness_comparison_partB.jpg", path = fig_dir, plot = g_combine,
       height = 9, width = 12, units = "cm", dpi = 1000)


## Merge plots
g_combine_all <- plot_grid(plotlist = c(list(g_optimized_fitness_comparison), optimized_fitness_comparison_trait_plotlist$plot), 
                           nrow = 1, align = "h", axis = "tblr", rel_widths = c(0.4, 0.31, 0.29), labels = letters[1:2], label_size = 10)


ggsave(filename = "figure3_optimized_fitness_comparison.jpg", path = fig_dir, plot = g_combine_all,
       height = 11, width = 17.8, units = "cm", dpi = 1000)





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
# 1 mvn           21         175    25       401    19
# 2 wrn           32         251    23       393    18
    

# Save
pheno_dat_summary %>%
  mutate(nursery = f_nursery_expand(nursery)) %>%
  rename(genotypes = line_name) %>%
  rename_all(str_to_title) %>%
  write_csv(x = ., path = file.path(fig_dir, "table1_data_summary.csv"))

## Data used in analysis
pheno_dat %>%
  left_join(., select(trial_metadata, trial, nursery, management)) %>%
  filter(!(nursery == "wrn" & management == "rainfed")) %>%
  group_by(trait, nursery, management) %>%
  summarize_at(vars(trial, location, year), n_distinct) %>%
  arrange(location) %>% as.data.frame()




    



# # Supplemental Figure 1: Proportion of variance explained by loadings --------
# 
# # Plot
# g_loadings_propVar <- loadings_prop_var %>%
#   mutate(trait_abb = abbreviate(str_add_space(trait), 2),
#          trait_abb = fct_reorder(trait_abb, loadings_var_prop, .fun = mean)) %>%
#   ggplot(aes(x = trait_abb, y = loadings_var_prop)) +
#   geom_boxplot() +
#   facet_grid(~ nursery, labeller = labeller(nursery = f_nursery_expand)) +
#   scale_x_discrete(name = NULL) +
#   scale_y_continuous(name = "Proportion of variance explained by loadings") +
#   theme_genetics(10) +
#   theme(panel.border = element_rect(fill = alpha("white", 0)), legend.position = "none")
# 
# # Save
# ggsave(filename = "figureS1_propvar_loadings.jpg", plot = g_loadings_propVar, path = fig_dir,
#        height = 3.5, width = 6, dpi = 1000)
# 
# 
# loadings_prop_var %>%
#   group_by(nursery, management, trait) %>%
#   summarize_at(vars(contains("var_prop")), mean) %>%
#   as.data.frame()




# Supplemental Figures 1: location precision (varR) -------------------

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
    tr <- unique(df1$trait)
    
    plot <- ggplot(df1, aes(x = x, y = log10(varY))) + 
      geom_jitter(width = 0.1, color = "grey85", size = 0.75) +
      geom_boxplot(alpha = 0, outlier.shape = NA, lwd = 0.25) +
      scale_y_continuous(name = expression('Variance of a genotype mean ['*log[10](italic(V[bar(Y)]))*']'), breaks = pretty, 
                         labels = function(x) ifelse(x > 1e4, formatC(x, format = "e", digits = 1), x)) +
      scale_x_discrete(name = NULL, labels = function(x) map_chr(str_split(x, "\\."), 3)) +
      labs(subtitle = paste0(c(f_trait_rename(tr), toupper(nur), str_to_title(mgmt)), collapse = ", ")) +
      theme_genetics(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), axis.text.y = element_text(size = 5), strip.placement = "outside")
    
    tibble(plot = list(plot), nRemoved = nrow(df) - nrow(df1))
    
  }) %>% ungroup()
    

## Save plots
precision_plot_list %>%
  split(.$nursery) %>%
  unname() %>%
  imap(~{
    
    nursery <- unique(.x$nursery)
    
    ## Modify plots
    plot_list1 <- map(.x$plot, ~. + theme(axis.title = element_blank()))
    
    # Add y axis grob
    y_axis <- grid::textGrob(label = expression('Variance of a genotype mean ['*log[10](italic(V[bar(Y)]))*']'), rot = 90)
    
    # Combine plots
    g_plot <- plot_grid(y_axis, add_sub(plot = plot_grid(plotlist = plot_list1, ncol = 4), label = "Location"), rel_widths = c(0.05, 1))
    ggsave(filename = paste0("figureS", 1+.y, "_precision_boxplot_by_location_", nursery, ".jpg"), plot = g_plot, path = fig_dir,
           height = 3 * ceiling(length(plot_list1) / 4), width = 12, dpi = 1000)
    
  })




# Supplemental Figures 4 & 5: location repeatability ----------------------


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
      geom_boxplot(alpha = 0, outlier.shape = NA) +
      scale_y_continuous(name = expression("Repeatability (within-location"~italic(rho[G])*")"), breaks = pretty) +
      scale_x_discrete(name = "Location") +
      labs(subtitle = paste0(c(f_trait_rename(tr), toupper(nur), str_to_title(mgmt)), collapse = ", ")) +
      theme_genetics(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside")
    
  }) %>% 
  ungroup()

## Save plots
location_cor_plot_list %>%
  split(.$nursery) %>%
  unname() %>%
  imap(~{
    
    nursery <- unique(.x$nursery)
    
    ## Modify plots
    plot_list1 <- map(.x$plot, ~. + theme(axis.title = element_blank()))
    
    # Add y axis grob
    y_axis <- grid::textGrob(label = expression("Repeatability (within-location"~italic(rho[G])*")"), rot = 90)
    
    # Combine plots
    g_plot <- plot_grid(y_axis, add_sub(plot = plot_grid(plotlist = plot_list1, ncol = 4), label = "Location"), rel_widths = c(0.05, 1))
    ggsave(filename = paste0("figureS", 3+.y, "_location_correlations_", nursery, ".jpg"), plot = g_plot, path = fig_dir,
           height = 3 * ceiling(length(plot_list1) / 4), width = 12, dpi = 1000)
    
  })




# Supplemental Figures 6 & 7: mega-environment assignment -----------------


## Save mega-environment plots
nursery_mega_environment_output_toplot <- nursery_mega_environment_output %>%
  mutate(delim = paste(nursery, management, sep = "_")) %>%
  split(.$delim)

## Iterate and save plots
for (i in seq_along(nursery_mega_environment_output_toplot)) {
  
  nursery_output <- nursery_mega_environment_output_toplot[[i]] %>%
    mutate(plot = map(plot, ~{
      .x$labels$subtitle <- f_trait_rename(.x$labels$subtitle)
      .x
    }))
  
  
  # Create a filename
  filename <- paste0("figureS", 5+i, "_mega_environments_", unique(nursery_output$nursery), "_",
                     unique(nursery_output$management), ".jpg")
  # Create a plot grid
  g_plot_grid <- plot_grid(plotlist = nursery_output$plot, ncol = 4)
  # Number of rows
  plot_rows <- ceiling(nrow(nursery_output) / 4)
  
  # Save
  ggsave(filename = filename, plot = g_plot_grid, path = fig_dir, width = 16, height = plot_rows * 4,
         dpi = 1000)
  
  
}



# Alternative supplemental figures

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
    
    # Adjust the subtitle
    g1$labels$subtitle <- f_trait_rename(g1$labels$subtitle)
    
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
for (i in seq_along(nursery_mega_environment_output_toplot)) {
  
  nursery_output <- nursery_mega_environment_output_toplot[[i]]
  
  # Create a filename
  filename <- paste0("figureS", 5+i, "_mega_environment_location_axes_", unique(nursery_output$nursery), "_",
                     unique(nursery_output$management), ".jpg")
  # Create a plot grid
  g_plot_grid <- plot_grid(plotlist = nursery_output$plot1, ncol = 4)
  # Number of rows
  plot_rows <- ceiling(nrow(nursery_output) / 4)
  
  # Save
  ggsave(filename = filename, plot = g_plot_grid, path = fig_dir, width = 16, height = plot_rows * 4,
         dpi = 1000)
  
  
}


location_performance_plot_data <- location_performance_plots %>% 
  mutate(plot_data = map(plot, "data")) %>% 
  unnest(plot_data) %>%
  select(nursery, management, trait, location, repre)

## Summarize angles for GP and GY
location_performance_plot_data %>%
  filter(trait %in% figure_traits) %>%
  arrange(nursery, management, trait, repre) %>%
  as.data.frame()
  



# Summarize the angles
location_performance_plots %>% 
  mutate(plot_data = map(plot, "data")) %>% 
  unnest(plot_data) %>%
  group_by(nursery, management, trait) %>%
  summarize_at(vars(repre), list(~min, ~max, ~mean)) %>%
  ungroup() %>%
  as.data.frame()





# Supplemental Figures 8 & 9: overall location performance ----------------

# Split
location_performance_plots1 <- location_performance_plots %>%
  mutate(delim = paste(nursery, management, sep = "_")) %>%
  split(.$delim)

## Iterate and save plots
for (i in seq_along(location_performance_plots1)) {
  
  df <- location_performance_plots1[[i]] %>%
    mutate(plot = map(plot, ~{
      .x$labels$subtitle <- f_trait_rename(.x$labels$subtitle)
      .x
    }))
  
  # Create a filename
  filename <- paste0("figureS", 7+i, "_location_overall_performance_", unique(df$nursery), "_",
                     unique(df$management), ".jpg")
  # Create a plot grid
  g_plot_grid <- plot_grid(plotlist = df$plot, ncol = 4)
  # Number of rows
  plot_rows <- ceiling(nrow(df) / 4)
  
  # Save
  ggsave(filename = filename, plot = g_plot_grid, path = fig_dir, width = 16, height = plot_rows * 3,
         dpi = 1000)
  
  
}





# Supplemental Figure 10: optimization results for all traits --------------

# Subset to plot
optimized_fitness_trait_scaled1 <- optimized_fitness_trait_scaled %>%
  filter(trait_weight_group == "equal") %>%
  # Merge method with number of optimal locations
  mutate(nLoc = ifelse(method == "Optimization", nOptimLoc, NA)) %>%
  group_by(nursery, management, loc_penalty, trait_weight_group) %>%
  mutate(nLoc = mean(nLoc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(loc_penalty1 = paste0(loc_penalty, " (n=", nLoc, ")"),
         method1 = paste0(method, " (n=", nOptimLoc, ")"))


# Create a list of colors
penalty_colors <- optimized_fitness_trait_scaled1 %>% distinct(loc_penalty, loc_penalty1) %>%
  mutate(i = as.numeric(as.factor(loc_penalty)), 
         color = map_chr(i, ~neyhart_palette("umn2")[-1:-2][.x])) %>%
  {setNames(object = .$color, nm = .$loc_penalty1)} %>%
  .[order(names(.))]

# Quantiles of random output
rand_quantiles <- c(0.10, 0.90)


## Split by nursery and create plots
optimized_fitness_comparison_trait_plotlist <- optimized_fitness_trait_scaled1 %>%
  group_by(nursery, management) %>%
  nest() %>%
  mutate(n = seq_len(nrow(.))) %>%
  group_by(nursery, management) %>%
  do(plot = {
    row <- .
    
    df <- row$data[[1]] %>%
      crossing(., select(row, nursery, management))
    n <- row$n
    
    # Create a custom grouping variable
    df1 <- df %>%
      mutate(method = as.character(method),
             method = ifelse(method == "BestRank", method1, method),
             method = fct_inorder(method)) %>%
      arrange(trait, method, loc_penalty) %>%
      mutate(group = ifelse(str_detect(method, "BestRank"), as.character(method), as.character(loc_penalty)),
             group = paste(trait, group, sep = "_"),
             group = fct_inorder(group))
    
    # Boxplots
    df1_boxplot <- filter(df1, method != "Optimization", !(method == "BestRank" & loc_penalty == last(loc_penalty)))  %>% 
      mutate(fill = ifelse(method == "Random", loc_penalty1, "bestRank"),
             fill = fct_inorder(fill))
    # Boxplot colors
    boxplot_colors <- ifelse(str_detect(distinct(arrange(df1_boxplot, fitness_component, trait, group), 
                                                 group, fitness_component, method)$method, "BestRank"), "white", "black")
    
    
    # Linerange
    df1_linerange <- filter(df1, method != "Optimization", !(method == "BestRank" & loc_penalty == last(loc_penalty)))  %>% 
      group_by(group, nursery, trait, loc_penalty1, method, fitness_component) %>%
      summarize_at(vars(value), list(lower = ~quantile(., rand_quantiles[1]), upper = ~quantile(., rand_quantiles[2]), 
                                     value = ~mean(.)), na.rm = TRUE) %>%
      ungroup() %>%
      mutate(fill = ifelse(method == "Random", loc_penalty1, "bestRank"),
             fill = fct_inorder(fill))
    
    # Points
    df1_points <- filter(df1, method != "Random") %>%
      mutate(color = ifelse(method == "Optimization", loc_penalty1, "bestRank"),
             color = fct_inorder(color))
    
    # Text annotation
    df1_text <- df1_points %>%
      filter(fitness_component == "varY", trait == "DiastaticPower", loc_penalty == first(loc_penalty)) %>%
      distinct(group, trait, value, method, fitness_component)
    
    df1_text1 <- df1_text %>%
      filter(method == "Optimization") %>%
      mutate(arrow_end_y = value + 0.02, arrow_end_x = 1.1, arrow_begin_y = 0.98, arrow_begin_x = 2)
    
    # Subset the colors
    penalty_colors_use <- penalty_colors[unique(df1$loc_penalty1)]
    
    
    # Plot
    ggplot(data = NULL, aes(x = trait, y = value, group = group)) +
      # Place an invisible point at 1
      geom_point(data = mutate(distinct(df1, nursery, trait, fitness_component, group), value = 1), 
                 shape = NA, key_glyph = "polygon", position = position_dodge(0.75)) +
      geom_boxplot(data = df1_boxplot, aes(fill = fill), position = position_dodge(0.75), alpha = 0.25, width = 0.5,
                   outlier.shape = NA, lwd = 0.25, color = boxplot_colors) +
      # geom_linerange(data = df1_linerange, aes(ymin = lower, ymax = upper, color = fill, alpha = "Random"), key_glyph = "rect", 
      #                lwd = 1.5, position = position_dodge(0.75)) +
      geom_point(data = df1_points, aes(color = color, shape = method), size = 1.5, position = position_dodge(0.75)) +
      facet_grid(fitness_component ~ nursery, scales = "free", switch = "y",
                 labeller = labeller(nursery = toupper, fitness_component = f_fitness_comp_replace)) +
      scale_y_continuous(name = NULL, breaks = pretty) +
      scale_x_discrete(name = NULL, labels = f_trait_rename, guide = guide_axis(check.overlap = TRUE, n.dodge = 2)) +
      scale_fill_manual(values = penalty_colors, guide = FALSE) +
      scale_color_manual(name = NULL, values = c("black", penalty_colors_use), breaks = names(penalty_colors_use), 
                         guide = guide_legend(nrow = 2)) +
      # scale_alpha_manual(name = NULL, values = 0.75, guide = guide_legend(order = 2, override.aes = list(fill = alpha("grey", 0.5)))) +
      scale_alpha_manual(guide = FALSE, values = 0.25) +
      scale_shape_discrete(name = NULL, labels = str_add_space, guide = guide_legend(order = 1, ncol = 1)) +
      # scale_shape_discrete(guide = FALSE, labels = str_add_space) +
      theme_acs2() +
      # theme(legend.background = element_rect(fill = alpha("white", 0)), legend.box.background = element_rect(fill = alpha("white", 0))) + 
      theme(strip.placement = "outside", legend.box = "horizontal", legend.spacing.y = unit(-0.1, "lines"), legend.direction = "horizontal",
            legend.position = "top", legend.justification = c(0, 0), legend.box.just = "right", legend.title = element_text(size = 8),
            legend.key.height = unit(0.5, "line"), legend.key.width = unit(0.5, "line"))

  }) %>% ungroup()

# Merge plots together
g_combine <- plot_grid(plotlist = optimized_fitness_comparison_trait_plotlist$plot, ncol = 1)


# Save
ggsave(filename = "figureS10_optimized_fitness_comparison_all_traits.jpg", path = fig_dir, plot = g_combine,
       height = 18, width = 18, units = "cm", dpi = 1000)





# Supplemental Figure 11: optimization results for all traits using trait weights ----

# Subset to plot
optimized_fitness_trait_scaled1 <- optimized_fitness_trait_scaled %>%
  filter(trait_weight_group == "weighted") %>%
  # Merge method with number of optimal locations
  mutate(nLoc = ifelse(method == "Optimization", nOptimLoc, NA)) %>%
  group_by(nursery, management, loc_penalty, trait_weight_group) %>%
  mutate(nLoc = mean(nLoc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(loc_penalty1 = paste0(loc_penalty, " (n=", nLoc, ")"),
         method1 = paste0(method, " (n=", nOptimLoc, ")"))

# Create a list of colors
penalty_colors <- optimized_fitness_trait_scaled1 %>% distinct(loc_penalty, loc_penalty1) %>%
  mutate(i = as.numeric(as.factor(loc_penalty)), 
         color = map_chr(i, ~neyhart_palette("umn2")[-1:-2][.x])) %>%
  {setNames(object = .$color, nm = .$loc_penalty1)} %>%
  .[order(names(.))]

# Quantiles of random output
rand_quantiles <- c(0.10, 0.90)


## Split by nursery and create plots
optimized_fitness_comparison_trait_plotlist <- optimized_fitness_trait_scaled1 %>%
  group_by(nursery, management) %>%
  nest() %>%
  mutate(n = seq_len(nrow(.))) %>%
  group_by(nursery, management) %>%
  do(plot = {
    row <- .
    
    df <- row$data[[1]] %>%
      crossing(., select(row, nursery, management))
    n <- row$n
    
    # Create a custom grouping variable
    df1 <- df %>%
      mutate(method = as.character(method),
             method = ifelse(method == "BestRank", method1, method),
             method = fct_inorder(method)) %>%
      arrange(trait, method, loc_penalty) %>%
      mutate(group = ifelse(str_detect(method, "BestRank"), as.character(method), as.character(loc_penalty)),
             group = paste(trait, group, sep = "_"),
             group = fct_inorder(group),
             trait = paste0(trait, " (", trait_weights$weighted[trait], ")"))
    
    # Boxplots
    df1_boxplot <- filter(df1, method != "Optimization", !(method == "BestRank" & loc_penalty == last(loc_penalty)))  %>% 
      mutate(fill = ifelse(method == "Random", loc_penalty1, "bestRank"),
             fill = fct_inorder(fill))
    # Boxplot colors
    boxplot_colors <- ifelse(str_detect(distinct(arrange(df1_boxplot, fitness_component, trait, group), 
                                                 group, fitness_component, method)$method, "BestRank"), "white", "black")
    
    
    # Linerange
    df1_linerange <- filter(df1, method != "Optimization", !(method == "BestRank" & loc_penalty == last(loc_penalty)))  %>% 
      group_by(group, nursery, trait, loc_penalty1, method, fitness_component) %>%
      summarize_at(vars(value), list(lower = ~quantile(., rand_quantiles[1]), upper = ~quantile(., rand_quantiles[2]), 
                                     value = ~mean(.)), na.rm = TRUE) %>%
      ungroup() %>%
      mutate(fill = ifelse(method == "Random", loc_penalty1, "bestRank"),
             fill = fct_inorder(fill))
    
    # Points
    df1_points <- filter(df1, method != "Random") %>%
      mutate(color = ifelse(method == "Optimization", loc_penalty1, "bestRank"),
             color = fct_inorder(color))
    
    # Text annotation
    df1_text <- df1_points %>%
      filter(fitness_component == "varY", trait == "DiastaticPower", loc_penalty == first(loc_penalty)) %>%
      distinct(group, trait, value, method, fitness_component)
    
    df1_text1 <- df1_text %>%
      filter(method == "Optimization") %>%
      mutate(arrow_end_y = value + 0.02, arrow_end_x = 1.1, arrow_begin_y = 0.98, arrow_begin_x = 2)
    
    # Subset the colors
    penalty_colors_use <- penalty_colors[unique(df1$loc_penalty1)]
    
    
    # Plot
    ggplot(data = NULL, aes(x = trait, y = value, group = group)) +
      # Place an invisible point at 1
      geom_point(data = mutate(distinct(df1, nursery, trait, fitness_component, group), value = 1), 
                 shape = NA, key_glyph = "polygon", position = position_dodge(0.75)) +
      geom_boxplot(data = df1_boxplot, aes(fill = fill), position = position_dodge(0.75), alpha = 0.25, width = 0.5,
                   outlier.shape = NA, lwd = 0.25, color = boxplot_colors) +
      # geom_linerange(data = df1_linerange, aes(ymin = lower, ymax = upper, color = fill, alpha = "Random"), key_glyph = "rect", 
      #                lwd = 1.5, position = position_dodge(0.75)) +
      geom_point(data = df1_points, aes(color = color, shape = method), size = 1.5, position = position_dodge(0.75)) +
      facet_grid(fitness_component ~ nursery, scales = "free", switch = "y",
                 labeller = labeller(nursery = toupper, fitness_component = f_fitness_comp_replace)) +
      scale_y_continuous(name = NULL, breaks = pretty) +
      scale_x_discrete(name = NULL, labels = f_trait_rename, guide = guide_axis(check.overlap = TRUE, n.dodge = 2)) +
      scale_fill_manual(values = penalty_colors, guide = FALSE) +
      scale_color_manual(name = NULL, values = c("black", penalty_colors_use), breaks = names(penalty_colors_use), 
                         guide = guide_legend(nrow = 2)) +
      # scale_alpha_manual(name = NULL, values = 0.75, guide = guide_legend(order = 2, override.aes = list(fill = alpha("grey", 0.5)))) +
      scale_alpha_manual(guide = FALSE, values = 0.25) +
      scale_shape_discrete(name = NULL, labels = str_add_space, guide = guide_legend(order = 1, ncol = 1)) +
      # scale_shape_discrete(guide = FALSE, labels = str_add_space) +
      theme_acs2() +
      # theme(legend.background = element_rect(fill = alpha("white", 0)), legend.box.background = element_rect(fill = alpha("white", 0))) + 
      theme(strip.placement = "outside", legend.box = "horizontal", legend.spacing.y = unit(-0.1, "lines"), legend.direction = "horizontal",
            legend.position = "top", legend.justification = c(0, 0), legend.box.just = "right", legend.title = element_text(size = 8),
            legend.key.height = unit(0.5, "line"), legend.key.width = unit(0.5, "line"))
    
  }) %>% ungroup()

# Merge plots together
g_combine <- plot_grid(plotlist = optimized_fitness_comparison_trait_plotlist$plot, ncol = 1)


# Save
ggsave(filename = "figureS11_optimized_fitness_comparison_all_traits_weighted.jpg", path = fig_dir, plot = g_combine,
       height = 18, width = 18, units = "cm", dpi = 1000)





# Supplemental Table 1: number of environments, year, locations per trait --------

# Summarize data available for each nursery (exclude rainfed WRN)
pheno_dat_summary_trait <- pheno_dat %>%
  left_join(., select(trial_metadata, trial, nursery, management)) %>%
  filter(!(nursery == "wrn" & management == "rainfed")) %>%
  group_by(nursery, trait) %>%
  summarize_at(vars(location, environment, year, line_name), n_distinct) %>%
  ungroup()

# nursery location environment  year line_name trait
# 1 mvn           21         175    25       401    19
# 2 wrn           32         251    23       393    18


# Save
pheno_dat_summary_trait %>%
  filter(trait %in% traits_keep) %>%
  mutate(nursery = f_nursery_expand(nursery),
         trait = str_add_space(trait)) %>%
  arrange(nursery, desc(environment), trait) %>%
  rename(genotypes = line_name) %>%
  rename_all(str_to_title) %>%
  write_csv(x = ., path = file.path(fig_dir, "tableS1_data_summary_trait.csv"))  






# Supplemental Table 2: AIC comparison when using pedigree information --------

pedigree_model_comparison %>%
  rename(Unrelated = modelI_aic, Pedigree = modelA_aic, Trait = trait) %>%
  gather(model, aic, -nursery, -Trait) %>%
  arrange(nursery, desc(model)) %>%
  unite(model, nursery, model, sep = ".") %>%
  mutate(Trait = str_add_space(Trait), model = fct_inorder(model)) %>%
  spread(model, aic) %>%
  mutate_at(vars(-Trait), ~formatC(x = ., digits = 2, format = "f")) %>%
  write_csv(x = ., path = file.path(fig_dir, "tableS2_model_aic_comparison.csv"))




# Supplemental Table 3: number of mega-environments per trait -------------


## Summarize the number of ME per trait
me_summary_table <- nursery_mega_environment_output %>%
  unnest(me_data) %>%
  filter(type == "environment") %>%
  distinct(nursery, management, trait, mega_environment, nEnv) %>%
  group_by(nursery, management, trait) %>%
  summarize(nME = n_distinct(mega_environment), MinNe = min(nEnv), MaxNe = max(nEnv), MeanNe = mean(nEnv),
            TotalE = sum(nEnv)) %>%
  ungroup() %>%
  arrange(nursery, nME) %>%
  as.data.frame()

# Clean and save
me_summary_table %>%
  select(-management) %>%
  mutate(nursery = f_nursery_expand(nursery), trait = str_add_space(trait),
         MeanNe = format_numbers(MeanNe)) %>%
  rename_at(vars(nursery, trait), str_to_title) %>%
  write_csv(x = ., path = file.path(fig_dir, "tableS3_mega_environment_summary.csv"))




# Supplemental Table 4: location metadata -------------------------

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
  write_csv(x = ., path = file.path(fig_dir, "tableS4_location_metadata.csv"))








# Other figures -----------------------------------------------------------






## Map the selected locations for the penalty of 0.01
selected_locations <- optimized_locations_all_traits %>%
  filter(method == "optimization", env_penalty == 0.01) %>% 
  mutate(optim_loc = map(optim_loc, ~reduce(.x, union))) %>% 
  unnest(optim_loc)


# Map
g_map1 <- ggplot(data = north_america, aes(x = long, y = lat)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = canada, aes(group = group), fill = "grey95", color = "grey50", lwd = 0.75) + # Add canada
  geom_polygon(data = usa_state, aes(group = group), fill = "grey95", color = "grey50", lwd = 0.75) +
  geom_point(data = trial_info_toplot1, aes(color = nursery), size = 2.5) +
  geom_point(data = inner_join(trial_info_toplot1, select(selected_locations, location = optim_loc)), 
             aes(color = nursery, shape = "Optimized"), size = 3.5) +
  geom_label_repel(data = inner_join(trial_info_toplot1, select(selected_locations, location = optim_loc)), 
                   aes(color = nursery, label = location), fill = alpha("white", 0.5)) +
  scale_color_viridis_d(begin = 0.2, end = 0.8, name = NULL) +
  scale_shape_manual(name = NULL, values = 15) +
  coord_map(projection = "bonne", lat0 = mean(ylim), xlim = xlim, ylim = ylim) +
  theme_void(base_size = 14) +
  theme(legend.position = c(0.20, 0.90), legend.background = element_rect(fill = "white", color = "white"),
        legend.margin = margin(), legend.spacing = unit(0, "line"))



# Save the figure
ggsave(filename = "trial_location_map_optimized.jpg", plot = g_map1, path = fig_dir,
       width = 8, height = 5, dpi = 1000)




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
optimized_locations_all_traits %>%
  filter(method == "optimization", loc_penalty == 0.01, trait_weight_group == "equal") %>%
  mutate(location = map(optim_loc, "agro_loc")) %>%
  unnest(location) %>%
  as.data.frame() %>%
  left_join(., aggregate(year ~ location + nursery + management, data = trial_metadata, FUN = n_distinct)) %>% 
  group_by(nursery, management) %>% 
  summarize_at(vars(year), list(~min, ~max, ~mean))

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



