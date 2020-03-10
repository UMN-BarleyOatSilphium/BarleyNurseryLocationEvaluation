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

# Color pallete for environments and line names
term_colors <- set_names(wes_palette("Zissou1")[c(1,5)], "line_name", "environment")

# Function to rename terms
f_term_replace <- function(x) str_replace_all(x, c("line_name" = "Genotype", "environment" = "Environment"))


#######################
# Map of locations
#######################



## Plot locations

# Get the map data for canada
canada <- map_data("world", "Canada")

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
  group_by(nursery, location, lat, long) %>% 
  summarize(nYear = n_distinct(year))

# Create a box for x and y axes
xlim <- range(pretty(trial_info_toplot$long))
ylim <- range(pretty(trial_info_toplot$lat))

# Map
g_map <- ggplot(data = north_america, aes(x = long, y = lat)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = canada, aes(group = group), fill = "grey85", color = "grey50", lwd = 0.75) + # Add canada
  geom_polygon(data = usa_state, aes(group = group), fill = "grey85", color = "grey50", lwd = 0.75) +
  geom_point(data = trial_info_toplot, aes(color = nursery), size = 2.5) +
  scale_color_viridis_d(begin = 0.2, end = 0.8, name = "Nursery", labels = toupper) +
  coord_map(projection = "mercator", xlim = xlim, ylim = ylim) +
  theme_void(base_size = 14) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_rect(fill = "white"),
        legend.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "line"))


# Save the figure
ggsave(filename = "trial_location_map.jpg", plot = g_map, path = fig_dir,
       width = 8, height = 5, dpi = 1000)










#######################
# Biplot
#######################

## Plot the fitted AMMI or GGE models


# Load data
load(file.path(result_dir, "bilinear_model_fit.RData"))


## Common plot modifyer
g_plot_mod <- list(
  scale_color_manual(values = term_colors, labels = f_term_replace, name = NULL),
  scale_x_continuous(breaks = pretty),
  scale_y_continuous(breaks = pretty),
  theme_light()
)


## Plotting df
bilinear_toplot <- bilinear_fit %>%
  group_by(trait, model) %>%
  do(plotting_df = {
    row <- .
    
    ## Choose plot based on ammi versus gge
    if (unique(row$model) == "ammi") {
      
      ## Get the IPCA scores
      row$fit[[1]]$scores %>%
        map(as.data.frame) %>%
        imap_dfr(~rownames_to_column(.x, "term") %>% mutate(group = ifelse(.y == "Gscores", "line_name", "environment"))) %>%
        # Add effects
        left_join(., map_df(row$fit[[1]][c("Geffect", "Eeffect")], ~tibble(term = names(.), effect = .)))
      
    } else {
      
      ## Get the IPCA scores; divide by singular values
      row$fit[[1]]$scores %>%
        map(~.x / matrix(data = sqrt(row$fit[[1]]$svdE$d[seq_len(ncol(.x))]), nrow = nrow(.x), ncol = ncol(.x))) %>%
        map(as.data.frame) %>%
        imap_dfr(~rownames_to_column(.x, "term") %>% mutate(group = ifelse(.y == "Gscores", "line_name", "environment"))) 

    }
    
  }) %>% ungroup()
    


# Group by trait and model; plot
biplots <- bilinear_toplot %>%
  group_by(trait, model) %>%
  do(plot = {
    row <- .
    df <- row$plotting_df[[1]]
    
    ## Choose plot based on ammi versus gge
    if (row$model == "ammi") {
      
      g_plot <- df %>%
        ggplot(aes(x = effect, y = PC1, color = group)) +
        geom_point(data = subset(df, group == "line_name"), size = 0.5) +
        geom_point(data = subset(df, group == "environment"), size = 1.5) +
        xlab("Effect") + ylab("IPCA1") +
        labs(subtitle = paste0("AMMI biplot, ", str_add_space(row$trait))) +
        g_plot_mod
      
    } else {
      
      g_plot <- df %>%
        ggplot(aes(x = PC1, y = PC2, color = group)) +
        geom_point(data = subset(df, group == "line_name"), size = 0.5) +
        geom_point(data = subset(df, group == "environment"), size = 1.5) +
        xlab("PC1") + ylab("PC2") +
        labs(subtitle = paste0("GGE biplot, ", str_add_space(row$trait))) +
        g_plot_mod
      
    }
    
    g_plot
    
  })

    
    












