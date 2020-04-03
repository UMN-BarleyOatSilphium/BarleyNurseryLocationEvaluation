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
library(GGEBiplots)
library(paletteer)
library(ggrepel)
library(cowplot)



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


## Prediction accuracy
pheno_blup %>%
  group_by(trait, nursery, management) %>%
  summarize(acc = cor(ppv, pheno, use = "pairwise.complete.obs"))


#######################
# Visualize phenotypic balance
#######################


## Visualize balance of lines in trials
ge_tab_pred <- pheno_blup %>% 
  split(.$nursery) %>%
  # Complete trial-line_name cases
  map_df(~complete(droplevels(.), nursery, line_name, trial, trait, fill = list(pheno = NA, ppv = NA))) %>%
  mutate(obs = case_when(
    !is.na(pheno) ~ "phenotype",
    !is.na(ppv) ~ "prediction",
    TRUE ~ as.character(NA)
  )) %>%
  group_by(line_name, trait, nursery) %>%
  mutate(nPheno = sum(obs == "phenotype", na.rm = T) / n_distinct(trial),
         nMissing = sum(is.na(obs)) / n_distinct(trial)) %>%
  ungroup() %>% 
  arrange(trait, trial, desc(nPheno), nMissing) %>%
  mutate_at(vars(line_name, trial), fct_inorder)



freq_df_plot <- ge_tab_pred %>%
  # Convert line_name and location to numeric
  mutate_at(vars(line_name, trial), list(num = ~fct_inseq(as.factor(as.numeric(.)))))


g_pred_freq <- freq_df_plot %>%
  filter(!is.na(obs)) %>%
  # filter(nursery == "mvn", trait == "GrainProtein") %>%
  ggplot(aes(x = trial_num, y = line_name_num, fill = obs)) +
  geom_tile() +
  scale_y_discrete(name = "Genotype") +
  scale_x_discrete(name = "Trial") +
  scale_fill_viridis_d(na.value = "white", begin = 0.2, end = 0.8,
                       guide = guide_legend(override.aes = list(color = "black"))) +
  facet_wrap(~ nursery + trait, scales = "free") + 
  theme_presentation2(base_size = 14) +
  theme(axis.text = element_blank(), legend.position = "bottom")

# Save
ggsave(filename = "genotype_trial_contingency_example.jpg", path = fig_dir, plot = g_pred_freq,
       height = 8, width = 5, dpi = 500)



## Table of phenotypic observations, predictions, and still missings
freq_df_plot %>%
  group_by(trait, nursery, obs) %>% 
  summarize(n = n()) %>% 
  spread(obs, n)


# Remove line names with any missing data
pheno_blup_tomodel <- pheno_blup %>%
  group_by(trait, nursery) %>%
  do({
    df <- .
    
    select(df, trial, line_name, ppv) %>%
      droplevels() %>%
      complete(trial, line_name, fill = list(ppv = NA)) %>%
      group_by(line_name) %>%
      filter(all(!is.na(ppv))) %>%
      # Take average of duplicate line_name-trial combinations
      group_by(line_name, trial) %>%
      summarize(ppv = mean(ppv)) %>%
      ungroup()
    
  }) %>% ungroup()


## Summarize the number of genotype-environment combinations
pheno_blup_tomodel %>%
  group_by(trait, nursery) %>% 
  summarize_at(vars(line_name, trial), n_distinct)



  

# # Re-plot contingencies
# g_pred_freq1 <- pheno_blup_tomodel %>%
#   # Convert line_name and location to numeric
#   mutate_at(vars(line_name, trial), as.factor) %>%
#   mutate_at(vars(line_name, trial), list(num = ~fct_inseq(as.factor(as.numeric(.))))) %>%
#   ggplot(aes(x = trial_num, y = line_name_num, fill = obs)) +
#   geom_tile() +
#   scale_y_discrete(name = "Genotype") +
#   scale_x_discrete(name = "Trial") +
#   scale_fill_viridis_d(na.value = "white", begin = 0.2, end = 0.8,
#                        guide = guide_legend(override.aes = list(color = "black"))) +
#   facet_wrap(~ nursery + trait, scales = "free") + 
#   theme_presentation2(base_size = 14) +
#   theme(axis.text = element_blank(), legend.position = "bottom")
# 
# # Save
# ggsave(filename = "genotype_trial_contingency_example.jpg", path = fig_dir, plot = g_pred_freq,
#        height = 8, width = 5, dpi = 500)




#######################
# Summarize error
#######################

## Summarize error variance by locations
met_var_comp <- met_mm_out %>% 
  unnest(var_comp) %>%
  filter(str_detect(component, "units")) %>%
  rename(environment = component) %>%
  mutate(environment = str_remove_all(environment, ":units"))

# Bind with trial metadata
met_var_comp1 <- met_var_comp %>%
  left_join(trial_metadata) %>%
  filter(variance != 0) %>%
  mutate(variance = ifelse(variance <= 0, NA, variance),
         sdev = sqrt(variance))

## Number of observations per location/trait
met_var_comp_plot <- met_var_comp1 %>% 
  group_by(trait, location) %>% 
  summarize(nE = n_distinct(environment)) %>%
  left_join(met_var_comp1, .)

# Plot error variance by location and trait
g_met_varR <- met_var_comp1 %>%
  ggplot(aes(location, y = sdev, fill = trait, text = environment)) +
  geom_jitter(aes(color = trait), width = 0.25) +
  # geom_boxplot(width = 0.5, alpha = 0.5) +
  scale_y_continuous(name = expression("Residual std. dev. ("~sigma[R]~")"), breaks = pretty, trans = "identity") +
  facet_grid(trait ~ nursery, scales = "free", switch = "y") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Plotly
plotly::ggplotly(g_met_varR)

# Save
ggsave(filename = "nursery_example_location_varR.jpg", path = fig_dir, plot = g_met_varR,
       height = 6, width = 8, dpi = 500)



## Plot residual SD across environments between traits
## Plot as pairs
resid_sd_plot_list <- met_var_comp1 %>%
  group_by(nursery, management) %>%
  do(pairs_plot = {
    dat <- .
    traits <- unique(dat$trait)
    
    g_plot <- dat %>% 
      select(environment, location, year, trait, sdev) %>%
      spread(trait, sdev) %>%
      ggpairs(aes(color = location), columns = traits,
              diag = list(continuous = wrap("densityDiag", alpha = 0.5))) +
      labs(subtitle = paste(toupper(unique(dat$nursery)), str_to_title(unique(dat$management)), sep = ", ")) +
      theme_presentation2(base_size = 10)
    
    # Return the plot
    g_plot
    
  }) %>% ungroup()

## Save plots
for (i in seq_len(nrow(resid_sd_plot_list))) {
  filename <- paste0(paste0(unlist(resid_sd_plot_list[i,1:2]), collapse = "_"), "_resid_var_trait_pairs.jpg")
  ggsave(filename = filename, plot = resid_sd_plot_list$pairs_plot[[i]], path = fig_dir,
         height = 5, width = 7, dpi = 1000)
  
  # Convert to plotly
  g_plotly <- plotly::ggplotly(p = resid_sd_plot_list$pairs_plot[[i]])
  # Save as HTML widget
  htmlwidgets::saveWidget(widget = plotly::as_widget(g_plotly),
                          file = file.path(fig_dir, gsub(pattern = "jpg", replacement = "html", x = filename)))
  
}











#######################
# Mega-environment analysis
#######################


## Fit a GGE model per trait ##
bilinear_fit <- pheno_blup_tomodel %>%
  group_by(trait, nursery) %>%
  do({
    df <- .
    df1 <- droplevels(df)
    
    # Fit the model
    ammi <- bilinear(x = df1, G = "line_name", E = "trial", y = "ppv", test = "bootstrap", B = 1, 
                    model = "AMMI")

    # Fit the GGE model
    gge <- bilinear(x = df1, G = "line_name", E = "trial", y = "ppv", test = "bootstrap", B = 1, 
                    model = "GGE")
    
    # return both
    tibble(model = c("ammi", "gge"), fit = list(ammi, gge))
    
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
  ggplot(aes(x = PC_num, y = propvar_SS, shape = model)) +
  geom_point() + 
  facet_wrap(~ nursery + trait, scales = "free_x")

## Determine significant PCs using "elbow" method
# Tolerance for difference in variance explained
tol <- 0.03

bilinear_sig_PCs <- bilinear_scores %>%
  unnest(varprop) %>%
  #
  mutate(propvar = propvar_SS) %>%
  #
  arrange(trait, model, PC_num) %>%
  group_by(trait, model, nursery) %>%
  do({
    mutate(., propvar_diff = c(abs(diff(propvar)), 0), stop = which.min(propvar_diff >= tol), 
           nPC = stop - 1)
  }) %>% ungroup()

## Summary df of number of sig PCs
(bilinear_sig_PCs_summ <- bilinear_sig_PCs %>%
  group_by(trait, model, nursery) %>%
  filter(PC_num %in% seq(1, unique(nPC))) %>%
  summarize(total_propvar = sum(propvar), nPC = unique(nPC)))


## Fit a model for that number of PCS
bilinear_fitN_fit <- bilinear_fit %>%
  left_join(., bilinear_sig_PCs_summ) %>%
  group_by(trait, model, nursery) %>%
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
    
    ## Get the singular value decomposition
    svd <- fitted$svdE
    
    g_scores <- svd$u %*% diag(sqrt(svd$d))
    dimnames(g_scores) <- list(g_effects$line_name, paste0("PC", seq_len(ncol(g_scores))))
    
    e_scores <- svd$v %*% diag(sqrt(svd$d))
    dimnames(e_scores) <- list(e_effects$environment, paste0("PC", seq_len(ncol(e_scores))))
    
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
    
  }) %>% ungroup() %>%
  left_join(., nest(group_by(select(bilinear_sig_PCs, trait, nursery, model, PC, propvar_SS), trait, nursery, model), .key = "PC_summ"))



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
  group_by(trait, model, nursery) %>%
  do({
    .x <- .
    
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
  group_by(trait, model, nursery) %>%
  do({
    group_by(., environment) %>% 
      filter(rank == min(rank)) %>% 
      ungroup() %>%
      mutate(me = as.numeric(as.factor(line_name)))
  }) %>%
  ungroup() %>%
  select(-mu, -y_hat, -dir, -rank)


bilinear_me %>%
  distinct(trait, model, nursery, me)



## Save the model output
save("bilinear_fit", "bilinear_fitN_fit", "bilinear_me", file = file.path(result_dir, "bilinear_model_fit.RData"))





#######################
# Model phenotypes
#######################


## Calculate correlations based on phenotypic data
location_gen_cor_pheno <- inner_join(pheno_dat, distinct(pheno_blup_tomodel, nursery, trait, trial, line_name)) %>%
  group_by(trait, nursery, line_name, location, environment) %>%
  summarize_at(vars(value), mean) %>%
  group_by(trait, nursery, location) %>%
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
  group_by(trait, location, nursery) %>%
  summarize_at(vars(corr, nCorr), mean)


# Plot
g_location_cor <- location_adj_year_cor %>%
  ggplot(aes(x = location, y = corr)) +
  geom_jitter(width = 0.25) +
  geom_boxplot(alpha = 0.5, fill = "grey85", width = 0.5) +
  facet_grid(trait ~ nursery, switch = "y", scales = "free_x", space = "free_x", drop = TRUE) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
# Save
ggsave(filename = "location_correlation_example.jpg", plot = g_location_cor, path = fig_dir,
       height = 5, width = 5, dpi = 1000)












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
  group_by(trait, model, nursery) %>%
  do(plotting_df = {
    row <- .
    
    ## Choose plot based on ammi versus gge
    if (unique(row$model) == "ammi") {
      
      ## Get the IPCA scores
      row$fit[[1]]$scores %>%
        map(as.data.frame) %>%
        imap_dfr(~rownames_to_column(.x, "term") %>% mutate(group = ifelse(.y == "Gscores", "line_name", "environment"))) %>%
        # Add effects
        left_join(., map_df(row$fit[[1]][c("Geffect", "Eeffect")], ~tibble(term = names(.), effect = .))) %>%
        left_join(., distinct(trial_metadata, trial, location, environment), by = c("term" = "trial"))
      
    } else {
      
      ## Get the IPCA scores; divide by singular values
      row$fit[[1]]$scores %>%
        map(~.x / matrix(data = sqrt(row$fit[[1]]$svdE$d[seq_len(ncol(.x))]), nrow = nrow(.x), ncol = ncol(.x))) %>%
        map(as.data.frame) %>%
        imap_dfr(~rownames_to_column(.x, "term") %>% mutate(group = ifelse(.y == "Gscores", "line_name", "environment"))) %>%
        left_join(., distinct(trial_metadata, trial, location, environment), by = c("term" = "trial"))
      
    }
    
  }) %>% ungroup() %>%
  # Add the PC summary
  left_join(., select(bilinear_fitN_fit, trait, model, nursery, PC_summ))



# Group by row plot
biplots <- bilinear_toplot %>%
  group_by(trait, model, nursery) %>%
  do(plot = {
    row <- .
    df <- row$plotting_df[[1]]
    PC_summ <- row$PC_summ[[1]]
    
    ## Choose plot based on ammi versus gge
    if (row$model == "ammi") {
      
      g_plot <- df %>%
        ggplot(aes(x = effect, y = PC1, color = group)) +
        geom_point(data = subset(df, group == "line_name"), size = 0.5, shape = 3) +
        geom_text(data = subset(df, group == "environment"), aes(label = environment), size = 2) +
        xlab("Effect") + 
        ylab(paste0("IPCA1 (", round(subset(PC_summ, PC == "PC1", propvar_SS, drop = TRUE), 2) * 100, "%)")) + 
        labs(subtitle = paste0("AMMI biplot, ", str_add_space(row$trait))) +
        g_plot_mod
      
    } else {
      
      g_plot <- df %>%
        ggplot(aes(x = PC1, y = PC2, color = group)) +
        geom_point(data = subset(df, group == "line_name"), size = 0.5, shape = 3) +
        geom_text(data = subset(df, group == "environment"), aes(label = environment), size = 2) +
        # geom_text_repel(data = subset(df, group == "environment"), aes(label = environment)) +
        # geom_point(data = subset(df, group == "environment"), size = 1.5) +
        xlab(paste0("PCA1 (", round(subset(PC_summ, PC == "PC1", propvar_SS, drop = TRUE), 2) * 100, "%)")) + 
        ylab(paste0("PCA2 (", round(subset(PC_summ, PC == "PC2", propvar_SS, drop = TRUE), 2) * 100, "%)")) + 
        labs(subtitle = paste0("GGE biplot, ", str_add_space(row$trait))) +
        g_plot_mod
      
    }
    
    g_plot
    
  })

## colors for the lines
segments_colors <- paletteer_d("ggsci", "uniform_startrek")


## Plot GGE-GGL biplots
# Group by row plot
ggl_biplots <- bilinear_toplot %>%
  filter(model == "gge") %>%
  group_by(trait, model, nursery) %>%
  do(plot = {
    row <- .
    df <- row$plotting_df[[1]]
    PC_summ <- row$PC_summ[[1]]
    
    
    ## Calculate average PCs for location
    location_df <- aggregate(cbind(PC1, PC2) ~ location + group, data = df, FUN = mean, subset = group == "environment") %>%
      rename_at(vars(PC1, PC2), ~paste0("location_", .)) %>%
      mutate(x0 = 0, x1 = location_PC1, y0 = 0, y1 = location_PC2)
    # Add location information to df
    df1 <- subset(df, group == "environment") %>% 
      left_join(., location_df) %>%
      mutate(x0 = location_PC1, x1 = PC1, y0 = location_PC2, y1 = PC2)
    
    ## GGL-GGE biplot
    g_plot <- df %>%
      ggplot(aes(x = PC1, y = PC2, color = group)) +
      geom_point(data = subset(df, group == "line_name"), size = 0.5, shape = 3) +
      geom_point(data = df1, color = segments_colors[1], size = 0.5) +
      geom_segment(data = location_df, aes(x = x0, y = y0, xend = x1, yend = y1), color = segments_colors[2]) +
      geom_segment(data = df1, aes(x = x0, y = y0, xend = x1, yend = y1), color = segments_colors[1]) +
      geom_text(data = location_df, aes(label = location, x = location_PC1, y = location_PC2), size = 4) +
      xlab(paste0("PCA1 (", round(subset(PC_summ, PC == "PC1", propvar_SS, drop = TRUE), 2) * 100, "%)")) + 
      ylab(paste0("PCA2 (", round(subset(PC_summ, PC == "PC2", propvar_SS, drop = TRUE), 2) * 100, "%)")) + 
      labs(subtitle = paste0("GGE biplot, ", str_add_space(row$trait), ", ", toupper(row$nursery))) +
      g_plot_mod +
      scale_color_discrete(guide = FALSE) 
    
    g_plot
    
  }) %>% ungroup()


# Combine plots
gge_plot_combine <- plot_grid(plotlist = ggl_biplots$plot, ncol = 2, align = "hv")

ggsave(filename = "ggl_gge_biplot_examples.jpg", plot = gge_plot_combine, path = fig_dir,
       height = 10, width = 12, dpi = 1000)




## Remove genotype points
ggl_biplots1 <- bilinear_toplot %>%
  filter(model == "gge") %>%
  group_by(trait, model, nursery) %>%
  do(plot = {
    row <- .
    df <- row$plotting_df[[1]]
    PC_summ <- row$PC_summ[[1]]
    
    
    ## Calculate average PCs for location
    location_df <- aggregate(cbind(PC1, PC2) ~ location + group, data = df, FUN = mean, subset = group == "environment") %>%
      rename_at(vars(PC1, PC2), ~paste0("location_", .)) %>%
      mutate(x0 = 0, x1 = location_PC1, y0 = 0, y1 = location_PC2)
    # Add location information to df
    df1 <- subset(df, group == "environment") %>% 
      left_join(., location_df) %>%
      mutate(x0 = location_PC1, x1 = PC1, y0 = location_PC2, y1 = PC2)
    
    ## GGL-GGE biplot
    g_plot <- df %>%
      ggplot(aes(x = PC1, y = PC2, color = group)) +
      # geom_point(data = subset(df, group == "line_name"), size = 0.5, shape = 3) +
      geom_point(data = df1, color = segments_colors[1], size = 0.5) +
      geom_segment(data = location_df, aes(x = x0, y = y0, xend = x1, yend = y1), color = segments_colors[2]) +
      geom_segment(data = df1, aes(x = x0, y = y0, xend = x1, yend = y1), color = segments_colors[1]) +
      geom_text(data = location_df, aes(label = location, x = location_PC1, y = location_PC2), size = 4) +
      xlab(paste0("PCA1 (", round(subset(PC_summ, PC == "PC1", propvar_SS, drop = TRUE), 2) * 100, "%)")) + 
      ylab(paste0("PCA2 (", round(subset(PC_summ, PC == "PC2", propvar_SS, drop = TRUE), 2) * 100, "%)")) + 
      labs(subtitle = paste0("GGE biplot, ", str_add_space(row$trait), ", ", toupper(row$nursery))) +
      g_plot_mod +
      scale_color_discrete(guide = FALSE) 
    
    g_plot
    
  }) %>% ungroup()


# Combine plots
gge_plot_combine <- plot_grid(plotlist = ggl_biplots1$plot, ncol = 2, align = "hv")

ggsave(filename = "ggl_gge_biplot_examples_nogenotypes.jpg", plot = gge_plot_combine, path = fig_dir,
       height = 10, width = 12, dpi = 1000)



