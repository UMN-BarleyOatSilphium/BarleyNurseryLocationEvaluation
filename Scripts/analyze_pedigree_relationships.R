## Barley Nursery Analysis
## 
## Analyze pedigree relationships
## 
## Author: Jeff Neyhart
## 


# Base script
proj_dir <- getwd()
source(file.path(proj_dir, "startup.R"))

# Additional packages
library(paletteer)
library(ggrepel)
library(cowplot)
library(lme4)
library(kinship2)



# Color pallete for environments and line names
term_colors <- set_names(paletteer_d(package = "wesanderson", palette = "Zissou1")[c(1,5)], "line_name", "environment")

# Function to rename terms
f_term_replace <- function(x) str_replace_all(x, c("line_name" = "Genotype", "environment" = "Environment"))




# Calculate pedigree depth ------------------------------------------------

pedigree_table1 <- pedigree_table %>%
  select(name, par1, par2, generation) %>%
  mutate(depth = kindepth(id = name, dad.id = par1, mom.id = par2))

# Distribution of depth
plot(table(pedigree_table1$depth))


























