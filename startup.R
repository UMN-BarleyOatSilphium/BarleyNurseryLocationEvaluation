## S2MET Prediction Models source script
## 
## A script that automatically loads the data relevant for the S2MET project

library(sommer)
library(tidyverse)
library(lubridate)
library(readxl)
library(neyhart)
library(broom)
library(forcats)

## Directories
## This the project root
proj_dir <- here::here()

## Google drive directory
gdrive_dir <- "C:/GoogleDrive"

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")


######
# MSI Source starts here
######


# Source the project functions
source(file.path(proj_dir, "functions.R"))


## 

## Load data

# Phenotype
pheno_dat <- read.csv(file = file.path(data_dir, "nursery_phenotype_data_use.csv"), stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  # Copy year
  mutate(year1 = year) %>%
  mutate_at(vars(trial, location, year, environment, line_name), as.factor) %>%
  # Tidy
  gather(trait, value, matches("[A-Z]", ignore.case = FALSE)) %>% 
  filter(!is.na(value)) %>%
  # Calculate line means for all trait/trial combos
  group_by_at(vars(-value)) %>%
  summarize(value = mean(value)) %>%
  ungroup()


# Trial metadata
trial_metadata <- read_csv(file = file.path(data_dir, "nursery_trial_metadata_use.csv"))

# Trial trait metadata
trial_trait_metadata <- read_csv(file = file.path(data_dir, "nursery_trial_trait_metadata_use.csv"))

# Line metadata 
line_metadata <- read_csv(file = file.path(data_dir, "nursery_entry_metadata_use.csv"))

# Pedigree relationship matrix
load(file = file.path(data_dir, "nursery_pedigree_relmat.RData"))

# Breakdown traits by agronomic or quality
all_traits <- sort(unique(pheno_dat$trait))
agro_traits <- c("GrainYield", "HeadingDate", "Lodging", "MaturityDate", "PlantHeight", "TestWeight")
quality_traits <- setdiff(all_traits, agro_traits)


## Subset more relevant traits
traits_keep <- c("BetaGlucan", "DiastaticPower", "FreeAminoNitrogen", "GrainProtein", "GrainYield", "HeadingDate",
                 "MaltExtract", "PlantHeight", "PlumpGrain", "SolubleProteinTotalProtein",  "TestWeight")

# Vector of traits where positive values are preferable
pos_val_traits <- c("DiastaticPower", "GrainYield", "MaltExtract", "PlumpGrain", "TestWeight")
# Vector of traits where negative values are preferable
neg_val_traits <- c("BetaGlucan", "FreeAminoNitrogen", "GrainProtein", "HeadingDate", "PlantHeight", "SolubleProteinTotalProtein")


# Function to rename terms
f_term_replace <- function(x) str_replace_all(x, c("line_name" = "Genotype", "environment" = "Environment"))

# Function to expand nursery abbreviations
f_nursery_expand <- function(x) c("mvn" = "Mississippi Valley Nursery", "wrn" = "Western Regional Nursery")[x]

# Rename fitness components
f_fitness_comp_replace <- function(x) 
  str_replace_all(x, c("varY" = "Precision", "repAvg" = "Repeatability", "reprAvg" = "Representativeness"))

# Function to rename traits
f_trait_rename <- function(x) str_replace(str_add_space(x), "Protein Total", "Protein/Total")
# Function to abbreviate traits
f_trait_abbr <- function(x) str_replace(abbreviate(f_trait_rename(x), 2), "SPP", "S/T")

