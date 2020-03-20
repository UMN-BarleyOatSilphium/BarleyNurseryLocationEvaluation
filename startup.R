## S2MET Prediction Models source script
## 
## A script that automatically loads the data relevant for the S2MET project

library(sommer)
library(tidyverse)
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
trial_metadata <- read_csv(file = file.path(data_dir, "nursery_trial_metadata_use1.csv")) %>%
  # Dryland == rainfed
  mutate(management = ifelse(management == "irrigated", "irrigated", "rainfed"))

# Line metadata
line_metadata <- read_csv(file = file.path(data_dir, "nursery_entry_metadata_use.csv"))

# Pedigree relationship matrix
load(file = file.path(data_dir, "nursery_pedigree_relmat.RData"))


