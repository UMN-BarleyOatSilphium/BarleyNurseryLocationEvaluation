## Code to harmonize all tidy from the historical barley nurseries
## 
## This script will do the following:
## 1. Use corrected line names from the pedigree parsing and match line names
## in the trial data
## 2. Verify pedigree and trial data line names match
## 
## Author: Jeff Neyhart
## 
## 


# Packages to use
packages <- c("tidyverse", "readxl", "lubridate", "measurements")
invisible(lapply(packages, library, character.only = TRUE))

## Source a script for relevant functions
source(file.path(getwd(), "functions.R"))

## Use my package for adding spaces to strings
library(neyhart)


# Trait metadata
trait_metadata <- read_csv("C:/GoogleDrive/BarleyLab/Breeding/PhenotypicData/Metadata/trait_metadata.csv") %>%
  rename_all(tolower) %>%
  select(trait, nickname, class, min, max, unit, unit.info)


# Directories containing data
shared_drive_dir <- c(paste0(LETTERS, ":/AGRO/BARLEY_LAB"), paste0(LETTERS, ":/BARLEY_LAB")) %>%
  subset(., dir.exists(.)) %>%
  first()

# Project directory
proj_dir <- getwd()


nursery_dir <- file.path(shared_drive_dir, "/Breeding/BARLEY/HistoricalNurseryData/")
raw_dir <- file.path(nursery_dir, "RawData")
agro_dir <- file.path(raw_dir, "Agronomic")
maltq_dir <- file.path(raw_dir, "Quality")
tidy_dir <- file.path(nursery_dir, "TidyData")
extr_dir <- file.path(nursery_dir, "ExtractedData")

# Directory of data for use in this project
data_dir <- file.path(proj_dir, "Data")


#######################
# Link original and corrected line names
#######################


## Load the entry metadata
load(file.path(extr_dir, "nursery_raw_entry_metadata.RData"))

## Load the corrected pedigree information
nursery_entry_data_clean <- read_tsv(file = file.path(tidy_dir, "nursery_line_info_clean_standardized_step9.txt")) %>%
  rename_all(tolower) %>%
  arrange(name) %>%
  mutate(aliases = str_split(aliases, ","),
         line_name_use = name)


## Load the nursery phenotype data
nursery_pheno_data <- read_csv(file = file.path(tidy_dir, "nursery_trait_data_tidy.csv"), na = c("", "NA")) %>%
  spread(trait, value) %>%
  mutate(line_name = str_trim(line_name))


## Distinct pheno line names
nursery_pheno_lines <- nursery_pheno_data %>%
  distinct(line_name) %>%
  arrange(line_name)

## Match these lines to line names in the nursery_entry_data_clean df
nursery_pheno_lines_lineNameMatch <- nursery_pheno_lines %>%
  left_join(., nursery_entry_data_clean, by = c("line_name" = "name"))

## Match no-matches from this first round to aliases in the nursery_entry_data_clean df
nursery_pheno_lines_aliasMatch <- nursery_pheno_lines_lineNameMatch %>%
  filter(is.na(line_name_use)) %>%
  select(line_name) %>%
  mutate(out = map(line_name, function(nm) filter(nursery_entry_data_clean, map_lgl(aliases, ~nm %in% .))))


## Find lines that didn't match either
nursery_pheno_lines_aliasMatch %>% 
  filter(map_lgl(out, ~nrow(.) == 0))
  
  
## Only 3; drop them? Sure


## Combine matching line names; rename if necessary; save this metadata
nursery_entry_metadata_use <- bind_rows(
  filter(nursery_pheno_lines_lineNameMatch, !is.na(line_name_use)),
  unnest(nursery_pheno_lines_aliasMatch)
) %>%
  select(line_name, pedigree, line_name_use, aliases) %>%
  mutate(aliases = map_chr(aliases, ~paste0(., collapse = ", ")))


## Extract other metadata information (i.e. row type, use)
nursery_entry_metadata_use1 <- nursery_entry_metadata_use %>% 
  left_join(., distinct(nursery_entry_data_expanded, line_name, row_type, end_use)) %>%
  group_by(line_name)

nursery_entry_metadata_use2 <- nursery_entry_metadata_use1 %>%
  filter(n() == 1) %>%
  bind_rows(., filter(nursery_entry_metadata_use1, n() > 1) %>% mutate(end_use = "MALT") %>% distinct()) %>%
  ungroup()


## Rename entries in the pheno df
nursery_pheno_data1 <- nursery_pheno_data %>% 
  inner_join(., select(nursery_entry_metadata_use2, line_name, line_name_use)) %>% 
  select(-line_name) %>% 
  select(trial:environment, line_name = line_name_use, names(.))



## Write both df to csv's
nursery_entry_metadata_use2 %>% 
  rename(alias = line_name) %>%
  select(line_name = line_name_use, names(.)) %>%
  write_csv(x = ., path = file.path(data_dir, "nursery_entry_metadata_use.csv"), na = "")

nursery_pheno_data1 %>%
  write_csv(x = ., path = file.path(data_dir, "nursery_phenotype_data_use.csv"), na = "")

# Copy the trial metadata
file.copy(from = file.path(tidy_dir, "nursery_trial_metadata_use.csv"), 
          to = file.path(data_dir, "nursery_trial_metadata_use.csv"), overwrite = TRUE)



##





