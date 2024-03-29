## Code to tidy the agronomic and quality data from historical barley nurseries
## 
## Author: Jeff Neyhart
## Date: 8 October 2019
## 
## 

## This script needs to do the following:
## 1. Read in and parse the malting quality data
## 2. Combine MVN agronomic data with the environment codes
## 3. Reformat and tidy the WRN agronomic data
## 4. Combine trait data
## 5. Clean up the pedigree data and search for additional alias and pedigrees from T3



# Packages to use
packages <- c("tidyverse", "readxl", "lubridate", "measurements", "jsonlite", "neyhart")
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




##############################
### Tidy trait data
##############################

# Load extracted maltq data
load(file.path(extr_dir, "maltq_data_extraction.RData"))

## Load extracted data from WRN and MVN agronomic reports
load(file.path(extr_dir, "wrn_extracted_metadata.RData"))
load(file.path(extr_dir, "mvn_extracted_metadata.RData"))



### MVN ###

## Extract trait data from the list
mvn_agro_data <- mvn_extracted_data_list %>%
  map(`[`, c("year", "trait_data")) %>%
  transpose() %>%
  # Add year
  pmap_df(~mutate(.y, year = .x)) %>%
  # Remove number and table name
  select(-number, -table_name) %>%
  # Convert some strings to title case
  mutate_at(vars(trait, location), str_to_title) %>%
  # Rename traits
  mutate(trait = case_when(
    trait == "Yield" ~ "GrainYield",
    trait == "Height" ~ "PlantHeight",
    trait == "Heading" ~ "HeadingDate",
    trait == "Ripe" ~ "MaturityDate",
    TRUE ~ str_remove_all(trait, " ")
  ))


## Custom location renaming df
## This should be consistent for trait data and station data
mvn_location_rename <- read_excel(path = file.path(extr_dir, "mvn_extracted_location_rename.xlsx"), sheet = "data_rename")


# Rename locations
mvn_agro_data1 <- mvn_agro_data %>%
  mutate(location = map_chr(location, ~subset(mvn_location_rename, str_detect(location, .), rename, drop = TRUE))) %>%
  left_join(., select(mvn_location_rename, location = rename, state))


## Row bind and tidy
mvn_agro_data_df <- mvn_agro_data1 %>%
  ## Parse all columns
  mutate_if(is.character, parse_guess) %>%
  ## Add nursery and management
  mutate(nursery = "mvn", management = "rainfed",
         location = str_to_title(location) %>% str_remove_all(., " ")) %>%
  filter(!is.na(value)) %>%
  select(line_name, location, year, nursery, management, trait, value) %>%
  arrange(year, location, trait, line_name)



## Load station data renaming df
mvn_location_rename <- read_excel(path = file.path(extr_dir, "mvn_extracted_location_rename.xlsx"), sheet = "station_rename")

# Read in the station data
mvn_station_data <- read_csv(file = file.path(extr_dir, "mvn_extracted_station_info.csv"))

# Rename locations in the station data df
# Parse some other columns
mvn_station_data_df <- mvn_station_data %>%
  mutate(location = str_to_title(location) %>% str_remove_all(., " "),
         year = as.numeric(str_extract(report, "[0-9]{4}")),
         management = "rainfed",
         nursery = "mvn") %>%
  left_join(., mvn_location_rename) %>%
  select(-location) %>%
  rename(location = rename) %>%
  mutate_at(vars(contains("date")), ymd) %>%
  select(nursery, location, state, year, management, names(.), -contains("number")) %>%
  filter_at(vars(planting_date:lsd_bua), any_vars(!is.na(.))) %>%
  # If replications is NA, assume 3
  mutate(replications = ifelse(is.na(replications), 3, replications),
         # Convert lsd_bua to kgha
         lsd_kgha1 = lsd_bua * 53.78595,
         lsd_kgha = coalesce(lsd_kgha, lsd_kgha1)) %>%
  rename(lsd_yield = lsd_kgha) %>% select(-lsd_kgha1, -lsd_bua)
         

# Create a trait metadata table; trait name and replications, lsd (if applicable)
mvn_agro_trait_metadata <- mvn_agro_data_df %>% 
  distinct_at(vars(-line_name, -value)) %>%
  left_join(., select(mvn_station_data_df, nursery, location, year, management, replications, lsd_yield)) %>%
  # Full join with complete list of trial/trait combinations for the nursery
  full_join(., distinct(mvn_agro_data_df, location, year, nursery, management, trait)) %>%
  mutate(lsd = ifelse(trait == "GrainYield", lsd_yield, NA),
         replications = ifelse(is.na(replications), 3, replications)) %>%
  select(-lsd_yield)
  
mvn_station_data_df <- select(mvn_station_data_df, -lsd_yield)







### WRN ###

## Extract trait data from the list
wrn_agro_data <- wrn_extracted_data_list %>%
  map(`[`, c("year", "trait_data")) %>%
  transpose() %>%
  # Add year
  pmap_df(~mutate(.y, year = .x)) %>%
  # Remove number or ent.
  select(-number, -ent.) %>%
  # Check traits
  mutate(trait = ifelse(is.na(trait), ifelse(str_detect(table_name, "YIELD"), "GRAIN YIELD", 
                                             ifelse(str_detect(table_name, "HEIGHT"), "PLANT HEIGHT", "NA" )), trait))


## Load trait data renaming df
wrn_location_rename <- read_excel(path = file.path(extr_dir, "wrn_extracted_location_rename.xlsx"), sheet = "data_rename")

# Rename locations in the trait data
wrn_agro_data1 <- wrn_agro_data %>%
  mutate(location = ifelse(str_detect(location, " "), str_remove_all(str_to_title(location), " "), location)) %>%
  left_join(., wrn_location_rename) %>%
  select(-location) %>%
  rename(location = rename)
  

# Modify elements
wrn_agro_data_df <- wrn_agro_data1 %>%
  # Coalesce line name, cultivar, or selection
  mutate(line_name = coalesce(selection, cultivar., line_name),
         # Extract dryland or irrigated
         management = case_when(
           !is.na(management) ~ management,
           str_detect(table_name, "DRYLAND") ~ "rainfed",
           TRUE ~ "irrigated"
         ),
         # Remove spaces and capital letters from trait
         trait = str_add_space(trait) %>% str_to_title() %>% str_remove_all(" "),
         # Add nursery
         nursery = "wrn") %>%
  # Reparse
  mutate_if(is.character, parse_guess) %>%
  # There is a NA in line name for test weight in 2001 - this is CDC Copeland
  mutate(line_name = ifelse(year == 2001 & trait == "TestWeight" & is.na(line_name), "CDC Copeland", line_name)) %>%
  filter(!is.na(value)) %>%
  select(line_name, location, state, year, nursery, management, trait, value) %>%
  arrange(year, location, trait, line_name)


# Read in the station data
wrn_station_data <- read_csv(file = file.path(extr_dir, "wrn_extracted_station_info.csv"))
# Read in the station trait metadata
wrn_trait_station_data <- read_csv(file = file.path(extr_dir, "wrn_extracted_trial_location_metadata.csv")) %>%
  mutate(location = ifelse(str_detect(location, " "), str_remove_all(str_to_title(location), " "), location),
         year = str_extract(report, "[0-9]{2,4}"),
         year = ifelse(nchar(year) == 2, paste0("20", year), year) %>% as.numeric())


# Rename locations in the station data df
# Parse some other columns
wrn_station_data_df <- wrn_station_data %>%
  mutate(location = ifelse(str_detect(location, " "), str_remove_all(str_to_title(location), " "), location),
         year = str_extract(report, "[0-9]{2,4}"),
         year = ifelse(nchar(year) == 2, paste0("20", year), year) %>% as.numeric()) %>%
  left_join(., wrn_location_rename) %>%
  select(-location) %>%
  rename(location = rename) %>%
  left_join(., distinct(wrn_agro_data_df, location, year, nursery, management)) %>%
  mutate_at(vars(contains("date")), ymd) %>%
  filter_at(vars(nursery, management), any_vars(!is.na(.)))

# Create a trait metadata table; trait name and replications, lsd (if applicable)
wrn_agro_trait_metadata <- wrn_trait_station_data %>% 
  left_join(., wrn_location_rename) %>%
  select(-location) %>%
  rename(location = rename) %>%
  filter(!is.na(location)) %>%
  select(location, year, nursery, management, trait, replications, lsd) %>%
  # Full join with complete list of trial/trait combinations for the nursery
  full_join(., distinct(wrn_agro_data_df, location, year, nursery, management, trait)) %>%
  mutate(replications = ifelse(is.na(replications), 3, replications))



## Combine agronomic data
agro_data_tidy <- bind_rows(mvn_agro_data_df, wrn_agro_data_df) %>%
  mutate_if(is.character, parse_guess) %>%
  # Create trial names
  mutate(trial = paste0(toupper(nursery), "_", year, "_", location, "_", str_to_title(management)),
         # Remove numbers from location
         location = str_remove(location, "[0-9]$")) %>%
  # Remove state
  select(-state)

## Combine station information
agro_station_data_tidy <- bind_rows(mvn_station_data_df, wrn_station_data_df) %>%
  arrange(nursery, year, location) %>%
  # Create trial names
  mutate(trial = paste0(toupper(nursery), "_", year, "_", location, "_", str_to_title(management)),
         # Remove numbers from location
         location = str_remove(location, "[0-9]$"),
         # Brandon in is Manitoba
         state = ifelse(location == "Brandon", "MB", state))
         



### Malting Quality ###

maltq_trait_data <- maltq_data_tidy1 %>%
  mutate(nursery = tolower(nursery)) %>%
  mutate_if(is.character, parse_guess) %>%
  ## All wrn malt quality trials were irrigated
  mutate(management = ifelse(nursery == "mvn", "rainfed", "irrigated")) %>%
  select(line_name, location, year, nursery, management, trait, value)


## Extract location data
maltq_station_data <- maltq_trait_data %>%
  distinct(nursery, location, year, management)


## Load location renaming DF
maltq_location_rename <- read_excel(path = file.path(extr_dir, "maltq_extracted_location_rename.xlsx"), sheet = "station_rename")


## Rename location in the trait data df
maltq_trait_df <- maltq_trait_data %>%
  left_join(., maltq_location_rename) %>%
  select(-location, -state) %>%
  rename(location = rename) %>%
  # Create trial names
  mutate(trial = paste0(toupper(nursery), "_", year, "_", location, "_", str_to_title(management)),
         # Remove numbers from location
         location = str_remove(location, "[0-9]$"))

## Rename location in the station data df
maltq_station_df <- maltq_station_data %>%
  left_join(., maltq_location_rename) %>%
  select(-location) %>%
  rename(location = rename) %>%
  # Create trial names
  mutate(trial = paste0(toupper(nursery), "_", year, "_", location, "_", str_to_title(management)),
         # Remove numbers from location
         location = str_remove(location, "[0-9]$"),
         # Brandon in is Manitoba
         state = ifelse(location == "Brandon", "MB", state))


## No need to edit the location strings for capitalization or sentence case





## Combine AGRO and MALTQ trait information
nursery_data_tidy <- bind_rows(agro_data_tidy, maltq_trait_df) %>%
  select(trial, location, year, nursery, management, trait, line_name, value) %>%
  distinct()

## Pull out trial metadata to merge
nursery_trial_metadata1 <- nursery_data_tidy  %>% 
  distinct(trial, location, year, nursery, management)

## Add metadata information from the station dfs
nursery_trial_metadata2 <- nursery_trial_metadata1 %>%
  left_join(., full_join(agro_station_data_tidy, maltq_station_df)) %>%
  arrange(year, location, nursery)

nursery_trial_metadata_tidy <- nursery_trial_metadata2 %>%
  select(-state) %>%
  # Add state information for those missing
  left_join(., filter(distinct(nursery_trial_metadata2, location, state), !is.na(state)), by = "location") %>%
  mutate(state = ifelse(location == "Kernen", "SK", state)) %>%
  arrange(trial, location, year, nursery, management, planting_date, harvest_date) %>%
  # Remove duplicates if planting date is available
  group_by(trial, location, year, nursery, management) %>%
  slice(1) %>% ungroup() %>%
  select(-replications)

# Create a trait/trial metadata df
nursery_trait_metadata_tidy <- bind_rows(mvn_agro_trait_metadata, wrn_agro_trait_metadata) %>%
  # Add maltq station data
  bind_rows(., mutate(distinct_at(maltq_trait_df, vars(location, year, nursery, management, trait)), replications = 1)) %>%
  # Create trial names
  mutate(trial = paste0(toupper(nursery), "_", year, "_", location, "_", str_to_title(management)),
         # Remove numbers from location
         location = str_remove(location, "[0-9]$")) %>%
  left_join(distinct(nursery_trial_metadata_tidy, trial, location, year, nursery, management), .)
  




##############################
### Environment metadata
##############################

# The output of this section should be a data.frame with: 
# 1. location
# 2. year
# 3. environment code
# 4. nursery 
# 5. management
# 6. planting date (if present)
# 7. harvest date (if present)
# 8. latitude (to be guessed by GEMS)
# 10. longitude (to be guessed by GEMS)

## Make a common trial key; change some location names
nursery_trial_metadata_key <- nursery_trial_metadata_tidy %>%
  distinct(trial, location, state, year, nursery, management) %>%
  # Make location corrections
  mutate(location = case_when(
    location == "Gennesse" ~ "Genesee", 
    location == "Moorehead" ~ "Moorhead",
    location == "Stirling" ~ "Sterling",
    location == "Tammany" ~ "Lewiston",
    # Kernen and Sasktatoon are the same locations
    location == "Kernen" ~ "Saskatoon",
    TRUE ~ location)
  ) %>%
  # Remove spaces from location, create environment code and trial name and RorI for rain/irr
  mutate(mgmt = ifelse(str_detect(trial, "Rainfed"), "R", "I"),
         trial_new = paste0(toupper(nursery), "_", year, "_", location, "_", str_to_title(management)),
         environment = paste0(toupper(abbreviate(location, 3)), str_sub(year, 3, 4), mgmt)) %>%
  select(trial, trial_new, location, state, year, environment, nursery, management) %>%
  distinct() 


# Create a new metadata df
nursery_trial_metadata_tidy1 <- nursery_trial_metadata_tidy %>%
  select(-location, -year, -nursery, -management, -state) %>%
  left_join(., nursery_trial_metadata_key) %>%
  select(trial = trial_new, location, state, year, environment, nursery, management, names(.), -trial)


## Add environment code to nursery_data_tidy
nursery_data_tidy1 <- nursery_data_tidy %>%
  select(-location, -year, -nursery, -management) %>%
  left_join(., nursery_trial_metadata_key) %>%
  select(trial = trial_new, trait, line_name, value)

# And to the trait metadata
nursery_trait_metadata_tidy1 <- nursery_trait_metadata_tidy %>%
  select(-location, -year, -nursery, -management) %>%
  left_join(., nursery_trial_metadata_key) %>%
  select(trial = trial_new, trait, replications, lsd)



## Determine approximate lat/long coordinates of locations ##

# String of us states and canadian provinces
# State names and Canadian provinces to remove
province_df <- rvest::html_table(rvest::html_session("https://www12.statcan.gc.ca/census-recensement/2011/ref/dict/table-tableau/table-tableau-8-eng.cfm"))[[1]][-1,]
province.names <- province_df$`Province/Territory`
province.abb <- province_df$`Internationally approved alpha code (Source: Canada Post)`

state.prov.abb <- c(state.abb, province.abb)
state.prov.names <- c(state.name, province.names)




# Use open streen map to find lat/long coordinates for each site

location_meta <- nursery_trial_metadata_tidy1 %>%
  distinct(location, state) %>%
  arrange(location) %>%
  mutate(out = list(NULL))

# Base url for sending queries
base_url <- "https://nominatim.openstreetmap.org/search.php?"

# Loop over locations
for (i in seq_len(nrow(location_meta))) {
  
  # Convert city and state into long form
  city <- location_meta$location[i]
  city <- str_replace_all(ifelse(city == "McVille", city, str_add_space(city)), " ", "%20")
  state <- str_replace_all(state.prov.names[state.prov.abb == location_meta$state[i]], " ", "%20")
  country <- ifelse(location_meta$state[i] %in% state.abb, "US", "CA")
  
  # Create the query string
  query_url <- paste0(
    base_url,
    "city=", city,
    "&state=", state,
    "&countrycodes=", country,
    "&limit=9&format=json"
  )
  
  if (city == "Kernen") {
    # Convert lat/long to tibble
    location_meta$out[[i]] <- tibble(lat = "52.150705", long = "-106.543686")
    
    
  } else {
    
    # Send query and read json output
    query_out <- read_json(path = query_url)
    
    # Convert lat/long to tibble
    location_meta$out[[i]] <- tibble(lat = query_out[[1]]$lat, long = query_out[[1]]$lon)
    
  }

}


location_meta1 <- location_meta %>% 
  unnest(out) %>% 
  mutate_if(is.character, parse_guess)

## Recombine with trial metadata, save
nursery_trial_metadata_tidy2 <- left_join(nursery_trial_metadata_tidy1, location_meta1)


## Write a table of locations as a CSV
write_csv(x = nursery_trial_metadata_tidy2, path = file.path(tidy_dir, "nursery_trial_metadata_use.csv"), na = "")




##############################
### Combine and validate trait data
##############################



## Tasks:
## 1. Clean up line names
## 2. Rename traits; make sure values make sense for those traits
## 3. Rename trials using the new location information from above

nursery_trial_metadata_tidy3 <- nursery_trial_metadata_tidy2 %>%
  mutate(trial_new = paste0(toupper(nursery), "_", year, "_", location, "_", str_to_title(management)))

# Combine maltq and agro trait data
nursery_trait_data <- nursery_data_tidy1 %>%
  # Rename trials, location, environment, management, etc
  select(trial, trait, line_name, value) %>%
  left_join(., nursery_trial_metadata_tidy3) %>%
  select(trial = trial_new, nursery, location, year, management, environment, trait, line_name, value)

nursery_trait_data %>% group_by_at(vars(-value)) %>% filter(n() > 1)  


## Look for missing trials
nursery_trait_data %>% filter(is.na(trial)) %>% distinct(trial, location, year) %>% as.data.frame()


# Clean up line names manually; save originals
nursery_trait_data_originals <- nursery_trait_data %>%
  select(line_name_original = line_name) %>% distinct() %>%
  # Parse; remove white space
  mutate_all(~parse_guess(x = ., trim_ws = TRUE)) %>%
  # Capitalize; extract alias
  mutate(line_name_edited = toupper(line_name_original),
         # Remove preceding numbers
         line_name_edited = str_remove(line_name_edited, "^[0-9]{1,2} "),
         # Remove _BWD
         line_name_edited = str_remove(line_name_edited, "_BWD$"),
         # Remove space and Check
         line_name_edited = str_remove_all(line_name_edited, " ") %>% str_remove(., "\\,CHECK"),
         # Remove preceding hyphens
         line_name_edited = str_remove(line_name_edited, "^-"),
         alias = str_extract_all(line_name_edited, "\\(.*\\)") %>% map(., ~str_remove_all(., "\\(|\\)")),
         line_name_edited = str_remove(line_name_edited, "\\(.*\\)") %>% str_remove("\\(.*"),
         line_name_edited = str_trim(line_name_edited)) %>%
  # Filter out line names with /
  filter(str_detect(line_name_edited, "/", negate = TRUE))

## Plug these through pedTools to identify duplicates/aliases
nursery_trait_data_originals %>%
  mutate(alias = map_chr(alias, ~paste0(., collapse = ", ")) %>% parse_character(),
         pedigree = NA) %>%
  select(line_name_edited, pedigree, alias) %>%
  write_tsv(x = ., path = file.path(extr_dir, "nursery_line_names_for_pedTools.txt"), na = "")


## Reload the output of pedtools
nursery_line_names_renamed <- read_tsv(file = file.path(proj_dir, "Scripts/pedtools/barleyNurseryLines/nursery_line_names_renamed_step8.txt")) %>%
  rename_all(~str_remove_all(., "\\(.*\\)")) %>%
  rename_all(tolower) %>%
  ## Split alias
  mutate(aliases = str_split(string = alias, pattern = ",")) %>%
  unnest(aliases) %>%
  distinct(aliases, corrected_name) %>%
  filter(aliases != corrected_name)


## Change line names
nursery_trait_data1 <- nursery_trait_data %>%
  left_join(., distinct(nursery_trait_data_originals, line_name_original, line_name_edited), by = c("line_name" = "line_name_original")) %>%
  # mutate(line_name = ifelse(is.na(corrected_name), line_name, corrected_name))
  select(trial, location, year, environment, trait, line_name = line_name_edited, value) %>%
  distinct() 
# 
# %>%
#   left_join(., nursery_line_names_renamed, by = c("line_name" = "aliases")) %>%
#   mutate(line_name = ifelse(is.na(corrected_name), line_name, corrected_name)) %>%
#   select(-corrected_name) %>%
#   distinct()




### Validate trait names and range
## Modify trait names

# Pull unique traits
unique_traits <- nursery_trait_data1 %>% 
  distinct(trait) %>% 
  arrange(trait) %>%
  rename(trait_old = trait) %>%
  # Find matching ontologies by nickname
  left_join(., trait_metadata, by = c("trait_old" = "nickname"))

# Rename traits with missing matches
trait_rename <- unique_traits %>%
  filter(is.na(trait)) 
# Create vector to rename traits
trait_rename_vct <- c("PlumpBarley" = "PlumpGrain", "ST" = "SolubleProteinTotalProtein", 
                      "ThinGrain" = "ThinGrains", "ThinBarley" = "ThinGrains", "Viscosity" = "WortViscosity")
# Rename
trait_rename1 <- trait_rename %>% 
  mutate(trait_rename = str_replace_all(trait_old, trait_rename_vct)) %>%
  select(trait_old, trait_rename) %>%
  filter(trait_old != trait_rename)

# Modify trait names in unqiue traits
unique_traits1 <- unique_traits %>% 
  left_join(., trait_rename1) %>% 
  mutate(trait_new = ifelse(is.na(trait), trait_rename, trait_old)) %>%
  select(trait_new, trait_old)

unique_traits2 <- unique_traits1 %>%
  # Find matching ontologies by nickname
  left_join(., trait_metadata, by = c("trait_new" = "nickname")) %>%
  filter(!is.na(trait)) %>%
  rename(nickname = trait_old) %>%
  # Parse min/max
  mutate_at(vars(min, max), parse_guess)

## Rename traits in the dataset
nursery_trait_data2 <- nursery_trait_data1 %>%
  left_join(., select(unique_traits1, trait = trait_old, trait_new)) %>%
  select(-trait) %>% rename(trait = trait_new) %>%
  filter(!is.na(trait))

# Rename in the trait metadata
nursery_trait_metadata_tidy2 <- nursery_trait_metadata_tidy1 %>%
  left_join(., select(unique_traits1, trait = trait_old, trait_new)) %>%
  select(-trait) %>% rename(trait = trait_new) %>%
  filter(!is.na(trait))





## Look for traits with obsene distributions
par(mfrow = c(4, 5))
for (tr in unique(nursery_trait_data2$trait)) {
  hist(subset(nursery_trait_data2, trait == tr, value, drop = T), main = tr,
       xlab = NULL)
}
par(mfrow = c(1,1))


# Vector of traits to investigate
traits_bad_units <- c("GrainYield", "TestWeight", "HeadingDate", "MaturityDate", "PlantHeight", "BetaGlucan", "PlumpGrain")



## One trait at a time
# Grain Yield

# Visualize
hist(subset(nursery_trait_data2, trait == "GrainYield" & value < 20, value, drop = T))
hist(subset(nursery_trait_data2, trait == "GrainYield" & value > 20 & value < 600, value, drop = T))
range(subset(nursery_trait_data2, trait == "GrainYield" & value > 20 & value < 600, value, drop = T))

  
# Look at distinct nuseries to understand
nursery_trait_data2 %>%
  filter(trait == "GrainYield") %>%
  # filter(value < 15) %>% # < 15 is likely Mg ha
  # filter(between(value, 15, 200)) %>% # > 15 and < 200 is likely bu ac
  filter(value < 200) %>%
  mutate(nursery = str_extract(trial, "^[A-Z]{3}")) %>%
  # distinct(nursery, year)
  group_by(trial) %>%
  summarize(n = n()) %>%
  arrange(desc(n))








## Edit grain yield
grain_yield_edit <- nursery_trait_data2 %>%
  filter(trait == "GrainYield") %>%
  mutate(value = case_when(
    value < 20 ~ value * 1000, # Mg ha to kg to ha,
    between(value, 25, 300) ~ value / (0.018587),
    TRUE ~ value
  ))

hist(grain_yield_edit$value)

# Heading Date
# If the trial is from MVN, heading date is in days after May 31
heading_date_edit <- nursery_trait_data2 %>% 
  filter(trait == "HeadingDate") %>%
  mutate(may31julian = yday(ymd(paste0(year, "0531"))),
         value = case_when(
           value < 100 ~ value + may31julian,
           TRUE ~ value
         )) %>%
  select(-may31julian)

hist(heading_date_edit$value)


# Maturity Date
maturity_date_edit <- nursery_trait_data2 %>% 
  filter(trait == "MaturityDate") %>%
  # Remove 0 value
  filter(value != 0) %>%
  mutate(may31julian = yday(ymd(paste0(year, "0531"))),
         value = value + may31julian) %>%
  select(-may31julian)

hist(maturity_date_edit$value)

# TestWeight
# if trial is MVN, units are kg/hl
# convert to g / l by multiplying by 10
# 
test_weight_edit <- nursery_trait_data2 %>% 
  filter(trait == "TestWeight") %>%
  mutate(value = case_when(
    str_detect(trial, "MVN") & value < 250 ~ value * 10, # kg/ha to g/l
    str_detect(trial, "WRN") & value < 200 ~ (value / 0.776815) * 10, # lb bu to g/l
    TRUE ~ value
  ))

hist(test_weight_edit$value)


# Plant Height
## List trials where the values are inches, not cm
trials_inches <- c(
  paste0("WRN_2017_", c("Aberdeen", "MosesLake", "Pullman", "Powell"), "_Irrigated"),
  paste0("WRN_2018_", c("Logan"), "_Irrigated")
)

plant_height_edit <- nursery_trait_data2 %>% 
  filter(trait == "PlantHeight") %>%
  mutate(value = ifelse(trial %in% trials_inches, value / 0.393701, value))

hist(plant_height_edit$value)


# BetaGlucan
# Visualize
hist(subset(nursery_trait_data2, trait == "BetaGlucan" & value > 2000, value, drop = T))


# Look at distinct nuseries to understand
nursery_trait_data2 %>%
  filter(trait == "BetaGlucan") %>%
  filter(value > 2000) %>%
  mutate(nursery = str_extract(trial, "^[A-Z]{3}")) %>%
  # distinct(nursery, year)
  group_by(trial) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

## This is real; remove these outliers
beta_glucan_edit <- nursery_trait_data2 %>% 
  filter(trait == "BetaGlucan", value <= 1500) # 1500 is the maximum defined by T3


# PlumpGrain

# Visualize
hist(subset(nursery_trait_data2, trait == "PlumpGrain", value, drop = T))


# Look at distinct nuseries to understand
nursery_trait_data2 %>%
  filter(trait == "PlumpGrain") %>%
  filter(value < 20) %>%
  group_by(trial) %>%
  summarize(n = n())

# Keep this - no edits

plump_grain_edit <- nursery_trait_data2 %>% 
  filter(trait == "PlumpGrain")






## Reassemble the data.frame
nursery_trait_data3 <- nursery_trait_data2 %>%
  filter(! trait %in% traits_bad_units) %>%
  bind_rows(., grain_yield_edit, heading_date_edit, maturity_date_edit, test_weight_edit, plant_height_edit,
            beta_glucan_edit, plump_grain_edit)

## Find duplicates
nursery_trait_data3_dups <- nursery_trait_data3 %>%
  group_by(trial, location, year, environment, trait, line_name) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  ungroup() %>%
  arrange(line_name)


## Aggregate all line_name-trial-trait combinations
nursery_trait_data4 <- nursery_trait_data3 %>%
  group_by_at(vars(-value)) %>%
  summarize(value = mean(value)) %>%
  ungroup() %>%
  filter(!is.na(line_name))




## Visualize again
par(mfrow = c(4, 5))
for (tr in unique(nursery_trait_data3$trait)) {
  hist(subset(nursery_trait_data3, trait == tr, value, drop = T), main = tr,
       xlab = NULL)
}
par(mfrow = c(1,1))



## Save the tidy trait data
write_csv(x = nursery_trait_data4, path = file.path(tidy_dir, "nursery_trait_data_tidy.csv"), na = "")

## Wide format
nursery_trait_data4 %>% 
  spread(trait, value) %>%
  write_csv(path = file.path(tidy_dir, "nursery_trait_data_wide.csv"), na = "")


# Save the trait trial metadata
write_csv(x = nursery_trait_metadata_tidy2, path = file.path(tidy_dir, "nursery_trial_trait_metadata_use.csv"), na = "")




