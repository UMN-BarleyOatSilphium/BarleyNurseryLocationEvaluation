## Barley Nursery Analysis
## 
## Miscellaneous cleaning script
## 
## Author: Jeff Neyhart
## Last modified: 26 Feb. 2020
## 


# Base script
proj_dir <- getwd()
source(file.path(proj_dir, "startup.R"))

# Other packages
library(jsonlite)
library(neyhart)

# String of us states and canadian provinces
# State names and Canadian provinces to remove
province_df <- rvest::html_table(rvest::html_session("https://www12.statcan.gc.ca/census-recensement/2011/ref/dict/table-tableau/table-tableau-8-eng.cfm"))[[1]][-1,]
province.names <- province_df$`Province/Territory`
province.abb <- province_df$`Internationally approved alpha code (Source: Canada Post)`

state.prov.abb <- c(state.abb, province.abb)
state.prov.names <- c(state.name, province.names)


# Use open streen map to find lat/long coordinates for each site

location_meta <- trial_metadata %>%
  distinct(location, state) %>%
  arrange(location) %>%
  mutate(out = list(NULL))

# Base url for sending queries
base_url <- "https://nominatim.openstreetmap.org/search.php?"

# Loop over locations
i = 1
while (is.null(location_meta$out[[i]])) {
  
  # Convert city and state into long form
  city <- str_replace_all(str_add_space(location_meta$location[i]), " ", "%20")
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
 
  i = i + 1
}


location_meta1 <- location_meta %>% 
  unnest(out) %>% 
  mutate_if(is.character, parse_guess)

## Recombine with trial metadata, save
trial_metadata1 <- left_join(trial_metadata, location_meta1)

write_csv(x = trial_metadata1, path = file.path(data_dir, "nursery_trial_metadata_use1.csv"), na = "")










