## Code to extract quality data from both the MVN and WRN nurseries
## 
## This script will extract data from the excel files containing quality data
## for both the Mississippi Valley Nursery and the Western Regional
## Nursery 
## 
## Author: Jeff Neyhart
## 
## Change log:
## 
## 
## 


## Source a script for relevant functions
source(file.path(getwd(), "startup.R"))

## Use my package for adding spaces to strings
library(neyhart)


# Directories containing data
shared_drive_dir <- c(paste0(LETTERS, ":/AGRO/BARLEY_LAB"), paste0(LETTERS, ":/BARLEY_LAB")) %>%
  subset(., dir.exists(.)) %>%
  first()

nursery_dir <- file.path(shared_drive_dir, "/Breeding/BARLEY/HistoricalNurseryData/")
raw_dir <- file.path(nursery_dir, "RawData")
agro_dir <- file.path(raw_dir, "Agronomic")
maltq_dir <- file.path(raw_dir, "Quality")
tidy_dir <- file.path(nursery_dir, "TidyData")
extr_dir <- file.path(nursery_dir, "ExtractedData")




############################
### Malt Quality Data
############################


## List files
maltq_files <- list.files(path = maltq_dir, full.names = TRUE, pattern = "\\.[A-Za-z]{3,4}$") %>%
  str_subset(string = ., pattern = "~", negate = TRUE) %>%
  str_subset(string = ., pattern = "99MVNSAS", negate = TRUE) # Skip one sheet

# Iterate over the list of files
maltq_data_out <- map(maltq_files, ~{
  
  filename <- .
  # Print
  cat("\n\n", basename(filename), "\n")
  
  # Rename extension
  filename1 <- paste0(str_remove(filename, "\\.[A-Za-z]{3,4}$"), tolower(str_extract(filename, "\\.[A-Za-z]{3,4}$")))
  
  ## Extract year
  year_number <- parse_number(basename(filename1)) %>% 
    str_pad(string = ., width = ifelse(nchar(.) <= 2, 2, 4), side = "left", pad = "0")
  
  year <- ifelse(nchar(year_number) == 4, year_number, ifelse(year_number > str_sub(year(today()), 3, 4), paste0("19", year_number), paste0("20", year_number)))
  
  ## Is this MVN or WRN
  nursery <- str_extract(toupper(filename1), pattern = "MVN|MVBN|WRSBN|WBN")
  nursery <- ifelse(str_detect(nursery, "WRSBN|WBN"), "WRN", ifelse(str_detect(nursery, "MVN|MVBN"), "MVN", "NA"))
  
  # Sheet names
  sheets <- excel_sheets(path = filename)
  
  ## Attempt a first read of all sheets, identify the table (if present) that contains pedigree information
  ## Any sheets with "LabNo" is a data sheet
  suppressMessages(first_read <- map(sheets, ~read_excel(path = filename1, sheet = ., col_names = FALSE, n_max = 20, guess_max = 20) %>%
                                       map_chr(~case_when(
                                         any(str_detect(tolower(.), "entries|description"), na.rm = TRUE) ~ "entry_info",
                                         any(str_detect(tolower(.), "average|addition"), na.rm = TRUE) ~ "skip",
                                         any(str_detect(tolower(.), "labno|lab no|lab \\#"), na.rm = TRUE) ~ "datatable",
                                         TRUE ~ as.character(NA)
                                       ) )) %>%
                     map_chr(~ifelse(any(!is.na(.)), as.character(na.omit(.)), "skip")))
  
  
  # If any are true, read in that table
  any_entry_info <- which(first_read == "entry_info")
  
  if (length(any_entry_info) > 0) {
    entry_info_read <- suppressMessages(read_excel(path = filename1, sheet = any_entry_info, col_names = FALSE, trim_ws = TRUE)) %>%
      # Remove blank rows
      filter_all(any_vars(!is.na(.)))
    
    ## Determine which column contains what information
    entry_no <- which(apply(X = entry_info_read, MARGIN = 2, FUN = function(col) any(str_detect(col, "No\\."), na.rm = TRUE)))
    new_entry <- which(apply(X = entry_info_read, MARGIN = 2, FUN = function(col) any(str_detect(col, "^X$|^\\*$"), na.rm = TRUE)))
    cultivar <- which(apply(X = entry_info_read, MARGIN = 2, FUN = function(col) 
      any(str_detect(col, "Cultivar|Selection|Entry|Identification|Number|Name") & str_detect(col, "New", negate = TRUE), na.rm = TRUE)))
    
    # Remove the entry_no and new_entry columns from those detected for cultivar
    cultivar1 <- setdiff(cultivar, c(entry_no, new_entry))
    # Which of these are entry/cultivar/selection
    entry_or_id <- which(apply(X = entry_info_read[,cultivar1,drop = F], MARGIN = 2, FUN = function(col) 
      any(str_detect(col, "Entry|Identification|Number|Name"), na.rm = TRUE)))
    entry_or_id1 <- cultivar1[entry_or_id]
    # which of these is cultivar or selection
    cultivar2 <- setdiff(cultivar1, entry_or_id1)
    
    # pedigree should be the only column containing slashes
    pedigree <- which(apply(X = entry_info_read, MARGIN = 2, FUN = function(col) any(str_detect(col, "/"), na.rm = TRUE)))
    
    # Others
    row_type <- which(apply(X = entry_info_read, MARGIN = 2, FUN = function(col) any(str_detect(col, "Rowed|Type"), na.rm = TRUE)))
    end_use <- which(apply(X = entry_info_read, MARGIN = 2, FUN = function(col) any(str_detect(col, "Use|Grade"), na.rm = TRUE)))
    cooperator <- which(apply(X = entry_info_read, MARGIN = 2, FUN = function(col) any(str_detect(col, "Cooperator"), na.rm = TRUE)))
    
    pedigree1 <- setdiff(pedigree, c(end_use, cooperator))
    
    ## Gather information from the df
    entry_info <- entry_info_read %>%
      select(number = entry_no, line_name = entry_or_id1, breeding_name = cultivar2, pedigree = pedigree1, row_type = row_type, end_use = end_use) %>%
      filter(nchar(number) <= 2)
    
    if (! "breeding_name" %in% names(entry_info)) {entry_info$breeding_name <- NA}
    if (! "line_name" %in% names(entry_info)) {entry_info$line_name <- NA}
    
    
    entry_info1 <- entry_info %>%
      ## line_name becomes alias if breeding name is present
      mutate(alias = ifelse(!is.na(breeding_name), line_name, NA),
             line_name = ifelse(is.na(breeding_name), line_name, breeding_name),
             number = parse_number(number)) %>%
      select(-breeding_name)
    
  } else {
    entry_info1 <- NULL
    
  }
  
  
  ## Sheets with a space are locations
  location_sheets <- sheets[first_read == "datatable"]
  
  
  # Read and parse the data
  maltq_read <- suppressMessages(read_maltq(file = filename1, sheets = location_sheets))
  
  ## Extract data and store
  maltq_data <- maltq_read$maltq.data %>%
    rename_all(~str_remove_all(., "\\.") %>% make.names(., unique = TRUE)) %>%
    rename(Location = sheet) %>%
    mutate(Nursery = nursery,
           Year = as.character(year))
  
  ## Return a list
  list(maltq_data = maltq_data, entry_info = entry_info1)
  
})

## Add names
names(maltq_data_out) <- basename(maltq_files)



## Special parse for file "99MVNSAS"
filename <- list.files(maltq_dir, pattern = "99MVNSAS", full.names = TRUE) %>%
  str_subset(., "~", negate = TRUE)

## Read it in - only the second sheet
col_names <- c("Location", "Entry_number", head(names(maltq_data_out$WRSBN99.xlsx$maltq_data)[-1], -5))

data_99MVNSAS <- read_excel(path = filename, sheet = "4MVN", col_names = col_names) %>%
  mutate_at(.vars = setdiff(names(select_if(., is.character)), c("Location", "VarietyorSelection")), parse_number) %>%
  mutate(Year = "1999", Nursery = "MVN") # Add year and nursery

## Add this to the list
maltq_data_out[[basename(filename)]] <- list(maltq_data = data_99MVNSAS, entry_info = NULL)






## Tidy the quality data
maltq_data_tidy <- maltq_data_out %>%
  map_df("maltq_data") %>%
  # Some renaming
  rename(line_name = VarietyorSelection, location = Location, row_type = Rowed, nursery = Nursery, year = Year) %>%
  select(matches("^[a-z]", ignore.case = FALSE), names(.))


## Consolidate traits - only include important traits (no rank or score)
traits <- c("kernelweight" = "KernelWeight", "on664" = "PlumpGrain", "barleycolor" = "BarleyColor", 
            "maltextract" = "MaltExtract", "st" = "ST", "dp" = "DiastaticPower", "wortcolor" = "WortColor",
            "turb" = "Turbidity", "wortclarity" = "WortClarity", "barleyprotein" = "GrainProtein", "wortprotein" = "WortProtein",
            "betaglucan" = "BetaGlucan", "freeaminonitrogen" = "FreeAminoNitrogen", "fan" = "FreeAminoNitrogen", 
            "visc" = "Viscosity", "pH" = "WortpH")

# List to store results
tr_list <- vector("list", length(traits))



# Combine information by trait
for (tr in unique(traits)) {
  # Get the trait name pattern
  tr_pattern <- names(subset(traits, traits %in% tr))
  
  # Subset data that matches the name
  tr_data <- bind_cols(map(tr_pattern, ~select(maltq_data_tidy, matches(.))))
  
  
  ## Consolidate
  tr_list[[tr]] <- tr_data %>%
    # Use coalesce to find the first non-NA instance in a set of vectors
    transmute(temp = coalesce(!!!.)) %>%
    rename_all(~tr)
  
}


# State names and Canadian provinces to remove
province_df <- rvest::html_table(rvest::html_session("https://www12.statcan.gc.ca/census-recensement/2011/ref/dict/table-tableau/table-tableau-8-eng.cfm"))[[1]][-1,]
province.names <- province_df$`Province/Territory`
province.abb <- province_df$`Internationally approved alpha code (Source: Canada Post)`


## Write the unique location names to a csv
maltq_data_tidy$location %>% unique() %>% writeClipboard()

## Read in the cleaned location df
locations_parsed <- read_excel(path = file.path(extr_dir, "maltq_nursery_locations.xlsx"))




## Combine with the original data
maltq_data_tidy1 <- maltq_data_tidy %>%
  ## Merge with parsed location data
  left_join(., select(locations_parsed, original, town), by = c("location" = "original")) %>%
  mutate(location = town,
         # Create environment variable
         environment = paste0(toupper(abbreviate(location, 3)), str_sub(year, 3, 4))) %>%
  select(matches("^[a-z]", ignore.case = FALSE), -on664, -pH) %>%
  bind_cols(., bind_cols(tr_list)) %>%
  # Gather traits
  gather(trait, value, matches("^[A-Z]", ignore.case = FALSE)) %>%
  # Filter out missing values and line names
  filter(!is.na(value), !is.na(line_name)) %>%
  arrange(year, nursery, location)


## Save this data
save("maltq_data_out", "maltq_data_tidy1", file = file.path(extr_dir, "maltq_data_extraction.RData"))

