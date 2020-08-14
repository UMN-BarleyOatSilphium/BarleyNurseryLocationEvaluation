## Code to parse the pedigrees from the regional nursery datasets
## 
## Author: Jeff Neyhart
## Date modified: 13 November 2019
## 
## 


# Packages to use
packages <- c("tidyverse", "readxl", "lubridate", "reticulate")
invisible(lapply(packages, library, character.only = TRUE))

## Source a script for relevant functions
source(file.path(getwd(), "functions.R"))

## Use my package for adding spaces to strings
library(neyhart)


# Directories containing data
shared_drive_dir <- c(paste0(LETTERS, ":/AGRO/BARLEY_LAB"), paste0(LETTERS, ":/BARLEY_LAB")) %>%
  subset(., dir.exists(.)) %>%
  head(1)

proj_dir <- getwd()


nursery_dir <- file.path(shared_drive_dir, "/Breeding/BARLEY/HistoricalNurseryData/")
raw_dir <- file.path(nursery_dir, "RawData")
agro_dir <- file.path(raw_dir, "Agronomic")
maltq_dir <- file.path(raw_dir, "Quality")
tidy_dir <- file.path(nursery_dir, "TidyData")
extr_dir <- file.path(nursery_dir, "ExtractedData")

# Project-specific directories
data_dir <- file.path(proj_dir, "Data")

# Grade priority
grade_priority <- c("MALT", "FEED", "FOOD")



# Grab pedigrees and lines from the extracted nursery data lists ----------



# 1. Assemble metadata
# 2. Extract pedigrees from metadata
# 3. Extract all line names from phenotype data

# Load the extracted nursery data lists
load(file.path(extr_dir, "maltq_data_extraction.RData"))



### Extract and assemble entry metadata information ###

# Extract entry metadata from the maltq data
maltq_entry_metadata <- maltq_data_out %>% 
  map("entry_info") %>% 
  subset(., map_lgl(., ~!is.null(.))) %>%
  imap_dfr(~mutate(.x, report = str_remove(.y, "\\..*$"))) %>%
  select(-number)
  
# Extract entry metadata from MVN
mvn_entry_metadata <- read_csv(file = file.path(extr_dir, "mvn_extracted_entry_metadata.csv"))

# Extract entry metadata from WRN
wrn_entry_metadata <- read_csv(file = file.path(extr_dir, "wrn_extracted_entry_metadata.csv")) %>%
  select(-entry) %>%
  rename(entry = name)


## Combine and make some edits
all_nursery_entry_metadata <- bind_rows(maltq_entry_metadata, mvn_entry_metadata, wrn_entry_metadata) %>%
  # Coalesce columns
  mutate(line_name = coalesce(line_name, entry),
         row_type = coalesce(row_type, type),
         end_use = coalesce(end_use, grade)) %>%
  select(line_name, program, pedigree, alias, row_type, end_use, report) %>%
  filter(!is.na(line_name))


## Mutate
## 1. Everything to capitals, extract aliases, remove spaces
## 2. row_type to numeric
## 3. End use to malt/feed/food
all_nursery_entry_metadata1 <- all_nursery_entry_metadata %>%
  mutate(row_type = parse_number(row_type),
         end_use = case_when(
           str_detect(toupper(end_use), "M|MALT") ~ "MALT",
           str_detect(toupper(end_use), "FOOD") ~ "FOOD",
           str_detect(toupper(end_use), "F|FEED") ~ "FEED",
           TRUE ~ as.character(NA)
         )) %>%
  mutate_at(vars(line_name, program, pedigree, alias), ~str_to_upper(.) %>% str_remove_all(" ")) %>%
  # Extract aliases
  mutate(alias1 = str_extract(line_name, "\\(.*\\)") %>% str_remove_all("\\(|\\)"),
         line_name = str_remove_all(line_name, "\\(.*\\)") %>% str_trim(),
         aliases = map2(alias, alias1, ~na.omit(c(.x, .y)) %>% as.character)) %>%
  select(-alias, -alias1)


## Write to CSV; make manual edits to things like pedigrees as line names or weird names
all_nursery_entry_metadata1 %>%
  mutate(aliases = map_chr(aliases, ~paste0(., collapse = ",")),
         .id = seq(nrow(.)),
         line_nameCorrected = line_name, pedigreeCorrected = pedigree) %>%
  arrange(line_name, report) %>%
  write_csv(x = ., path = file.path(extr_dir, "nursery_entry_metadata_to_correct.csv"), na = "")




# New aliases based on this pedigree information --------------------------



## New aliases
# Manually curate new aliases
new_aliases_from_pedigree <- tribble(
  ~line_name, ~pedigree, ~aliases,
  "ST", NA, "GD-C2",
  "B8907758", "CLIVIA/9448-83", c(),
  "D1-72", "SHYRI/GALENA", c(),
  "ANT28-484", NA, c("GRIT"),
  "CAMINANT", "ANT28-484/BLENHEIM", c(),
  "ANT499", NA,  c("APEX"),
  "12697-94", NA, c("ANT643"),
  "ANT29-667",  NA, c("HARRINGTON"),
  "WA7758-89", "CLIVIA/9448-83", c(),
  "BARONESSE", "BUMPER/HAZEN//AZURE", c("WA7642-92"),
  "EBC-187", NA, c("WA7114-93"),
  "ND12200", NA, c(),
  "M141", NA, c("FEG175-57"),
  "QUEST", NA, c("FEG65-02", "M122"),
  "X06G10", "BARONESSE/HARRINGTON", c(),
  "05WA-344.1", "PMUT-422H/CDCCANDLE", c(),
  "LOGANSIB", "ND7085/BOWMANSIB//ND7556", c(),
  "LOGAN", "ND7085/BOWMANSIB//ND7556", c(),
  "WA7642-92", "7190-86/BARONESSE", c(),
  "WA7114-93", "9029-84/EBC-187", c(),
  "FEG18-20", "MNBRITE/SI4-29", c(),
  "MX", "M77-840/M79-20", c(),
  "MOREX", "M41/M79-20", c(),
  "GD2-18", "HARRINGTON/EXCEL", c(),
  "ELLICE", NA, c("PI503880"),
  "VA72-44-362", NA, c("HENRY", "72-44-362"),
  "WA8537-68", "WA7698-62/FOMA", c(),
  "WA7698-62", "BETZES/HEINESHANNA/PIROLINE", c(),
  "MN99-101", NA, c("M120"),
  "FEG73-13", NA, c("M123"),
  "MN00-61", NA, c("M121"),
  "WA2196-68", "LUTHER/HUDSON", c(),
  "WA2509-65", "ALPINE/SVALOF//WHITEWINTER/TRIPLEBEARDEDMARIOUT-305", c(),
  "92-42-46", NA, "LR-RES",
  "AIM", NA, c("RPH3"),
  "96-44-307", "CALLAO/89-42-10", c(),
  "89-42-10", "CI4979/MONROE", c(),
  "97B-283", "CALLAO/SC830366", c(),
  "SC830366", "MCN601/HRSN/2/GEMBLOUX/3/72-44-362", c(),
  "VA97B-415", "90-42-56/CALLAO//92-42-62", c(),
  "90-42-56", "BSY*2//A/A/3/79-44-167", c(),
  "92-42-62", "BSY//B/A", c(),
  "VA97B-233", "90-41-9/88-44-725", c(),
  "90-41-9", "JFS/A//BSY/RAPIDAN", c(),
  "88-44-725", "RABAT/MONROE//79-44-167", c(),
  "M00-62", NA, c("M125"),
  "WA6415-66", "WA3564/UNITAN", c(),
  "UTSDB1-1009", "WOODVALE // PRIMUS / S.D. 67-297", c("UT-S.D.B1-1009"),
  "UT81B306-1731", "STEPTOE/M27", c(),
  "WALKER", "STEPTOE/M27", c(),
  "BRACKEN", "WA3564/UNITAN", c(),
  "WA641566", "WA3564/UNITAN", c(),
  "UT75B65-532", "ID633019/WOODVALE", c(),
  "CI3515", "HARRISON/3/CEBADACAPA/WONG//AWNLETEDHUDSON", c("A"),
  "JEFFERSON", "HARRISON/3/CEBADACAPA/WONG//AWNLETEDHUDSON", c("A"),
  "THOROUGHBRED", "VA90-44-110/SC872143", c(),
  "VA90-44-110", "CI8618/SRY//SUX/3/SRY", c(),
  "SC872143", "VA75-42-45/SC793556//CI2457", c("PLAISANT"),
  "96-41-39", "NOMINI/85-41-3", c(),
  "VA75-42-45", NA, c("75-42-45"),
  "96B-315", "89-42-8/CALLAO", c(),
  "89-42-8", "CI4979/MONROE", c(),
  "SC871077", "75-42-45/SC793556//CI2457", c(),
  "90-42-56", "BSY*2//A/A/3/79-44-167", c(),
  "H-585", "VA75-42-45/SC793556//CI2457", c(),
  "VA97B-176", "89-42-8/CALLAO", c(),
  "VA98B-112", "89-42-8/VA99B-172", c(),
  "VA99B-172", "CALLAO//90-41-10/90-42-9", c(),
  "VA97B-398", "CALLAO/92-42-46", c(),
  "VA97B-275", "CALLAO/SC83036", c(),
  "VA95-42-58", "WYSOR/88-44-301", c(),
  "VA96-44-304", "CALLAO/89-41-6", c(),
  "VA97B-284", "CALLAO/SC830366", c(),
  "VA97B-178", "89-42-8/CALLAO", c(),
  "SC911490", "VA84-44-78/HOV/CC/WG/HSN/HRS/SC812412", c(),
  "VA96B-113", "PENNBAR66/83-44-59", c(),
  "VA96B-70", "CALLAO/PAMUNKEY", c(),
  "VA00B-214", "CMB82A-520/91-44-611//PAMUNKEY/3/CALLAO", c(),
  "CMB82A-520", "GLORIA_S_/SAIDA", c(),
  "VA96B-301", "CALLAO/89-41-6", c(),
  "VA00B-259", "CMB74A-333//90-44-90/90-42-22/3/92-42-46", c(),
  "CMB74A-333", "MANKER/ATHENAIS", c(),
  "VA96-44-321", "CALLAO/SC830366", c(),
  "VA99B-303", "CMB79-54/CALLAO//PAMUNKEY", c(),
  "H159-005", "H30-52/H30-11", c(),
  "H30-52", NA, c("MCGREGOR"),
  "H30-11", NA, c("MCDIARMID"),
  "VA99B-10", "85-41-3/CALLAO", c(),
  "VA00B-244", "90-42-45/90-44-1//CR366.13.2/3/92-42-46/4/CALLAO", c(),
  "90-44-1", "BOONE/HRY//79-44-167", c(),
  "VA99B-319", "CMB80A-56//90-42-47/CALLAO/3/PAMUNKEY", c(),
  "98306-6", NA, c("CI5098", "PI83"),
  "VA97B-398", "CALLAO/92-42-46", c(),
  "BARONESSE", NA, c("WPB-BZ-594-35"),
  "BZ585-85", "WAXBAR/TR-451", c(),
  "WA8611-90", "SM8618/16277-85", c(),
  "29-53", "LK6-44/EXCEL", c(),
  "75-15", "LK6-44/CONLON", c(),
  "75-26", "LK6-44/CONLON", c(),
  "75-35", "LK6-44/CONLON", c(),
  "FEG141-20", NA, c("M135"),
  "X04041-T32", "01WA-13862.3/RADIANT", c(),
  "X04041-T34", "01WA-13862.3/RADIANT", c(),
  "X04041-T86", "01WA-13862.3/RADIANT", c(),
  "X04055-T21", "01WA-12501.2/RADIANT", c()
)



## Read in the corrected data - combine with new aliases
all_nursery_entry_metadata2 <- read_excel(path = file.path(extr_dir, "nursery_entry_metadata_corrected_20200213.xlsx")) %>%
  filter(!is.na(line_nameCorrected)) %>%
  mutate_at(vars(contains("Corrected")), ~str_to_upper(.) %>% str_remove_all(" ")) %>%
  select(-line_name, -pedigree) %>%
  rename_at(vars(contains("Corrected")), ~str_remove(., "Corrected")) %>%
  bind_rows(., mutate(new_aliases_from_pedigree, aliases = map_chr(aliases, ~paste0(., collapse = ",")) %>% ifelse(. == "", NA, .)))



# Examine line names from the nursery trait data --------------------------


        

### Extract and assemble information from ALL phenotyped lines that are not in the pedigrees ###

## First load the phenotype data; find distinct line names
line_names_from_trait_data <- read_csv(file = file.path(tidy_dir, "nursery_trait_data_tidy.csv")) %>%
  distinct(line_name) %>%
  # Same corrections as above
  mutate(line_name = str_to_upper(line_name) %>% str_remove_all(" "))

# Filter out line names that are already in the pedigree information df (line names or aliases)
line_names_from_trait_data_exlcluded <- line_names_from_trait_data %>%
  filter(! line_name %in% all_nursery_entry_metadata2$line_name) %>%
  filter(! map_lgl(line_name, ~any(str_detect(string = all_nursery_entry_metadata2$aliases, pattern = .), na.rm = TRUE)))

## Add these lines to the metadata df
all_nursery_entry_metadata3 <- all_nursery_entry_metadata2 %>%
  bind_rows(., line_names_from_trait_data_exlcluded)






# Load reference information ----------------------------------------------




## Load the reference information from the shared drive
t3_pedigree_reference <- read_csv(file = file.path(shared_drive_dir, "Breeding/BARLEY/T3 uploading/T3LineRecords/t3_pedigree_reference.csv")) %>%
  mutate(aliases = map(aliases, ~str_split(., ", ")[[1]]),
         aliases = modify_if(aliases, ~all(is.na(.)), ~character(0)))





# Load and edit nursery line name information -----------------------------





## Combine agro and maltq entry info
## Make corrections:
## 1. Capitalize
## 2. Remove spaces
## 3. Extract aliases in parens
## 4. Make aliases a list
nursery_entry_info1 <- all_nursery_entry_metadata3 %>%
  select(-.id) %>%
  # Parse
  mutate_if(is.character, parse_guess) %>%
  # Capitalize, remove spaces
  mutate_if(is.character, ~str_remove_all(., " ") %>% str_remove_all(., "_") %>% str_trim() %>% toupper()) %>%
  mutate(line_name = str_replace_all(line_name, "1D", "ID")) %>%
  # Capitalize, remove spaces - do it again
  mutate_if(is.character, ~str_remove_all(., " ") %>% str_trim() %>% toupper()) %>%
  arrange(line_name)


## Some lines have duplicate names and pedigree, but non-matching row_types.
## Fix these
nursery_entry_info2 <- nursery_entry_info1 %>%
  mutate(row_type = case_when(
    line_name == "STEPTOE" ~ 2,
    line_name == "WA8611-90" ~ 2,
    line_name == "WA8625-90" ~ 2,
    line_name == "WA9792-90" ~ 2,
    TRUE ~ as.numeric(row_type)
  )) %>%
  # Correct some names in pedigrees
  mutate(pedigree = str_replace_all(string = pedigree, pattern = "MTT", "MT"),
         pedigree = ifelse(pedigree == "BARONESSE/3/CRYSTAL/KLAGES*3/PI366450", "BARONESSE/3/CRYSTAL//KLAGES*3/PI366450", 
                           pedigree))





## Examimne the pedigrees and pedigree aliases 
## to extract additional pedigrees and aliases
examine_pedigree <- nursery_entry_info2$pedigree
pedigree_repaired <- str_trim(examine_pedigree)

# Look for parentheses, brackets or braces
bracket_detect_which <-  str_which(examine_pedigree, "\\(|\\[|\\{")
bracket_detect <- examine_pedigree[bracket_detect_which]

unique(bracket_detect)



## Remove brackets, braces, and parentheses, and anything inside them
bracket_fix <- str_replace_all(bracket_detect, pattern = "\\(.*\\)", replacement = "") %>%
  str_replace_all(pattern = "\\[.*\\]", replacement = "") %>%
  str_replace_all(pattern = "\\{.*\\}", replacement = "")

## Add repaired pedigrees
pedigree_repaired[bracket_detect_which] <- bracket_fix



# Look for punctuation
punct_detect_which <-  str_which(pedigree_repaired, "\\!|\\#|\\$|\\%|\\&\\'|\\+|\\,|\\'|\\:|\\;|\\<|\\=|\\>|\\?|\\@|\\^|\\_|\\`")
punct_detect <- pedigree_repaired[punct_detect_which]

## Remove punctuation
punct_fix <- str_remove(punct_detect, "=STEPTOE/LARKER") %>%
  str_remove(., ",HVM74=KK") %>%
  str_remove_all(., ",F[0-9]") %>%
  str_replace_all(., "\\#", "_") %>%
  str_replace_all(., "\\+", "-") %>%
  str_remove_all(., "\\,|\\@|\\'") %>% # Remove apostrophes
  str_remove_all("HASEXCELROBUSTANDM47INTHE") %>%
  str_remove_all("\\[89-42-8|\\[CMB80A-56|\\[VA90-44-110:SEGHULLESS-96|\\[PENNBAR66")
punct_fix[punct_fix == ""] <- NA

## Add repaired pedigrees
pedigree_repaired[punct_detect_which] <- punct_fix





## Add fixed pedigrees back into the data.frame
nursery_entry_info3 <- nursery_entry_info2 %>%
  mutate(pedigree = pedigree_repaired)



# Find unique line_names.
# Then build up pedigree and accessory information
unique_line_name <- nursery_entry_info3 %>% 
  # Drop report
  select(-report) %>%
  group_by(line_name) %>%
  nest()

## Merge duplicate information per line name
unique_line_name1 <- unique_line_name %>%
  group_by(line_name) %>%
  do(merge_entry_info(x = .$data[[1]])) %>%
  ungroup()

## Examine remaining duplicates
unique_line_name1 %>% 
  group_by(line_name) %>%
  filter(n() > 1) %>%
  as.data.frame()


## allow those duplicates to proceed; pedtools will catch them.

# Rename the nursery information;
# rename the reference information too
nursery_entry_data_referenced <- unique_line_name1 %>%
  mutate(original = line_name, .id = seq(nrow(.)),
         # Fix /2/ in a pedigree
         pedigree = ifelse(str_detect(pedigree, "/2/"), str_replace_all(pedigree, "/2/", "//"), pedigree))




# Expand pedigrees --------------------------------------------------------


# Prepare the reference df
reference <- t3_pedigree_reference %>%
  filter(!is.na(pedigree)) %>%
  # Remove entries from composite crosses or 2 row male sterile pop
  filter(str_detect(pedigree, "CROSS|POPULATION", negate = TRUE)) %>%
  ## Add aliases as lines; this will be fixed in pedTools
  bind_rows(., unnest(filter(., !map_lgl(aliases, is_empty))) %>% select(line_name = aliases, pedigree) ) %>%
  select(-aliases)



## Expand the pedigrees
nursery_entry_data_expanded_list <- nursery_entry_data_referenced %>%
  group_by(line_name, .id) %>%
  nest() %>%
  mutate(out = list(NULL))

i <- min(which(map_lgl(nursery_entry_data_expanded_list$out, is.null)))

# Loop over individuals
for (i in seq(i, nrow(nursery_entry_data_expanded_list))) {

  print(i)
  
  input <- nursery_entry_data_expanded_list[i,]
  data <- input$data[[1]]
  epx <- expand_pedigree2(line_name = input$line_name, pedigree = data$pedigree, reference = reference)
  
  ## Add output to list
  nursery_entry_data_expanded_list$out[[i]] <- epx
  
}


# ## Debug
# line_name <- input$line_name
# pedigree <- data$pedigree
# # reference_temp <- reference
# reference <- reference1



## Unnest
nursery_entry_data_expanded <- nursery_entry_data_expanded_list %>% 
  unnest(data, out) %>%
  select(line_name, pedigree, original, aliases, row_type, end_use, accessory) %>%
  # Replace EXCEL* in pedigree
  mutate(pedigree = ifelse(pedigree == "EXCEL*/M60", "EXCEL/M60", pedigree)) %>%
  filter(!is.na(line_name))


# Find the unique accessory individuals and pull from the reference
accessory_reference_df <- reference %>%
  filter(line_name %in% reduce(nursery_entry_data_expanded$accessory, union))

# Pull these same lines from the entry database - search for both line names and 
# aliases
reference_data_use <- t3_pedigree_reference %>% 
  filter(line_name %in% accessory_reference_df$line_name | 
           map_lgl(aliases, ~any(. %in% accessory_reference_df$line_name) )) %>%
  select(line_name, pedigree, aliases)




## Save as RData - this will be used later to change names in the dataset
save("nursery_entry_data_expanded", "reference_data_use", 
     file = file.path(extr_dir, "nursery_raw_entry_metadata.RData"))




# Prepare data for PedTools -----------------------------------------------



## Combine nursery data and reference data
entry_info_pedtools_input <- bind_rows(
  reference_data_use,
  select(nursery_entry_data_expanded, line_name, pedigree, aliases)
)

## Make some changes that hamper success of pedtools
# Remove  WA7642-92 as alias of Baronesse
entry_info_pedtools_input1 <- entry_info_pedtools_input %>%
  mutate(aliases = modify_if(.x = aliases, .p = entry_info_pedtools_input$line_name == "BARONESSE", 
                             .f = ~setdiff(.x, "WA7642-92"))) %>%
  # Change ROBUST2* to ROBUST*2 using the regex: ([A-Z0-9]*)([0-9]{1,})(\\*) to \\1\\3\\2
  mutate_at(vars(line_name, pedigree), ~str_replace_all(., "ROBUST2\\*", "ROBUST\\*2")) %>%
  # Expand list of aliases into a string
  mutate(aliases = map_chr(aliases, ~paste0(., collapse = ", ")),
         aliases = parse_character(aliases)) %>%
  distinct() %>%
  rowid_to_column(var = ".id")
  

## Filter or correct failed pedigrees
entry_info_pedtools_input1_pedfail <- entry_info_pedtools_input1 %>%
  filter(!map_lgl(pedigree, check_pedigree))

## Correct some, remove others
# Make some corrections
entry_info_pedtools_input1_pedfail1 <- entry_info_pedtools_input1_pedfail %>%
  mutate(pedigree = case_when(
    pedigree == "00NZ304XCELLAR" ~ "00NZ304/CELLAR",
    pedigree == "00NZ304X85AB2323" ~ "00NZ304/85AB2323",
    pedigree == "UT91B706-A-259XBU585-82" ~ "UT91B706-A-259/BU585-82",
    pedigree == "UT91B706-A-259XDA587-170" ~ "UT91B706-A-259/DA587-170",
    pedigree == "PB196-2R-6123XPB197-2R-7090" ~ "PB196-2R-6123/PB197-2R-7090",
    # If pedigree contains no slashes, but aliases are present, just set to NA
    str_count(pedigree, "/") == 0 & !is.na(aliases) ~ as.character(NA),
    TRUE ~ pedigree
  )) %>%
  # Remove those without slashes
  filter(str_count(pedigree, "/") > 0 | is.na(pedigree))
  


## Merge lines with identical pedigrees
entry_info_pedtools_input1_dups <- entry_info_pedtools_input1 %>%
  # First filter out .ids removed in the previous step
  filter(!.id %in% setdiff(entry_info_pedtools_input1_pedfail$.id, entry_info_pedtools_input1_pedfail1$.id)) %>%
  group_by(line_name) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  arrange(line_name) %>%
  # Merge aliases by line_name/pedigree
  group_by(line_name, pedigree) %>%
  do(aliases = paste0(as.character(na.omit(.$aliases)), collapse = ", ")) %>%
  unnest(aliases) %>%
  mutate(aliases = ifelse(aliases == "", NA, aliases)) %>%
  # Merge aliases if a duplicate line name has an NA pedigree
  group_by(line_name) %>%
  do({
    df <- .
    if (any(is.na(df$pedigree))) {
      line_name <- unique(df$line_name)
      pedigree <- as.character(na.omit(unique(df$pedigree)))
      aliases <- str_split(string = df$aliases, pattern = ",") %>% map(str_trim) %>% unlist() %>% unique()
      
      df1 <- crossing(line_name, pedigree, aliases = paste0(aliases, collapse = ", "))
    } else {
      df1 <- df
    }
    
    # Filter out failed pedigrees; if all failed, pass
    good_ped <- map_lgl(df1$pedigree, check_pedigree)
    df2 <- if (any(good_ped)) filter(df1, good_ped) else df1
    
    filter(df2, nchar(pedigree) == max(nchar(pedigree)))
    
  })


# Remove dups and poor .ids from the previous df; combined with edited dups
entry_info_pedtools_input2 <- entry_info_pedtools_input1 %>%
  # First filter out .ids removed in the previous step
  filter(!.id %in% setdiff(entry_info_pedtools_input1_pedfail$.id, entry_info_pedtools_input1_pedfail1$.id)) %>%
  group_by(line_name) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  select(-.id) %>%
  bind_rows(., entry_info_pedtools_input1_dups) %>%
  arrange(line_name) %>%
  ## Remove the string "NA" from alias strings
  mutate(aliases = str_split(aliases, ",") %>% map(str_trim) %>% map(~as.character(na.omit(parse_character(.)))),
         aliases = map_chr(aliases, ~paste0(., collapse = ", ")),
         aliases = ifelse(aliases == "", NA, aliases))
 

write_tsv(x = entry_info_pedtools_input2, 
          path = file.path(proj_dir, "Scripts/pedtools/barleyNurseryPedigree/nursery_line_info_pedtools_input.txt"), na = "")




# Format for calculating coefficients of coancestry -----------------------




# Load pedigreeTools
library(pedigreeTools)


# Read in the parent extended pedigree information
# Do not read in the version with expanded backcrosses
# pedtools_out <- read_tsv(file = file.path(data_dir, "nursery_line_info_pedigree_parents.txt"))
pedtools_out <- read_tsv(file = file.path(tidy_dir, "nursery_line_info_pedigree_bc_parents_step12.txt"))


## Determine parent contribution
pedtools_out1 <- pedtools_out %>%
  mutate_at(vars(contains("Parent")), list(contribution = ~str_extract(., "\\*[0-9]{1,}$|^[0-9]{1,}\\*") %>% parse_number())) %>%
  # Edit parent names
  mutate_at(vars(contains("contribution")), ~ifelse(is.na(.), 1, .)) %>%
  mutate_at(vars(contains("contribution")), ~(1 - (0.5^.))) %>%
  mutate(Parent1_contribution = map2_dbl(Parent1_contribution, Parent2_contribution, ~ifelse(.x < .y, 1 - .y, .x)),
         Parent2_contribution = map2_dbl(Parent2_contribution, Parent1_contribution, ~ifelse(.x < .y, 1 - .y, .x)),
         # Add a vector of selfing; if the entry name is pedc1, selfing == 0
         selfing = ifelse(str_detect(Name, "pedc"), 0, 10))
         

# ## Sanity check
pedtools_out1 %>% filter(Parent1_contribution + Parent2_contribution != 1)



# # Remove UNKNOWN as a entry name
# pedtools_out2 <- pedtools_out1 %>%
#   filter(Name != "UNKNOWN")
# 
# ## Replace UNKNOWN with pedx numbers
# # First find the highest pedx number
# pedx_max1 <- select(pedtools_out2, Name, Parent1, Parent2) %>% unlist() %>% 
#   str_subset(., "pedx") %>% parse_number() %>% max()
# 
# # Replace UNKNOWN in parent1
# which_unknown <- which(pedtools_out2$Parent1 == "UNKNOWN")
# pedtools_out2$Parent1[which_unknown] <- paste0("pedx", pedx_max + seq_along(which_unknown))
# 
# # Recalculate max
# pedx_max2 <- select(pedtools_out2, Name, Parent1, Parent2) %>% unlist() %>% 
#   str_subset(., "pedx") %>% parse_number() %>% max()
# 
# # Replace UNK in parent2
# which_unknown <- which(pedtools_out2$Parent2 == "UNKNOWN")
# pedtools_out2$Parent2[which_unknown] <- paste0("pedx", pedx_max + seq_along(which_unknown))
# 
# pedx_max2 <- select(pedtools_out2, Name, Parent1, Parent2) %>% unlist() %>% 
#   str_subset(., "pedx") %>% parse_number() %>% max()
# 
# ## Add these new pedx lines to the df
# pedtools_out3 <- pedtools_out2 %>%
#   add_row(Name = paste0("pedx", seq(pedx_max1 + 1, pedx_max2))) %>%
#   # Edit original parent names
#   mutate_at(vars(Parent1, Parent2), list(edit = ~str_remove(., "\\*[0-9]{1,}$|^[0-9]{1,}\\*")))
#   


## Add the new pedx lines to the df
pedtools_out3 <- pedtools_out1 %>%
  # Edit original parent names
  mutate_at(vars(Parent1, Parent2), list(edit = ~str_remove(., "\\*[0-9]{1,}$|^[0-9]{1,}\\*")))


## Attempt to sort the pedigree
# First convert for reading in pedigreeTools
pedigreeToolsIn <- pedtools_out3 %>%
  select(label = Name, sire = Parent1_edit, dam = Parent2_edit) %>%
  as.data.frame()

# Sort
pedigreeToolsIn_sorted <- editPed(sire = pedigreeToolsIn$sire, dam = pedigreeToolsIn$dam,
                                  label = pedigreeToolsIn$label)

# Merge the sorted table with the original table
pedigree_table <- pedigreeToolsIn_sorted %>% 
  select(name = label, par1 = sire, par2 = dam, generation = gene) %>% 
  left_join(., pedtools_out3, by = c("name" = "Name", "par1" = "Parent1_edit", "par2" = "Parent2_edit"))


## Calculate A assuming selfing
pedigreeA_self <- getASelfing(ID = pedigree_table$name, Par1 = pedigree_table$par1, Par2 = pedigree_table$par2, 
                              nCycles = pedigree_table$selfing, nCyclesDefault = 10)

# Convert to matrix
pedigreeA_self_mat <- as.matrix(pedigreeA_self)


## Check a known backcross entry
pedtools_out3 %>% 
  # filter(str_detect(Pedigree, "\\*4")) %>% 
  filter(str_detect(Pedigree, "\\*2")) %>% 
  as.data.frame()
# *2 backcross
pedigreeA_self_mat["00ID1550", c("COLTER", "PMUT-422H")] / 2
# *4 backcross
pedigreeA_self_mat["01ST1587", c("WPB-BZ-594-35", "STARS9301B")] / 2


## Looks good!




## Convert the matrix to a df
pedigree_Amat_df <- pedigreeA_self_mat %>%
  as.data.frame() %>%
  rownames_to_column("line1") %>%
  gather(line2, relat, -line1) %>%
  # Remove lines that are pedx or pedc
  filter(str_detect(line1, "pedx|pedc", negate = TRUE)) %>%
  filter(str_detect(line2, "pedx|pedc", negate = TRUE)) %>%
  filter(str_detect(line1, "/", negate = TRUE)) %>%
  filter(str_detect(line2, "/", negate = TRUE)) %>%
  as_tibble()

## Convert back to a matrix
pedigree_Amat <- pedigree_Amat_df %>%
  as.data.frame() %>%
  spread(line2, relat) %>%
  column_to_rownames("line1") %>%
  as.matrix() %>%
  Matrix::Matrix(., sparse = TRUE)



## Some statistics:
## 1. How many lines have zero relationship to all other lines?
pedigree_Amat_temp <- pedigree_Amat
diag(pedigree_Amat_temp) <- 0 # Set zero

# summarize
unrelated <- apply(X = pedigree_Amat_temp, MARGIN = 1, FUN = function(line) all(line == 0))
sum(unrelated)

# 197

names(subset(unrelated, unrelated))

# Distribution of relationships




# Save both
save("pedigree_Amat", "pedigree_Amat_df", "pedigree_table", file = file.path(tidy_dir, "nursery_pedigree_relmat.RData"))





  
