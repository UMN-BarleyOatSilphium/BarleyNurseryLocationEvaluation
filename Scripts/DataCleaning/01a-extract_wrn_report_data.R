## Code to extract metadata from the WRN reports
## 
## Jeff Neyhart
## 

# Packages to use
packages <- c("data.table", "tidyverse", "lubridate", "readxl", "pdftools")
invisible(lapply(packages, library, character.only = TRUE))


# Trait metadata
trait_metadata <- read_csv("C:/GoogleDrive/BarleyLab/Breeding/PhenotypicData/Metadata/trait_metadata.csv")


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

## Source a script for relevant functions
source(file.path(getwd(), "source_functions.R"))


# Change working directory to BARLEY LAB directory
wrn_dir <- file.path(nursery_dir, "PDF_Reports/WesternBarleyNursery")
report_dir <- wrn_dir


### Tools for data extraction
# Renaming vector for pedigree columns
ped_col_rename <- c("Source" = "Program", "No" = "Number", "Number" = "Number", "Entry" = "Entry", 
                    "Parentage" = "Pedigree", "Type|TYPE" = "Type", "Grade" = "Grade")




## List all report files
pdf_files <- sort(list.files(report_dir, pattern = "WRBN", full.names = TRUE)) %>%
  str_subset(string = ., pattern = ".pdf")

# Create an empty list to store output
extracted_data_list <- list()


## Iterate over files
for (p in seq_along(pdf_files)) {
# for (p in seq(p, length(pdf_files))) {
  
  
  pdf_file <- pdf_files[p]
  
  ## Extract text
  pdf_txt <- pdf_text(pdf = pdf_file) %>%
    # convert everything to upper-case
    str_to_upper()
  
  # Get coordinates of text
  pdf_dat <- pdf_data(pdf = pdf_file)
  
  # Which is the contents page?
  toc_page <- which(str_detect(pdf_txt, "TABLE OF CONTENTS"))
  # Which are the general information pages?
  # Exclude pages up until the toc page
  gen_info_pages <- setdiff(which(str_detect(pdf_txt, "GENERAL INFORMATION")), seq(toc_page))

  
  ## Pattern matching to extract relevant items in the table of contents
  pattern <- "TABLE [0-9]*\\..*[0-9]{1,2}|GENERAL .*[0-9]{1,2}"

  ## Get the pages of all relevant content
  content_df <- str_extract_all(pdf_txt[[toc_page]], pattern = pattern)[[1]] %>%
    ## Split the table names and page numbers
    tibble(raw = .) %>%
    mutate(page_number = parse_number(str_extract(string = raw, pattern = "[0-9]{1,2}$")),
           table = str_remove_all(string = raw, pattern = "[0-9]{1,2}$") %>% str_trim(),
           pdf_page_number = page_number + toc_page) 
  
  # Pages to exclude in this search
  pages_to_omit <- union(seq(toc_page), gen_info_pages)
  # Make these pages empty for searching
  pdf_text_search <- pdf_txt
  pdf_text_search[pages_to_omit] <- replicate(length(pages_to_omit), "", simplify = FALSE)
  
  content_df <- content_df %>%
    ## Use a quick pattern match to determine if the TOC is correctly displaying page
    ## numbers. Use the most appropriate column
    filter(str_detect(table, "TABLE")) %>%
    mutate(tab_pattern = map_chr(table, ~paste0(str_extract(., "TABLE [0-9]*"), c("\\:", "\\."), collapse = "|")),
           correct_page = map_dbl(tab_pattern, ~str_which(string = pdf_text_search, pattern = .)[1])) %>%
    # If correct_page is NA, search using the trait name
    mutate(correct_page = pmap_dbl(select(., correct_page, table, tab_pattern), 
                               ~ifelse(is.na(..1), str_which(string = pdf_text_search, pattern = str_trim(str_remove(..2, ..3))), ..1)) ) %>%
    select(table, page_number, pdf_page_number, correct_page) %>%
    bind_rows(., subset(content_df, str_detect(table, "TABLE", negate = T), -raw)) %>%
    arrange(page_number)
  
  

  # which column is correct?
  page_column <- content_df %>% 
    summarize_at(vars(contains("number")), ~mean(. == correct_page, na.rm = TRUE)) %>% 
    which.max() %>% 
    names()

  
  
  #############################
  ## Year
  #############################
  
  report_year <- str_extract(pdf_txt[[2]], "[0-9]{4}")
  
  
  
  #############################
  ## Station Names
  #############################
  
  ## Stations are listed on the general information pages
  gen_info_pages <- subset(content_df, str_detect(table, "GENERAL"))[[page_column]]
  
  ## Iterate and collect station info and the nursery name (dryland or regular)
  station_info_list <- list()
  
  for (page in gen_info_pages) {
    
    page_text <- pdf_txt[[page]]
    
    # Extract the nursery name
    nursery_name <- str_extract(page_text, pattern = paste0(c(paste0(".*\\, ", report_year), paste0(report_year, ".*")), collapse = "|")) %>%
      str_trim() %>%
      str_remove(., paste0("\\, ", report_year))
    
    # Extract the locations/stations
    station_names <- page_text %>% 
      str_split(string = ., pattern = "GENERAL INFORMATION") %>% 
      .[[1]] %>% 
      head(1) %>%
      str_extract_all(string = ., "[A-Z ]{1,}\\, [A-Z]{2}") %>%
      .[[1]] %>%
      # Remove THE
      str_remove_all(., "THE") %>%
      str_trim() %>%
      subset(., str_count(., " ") <= 3) %>%
      unique()
    
    # Create location names
    station_names_df <- tibble(raw = station_names) %>%
      separate(data = ., col = raw, into = c("location", "state"), sep = ", ", remove = FALSE) %>%
      mutate(location = str_to_title(location))
      
    ## Add to the list
    station_info_list[[nursery_name]] <- station_names_df
    
  }
  
  
  ## Combine into a single DF
  station_info_df <- imap_dfr(station_info_list, ~mutate(.x, nursery = .y)) %>%
    filter(str_detect(raw, "DATA|AND", negate = TRUE))
  
  
  
  
  #############################
  ## Variety pedigree information
  #############################
  
  # Get the pages
  pedigree_pages <- subset(content_df, str_detect(table, "PLANTING LIST|ENTRY LIST"))[[page_column]]
  
  ## Iterate and collect pedigree information
  pedigree_info_list <- map(station_info_list, ~NULL)

  
  for (i in seq_along(pedigree_pages)) {
    
    page <- pedigree_pages[i]
    page_text <- pdf_txt[[page]]
    
    ## Get coordinates for this table
    page_dat <- pdf_dat[[page]]
    
    ## Find the column names
    col_coord <- page_dat %>% 
      filter(str_detect(text, paste0(names(ped_col_rename), collapse = "|"))) %>%
      # subset for the most common y
      filter(y %in% names(which.max(table(y))))
    

    
    # Determine suitable x column coordinates
    # Find text with x that falls between the x of the column names
    ## So, first find the ranges
    ymax <- subset(page_dat, str_detect(text, "[Nn]ew")) %>%
      subset(str_detect(tolower(text), "newdale", negate = TRUE)) %>% # Not newdale
      subset(str_detect(tolower(text), "^new|^\\*new")) %>% 
      pull(y) %>% 
      min()
    

    # Filter the rows beyond the columns
    page_dat_post_col <- page_dat %>% 
      # Move past the columns
      filter(y > unique(col_coord$y)) %>%
      # Stop at new entries
      filter(y < ymax)
    
    ## Since the data is more or less centered within a cell,
    ## we can use the average x coordinate of each column name (x + (x + width) / 2).
    ## Then for each column, we find data that overlaps with this average. We take the longest
    ## data string (or the column names if they are longer) and use these to demark the columns in
    ## the table
    
    ## Create a list to store the column cut coordinates
    col_data_list <- list()
    
    ## List over columns
    for (j in seq(nrow(col_coord))) {
      
      # Determine what text overlaps with that column
      # x1 is first column
      # x2 is the next column
      
      ## If j == 1, x1 is the min x
      if (j == 1) {
        x1 <- min(page_dat$x)
        
      } else {
        x1 <- col_coord$x[j]
        
      }
        
      # Width of the column name
      width1 <- col_coord$width[j]
      
      ## The text to target should be between(x1, x, x+width)
      dat_use <- page_dat_post_col %>% 
        filter(map2_lgl(x, width, ~any(seq(.x, .x + .y) %in% seq(x1, x1 + width1))))
      
      ## If dat_use is nrow 0, increase the margin of x1
      if (nrow(dat_use) == 0) {
        dat_use <- page_dat_post_col %>% 
          filter(map2_lgl(x, width, ~any(seq(.x, .x + .y) %in% seq(x1 - 35, x1 + width1))))
      }
      
      # Which is longer, the column name or the longest data string?
      # This is x2
      x2 <- pmax(x1 + width1, max(dat_use$width + dat_use$x))
      
      # Which is less, the left bound of the column name or the left bound of the 
      # most left element in the column
      # This is the new x1
      x1 <- pmin(x1, min(dat_use$x))
      
      
      ## Find data that falls within this range
      ## Add to the list
      col_data <- filter(page_dat_post_col, between(x, x1, x2)) %>%
        # split by y and concatenate
        split(.$y) %>%
        map_df(~tibble(y = unique(.x$y), text = paste0(.x$text, collapse = " "))) %>%
        ## Filter out rows that are just *
        filter(text != "*")
      
      
      col_data_list[[j]] <- col_data
      
    }
    
    # Filter the nrow = 0
    to_keep <- which(map_lgl(col_data_list, ~nrow(.) > 0))
    
    # Remove this from the data list and the column name list
    col_data_list <- col_data_list[to_keep]
    col_data_names <- col_coord$text[to_keep]
    
    ## Re-order the list so the longest column is first
    col_orders <- order(map_dbl(col_data_list, nrow), decreasing = TRUE)
    col_data_list <- col_data_list[col_orders]
    # Reorder the names, too
    col_data_names <- col_data_names[col_orders]
    
    
    ## Merge by y
    ## Use coordinate join with a tolerance for the y value
    col_data_df <- reduce(col_data_list, coord_join, by = "y", tol = 5) %>%
      select(-y) %>%
      # Add column names
      rename_all(~col_data_names) %>%
      # Remove punctuation
      rename_all(~str_remove_all(., "[:punct:]")) %>%
      # Rename
      rename_all(~str_replace_all(string = ., ped_col_rename)) %>% 
      ## Parse
      mutate_at(vars(which(names(.) != "Number")), parse_character)
    
    # If the number column is absent, add it
    if (!any(names(col_data_df) %in% "Number")) {
      col_data_df$Number <- as.character(seq(nrow(col_data_df)))
      
    }
    
    col_data_df <- col_data_df %>%
      mutate_at(vars(which(names(.) == "Number")), parse_number) %>%
      filter(!is.na(Number)) %>%
      arrange(Number)
    

    ## Add the pedigree to the list
    pedigree_info_list[[i]] <- col_data_df
    
  }
  

  ## Edit pedigree data
  ## Combine into a single DF
  pedigree_info_df <- pedigree_info_list %>%
    map(~{
      # Lowercase column names
      rename_all(.x, tolower) %>%
        mutate(entry = str_remove(entry, "\\*") %>% str_trim(),
               alias = str_extract(entry, "\\(.*\\)"),
               alias = str_remove_all(alias, "\\(|\\)"),
               entry = str_trim(str_remove(entry, "\\(.*\\)"))) %>%
        # Capitalize all
        mutate_if(is.character, toupper) %>%
        # Remove completely empty rows
        filter_all(any_vars(!. == ""))
    }) %>%
    imap_dfr(., ~mutate(.x, nursery = .y))
      
  
  
  
  
  #############################
  ## Extract trait data
  #############################
  
  # Trait pattern
  trait_pattern <- "GRAIN YIELD|TEST WEIGHT|PLANT HEIGHT|HEADING DATE|PLUMP BARLEY|THIN BARLEY"
  
  # Determine the tables to extract
  trait_content_df <- content_df %>%
    filter(str_detect(table, trait_pattern))
  
  
  # A list to store df
  trait_data_list <- list()
  
  ## Iterate over the content df
  for (i in seq(nrow(trait_content_df))) {

    # Get the page
    page <- trait_content_df[[page_column]][i]
    
    # Get the text and coordinates at this page
    page_text <- pdf_txt[[page]]
    page_dat <- pdf_dat[[page]]
    
    # What is the table name
    table_name <- str_extract(string = page_text, pattern = "TABLE [0-9]*.*")
      
    
    # ## Get the table - guess
    # tab <- extract_tables(file = pdf_file, pages = page)[[1]]
    
    ## Use the station information above to help determine columns
    station_col_names <- unique(station_info_df$location) %>% 
      subset(., . != "") %>%
      str_split(" ") %>% 
      map_chr(~.[which.max(length(.))]) %>% 
      paste0(collapse = "|")
    
    ## Look for relevant column names
    col_pattern <- paste0("Number|Entry|Ent|Selection|Cultivar|", station_col_names)
    
    ## Find the column names and coordinates
    col_coord <- page_dat %>%
      # Remove special characters
      filter(str_detect(str_to_title(text), col_pattern)) %>%
      filter(y < 200) %>%
      arrange(x)
      
    ## Choose number over entry, if both are present
    detect_col <- map(.x = c("Number", "Entry"), ~str_which(string = col_coord$text, pattern = .))
    if (length(unlist(detect_col)) > 1) {
      col_coord <- col_coord[-setdiff(unlist(detect_col), unlist(detect_col)[1]),,drop = FALSE]
    }      
    
    
    ## If the nummber of columns is not equal to the number of station locations + 2, try alternatives
    ## to search for column names
    col_search <- map_lgl(col_coord$text, ~str_detect(string = toupper(station_col_names), pattern = toupper(.x)))
    
    # The number of col_search non-matches should be 2
    if (sum(!col_search) != 2) {
      
      # Find alternative column identifiers
      col_coord_alt <- page_dat %>% 
        filter(str_detect(text, "C\\.V\\.|CV|Average|Mean")) %>%
        mutate(y = 0) %>%
        top_n(x = ., n = -1, wt = x) %>%
        head(1)
    
      ## Add this to col_coord, remove y
      col_coord <- bind_rows(col_coord, col_coord_alt) %>%
        arrange(x)
        
    }
    
    
    # Determine the max line based on summary values
    last_y <- page_dat %>%
      filter(str_detect(toupper(text), "LSD|CV|C\\.V\\.|AVERAGE|MEAN"),
             y > mean(y)) %>%
      pull(y)
    last_y <- ifelse(length(last_y) == 0, max(page_dat$y), min(last_y))
    
    page_dat_post_col <- page_dat %>%
      filter(y > max(col_coord$y)) %>%
      # Remove the "Average" rows
      filter(y < last_y)
    
    ## Since the data is more or less centered within a cell,
    ## we can use the average x coordinate of each column name (x + (x + width) / 2).
    ## Then for each column, we find data that overlaps with this average. We take the longest
    ## data string (or the column names if they are longer) and use these to demark the columns in
    ## the table
    
    col_data_list <- list()
    
    ## List over columns
    for (j in seq(nrow(col_coord))) {
      
      ## Determine the x value for the left hand cut
      ## x1 is left, x2 is right
      # If j == 1, the cut x is the minimum x
      if (j == 1) {
        x1 <- min(page_dat$x)
                
      } else {
        x1 <- col_coord$x[j]
        
      }
      
      ## To figure out the right hand side (x2), we need to know
      ## which is longer:
      ## 1) the length of the column name string
      ## 2) the length of the longest element string in that column
      
      # Extract the width of the column name
      width1 <- col_coord$width[j]
      
      ## Determine what text overlaps with that column
      ## The text to target should be between(x1, x, x+width)
      dat_use <- page_dat_post_col %>% 
        filter(map2_lgl(x, width, ~any(seq(.x, .x + .y) %in% seq(x1, x1 + width1))))
      
      # Which is longer, the column name or the longest data string?
      # This is x2
      x2 <- pmax(x1 + width1, max(dat_use$width + dat_use$x))
      
      # Which is less, x1 (above) or the leftmost x of text in this column?
      # This is the new x1 
      x1 <- pmin(x1, min(dat_use$x))
      
      ## Find data that falls within this range
      ## Add to the list
      col_data <- filter(page_dat_post_col, between(x, x1, x2)) %>%
        # split by y and concatenate
        split(.$y) %>%
        map_df(~tibble(y = unique(.x$y), text = paste0(.x$text, collapse = " "))) %>%
        # Add min and max x
        mutate(min_x = x1, max_x = x2)
      
      ## If j == 1, determine the average distance between rows
      ## Take half of this distance and round down
      if (j == 1) { 
        avg_row_dist <- floor(mean(diff(col_data$y)) / 2) 
      } 
      
      
      col_data_list[[j]] <- col_data
      
    }
    

    ## Merge by y
    ## Add some wiggle room using the average half distance between rows
    col_data_df <- col_data_list %>%
      # Remove x
      map(., ~select(., -min_x, -max_x)) %>%
      reduce(., coord_join, by = "y", tol = avg_row_dist / 2, remove.NA = T) %>%
      select(-y) %>%
      # Add column names
      rename_all(~make.names(names = col_coord$text, unique = TRUE)) %>%
      mutate_all(parse_guess) %>%
      filter_at(vars(1), ~!is.na(.))
    
    
    ## Check the columns again
    ## If the nummber of columns is not equal to the number of station locations + 2, try alternatives
    ## to search for column names
    col_search <- map_lgl(names(col_data_df), ~str_detect(string = toupper(station_col_names), pattern = toupper(.x)))
    
    # The number of col_search non-matches should be 2
    if (sum(!col_search) != 2) {
      
      ## Check if it's a duplication issue
      dup <- any(map_lgl(names(col_data_df)[!col_search], ~str_detect(string = ., station_col_names)))
      if (!dup) { 
      
        ## In this case, look for CV as the first column name, and try to extract numbers
        names(col_data_df)[1] <- "CV"
        
        # Extract numbers
        col_data_df <- col_data_df %>%
          mutate(number = str_extract(.[[1]], "^[0-9]*"),
                 selection = str_remove(.[[1]], "^[0-9]*") %>% str_trim() ) %>%
          select(-1) %>%
          select(number, selection, names(.)) %>%
          filter_at(vars(1), ~!is.na(.) & . != "")
      }
      
    }
    
    
    
    ## Extract location-specific metadata
    ## These are things like CV, LSD, mean, etc.
    
    # Get the min and max x coordinates from the trait data df
    col_data_coord <- col_data_list %>%
      map_df(~distinct(., min_x, max_x)) %>%
      bind_cols(col_coord, .)
    
    # Determine the min and max line based on summary values
    # Search for CV, C.V., LSD, Average, Mean
    metadata_search <- page_dat %>%
      filter(str_detect(toupper(text), "CV|C\\.V\\.|LSD|MEAN|AVERAGE"))
    
    if (max(page_dat_post_col$y) >= max(metadata_search$y)) {
      min_y_threshold <- subset(page_dat, str_detect(text, last(col_data_df[[2]])), y, drop = T) + 1
        
    } else {
      min_y_threshold <- max(page_dat_post_col$y) + 1
      
    }

    min_y <- page_dat %>%
      # This y must be >= the max y from the previous extraction
      # Unless this max is also the max of the page
      filter(y >= min_y_threshold) %>%
      pull(y) %>%
      min()
    
    # Max y is the max y of the entire page or the the max y before encountering astericks
    max_y <- page_dat %>% 
      filter(str_detect(toupper(text), "\\*"), y >= min_y) %>% 
      pull(y) %>% 
      max() %>%
      abs()
    # Take off 1 unit so it is excluded
    max_y <- pmin(max_y - 1, max(page_dat$y))

    page_dat_post_col <- page_dat %>%
      filter(between(y, min_y, max_y))
    
    ## Since the data is more or less centered within a cell,
    ## we can use the average x coordinate of each column name (x + (x + width) / 2).
    ## Then for each column, we find data that overlaps with this average. We take the longest
    ## data string (or the column names if they are longer) and use these to demark the columns in
    ## the table
    
    # Empty list to store extracted data
    col_metadata_list <- list()
    avg_row_dist <- NULL
    
    ## List over columns
    for (j in seq(nrow(col_data_coord))) {
      
      ## Get x1 and x2
      x1 <- col_data_coord$min_x[j]
      x2 <- col_data_coord$max_x[j]
      
      ## Find data that falls within this range
      ## Add to the list
      col_data <- filter(page_dat_post_col, between(x, x1, x2)) %>%
        # split by y and concatenate
        split(.$y) %>%
        map_df(~tibble(y = unique(.x$y), text = paste0(.x$text, collapse = " ")))
      
      # If avg_row_dist is not filled yet, fill it
      if (is.null(avg_row_dist) & nrow(col_data) > 0) {
        avg_row_dist <- floor(mean(diff(col_data$y)) / 2) 
      }
      
      col_metadata_list[[j]] <- col_data
      
    }
    
    # If avg_row_dist is NA, make it 0
    avg_row_dist <- ifelse(is.na(avg_row_dist), 0, avg_row_dist)
    
    
    ## Remove empty elements (only if first)
    element_to_keep <- which(map_lgl(col_metadata_list, ~nrow(.) > 0))
    
    ## Merge by y
    ## Add some wiggle room using the average half distance between rows
    col_metadata_df <- col_metadata_list[element_to_keep] %>%
      # Remove any lines in the first element that are just numbers
      modify_at(.x = ., .at = 1, ~filter(., str_detect(text, "05", negate = TRUE))) %>%
      reduce(., coord_join, by = "y", tol = avg_row_dist / 2, remove.NA = T) %>%
      select(-y) %>%
      # Add column names
      rename_all(~make.names(names = col_coord$text[element_to_keep], unique = TRUE)) %>%
      filter_at(vars(1), ~!is.na(.))

    ## Add to the list
    trait_data_list[[table_name]] <- list(data = col_data_df, metadata = col_metadata_df)
    
  }
  
  
  # Combine into a single df
  trait_data_df <- trait_data_list %>%
    map("data") %>%
    imap_dfr(~{
      # Determine character columns
      char_col <- summarize_all(.x, is.character) %>% unlist() %>% which()
      
      # Convert data from character to numeric
      mutate_at(.x, vars(setdiff(char_col, 2)), parse_number) %>%
        mutate(table_name = .y, trait = str_extract(.y, trait_pattern)) %>% 
        ## Rename "entry" as "number"
        rename_all(~str_replace_all(string = ., pattern = "Entry", replacement = "Number"))
      
    }) %>%
    rename_at(vars(1:2), tolower) %>%
    gather(location, value, -1, -2, -table_name, -trait) %>%
    filter(!is.na(value))
  
  
  ## Create renaming df to rename locations
  location_rename <- trait_data_df %>%
    distinct(location) %>%
    # Detect duplicates
    mutate(dup = str_detect(location, "\\.[0-9]{1,}"),
           location1 = map2_chr(location, dup, ~ifelse(.y, str_remove_all(.x, "\\.[0-9]{1,}"), .x)),
           location_new = map(location1, ~str_subset(string = unique(station_info_df$location), pattern = .))) %>%
    # Filter out length 0 vectors
    filter(map_lgl(location_new, ~length(.) > 0)) %>%
    unnest(location_new) %>%
    # Add dup number back in
    mutate(location_new = map2_chr(location, location_new, 
                                   ~ifelse(str_detect(.x, "\\.[0-9]{1,}"), paste0(.y, str_extract(.x, "\\.[0-9]{1,}")), .y))) %>%
    select(location, location_new)
  
  
  ## Combine dfs and rename locations
  trait_data_df1 <- trait_data_df %>%
    right_join(., location_rename) %>%
    select(-location) %>%
    rename(location = location_new)
  
  ## Metadata
  # Combine into a single df
  trait_metadata_df <- trait_data_list %>%
    map("metadata") %>%
    imap_dfr(~{
      
      ## Reshape
      gather(.x, location, value, -1) %>% 
        spread(1, value) %>% 
        mutate(table_name = .y,
               trait = str_extract(.y, trait_pattern)) %>%
        rename_all(make.names) %>%
        select(-starts_with("X")) %>%
        rename_all(~str_remove_all(string = ., pattern = "[:punct:]|[0-9]*") %>% str_trim())

    })
      

  
  #############################
  ## Assemble information in a list
  #############################
  
  extracted_data_list[[basename(pdf_file)]] <- list(
    filename = basename(pdf_file), year = report_year, 
    station_info = station_info_df,
    pedigree = pedigree_info_df, 
    trait_data = trait_data_df1,
    trait_metadata = trait_metadata_df
  )
  
  # Print
  cat("\nData extracted for WRBN for report year: ", report_year)
  
} # Close the loop


## Save this data
wrn_extracted_data_list <- extracted_data_list






##### Add additional reports ####

## Date formats for parsing
date_formats <- c("mdy", "%B %d", "%B %d, %Y", "j", "%b %d")


wrn_excel_report_files <- list.files(agro_dir, pattern = "WR", full.names = TRUE) %>%
  str_subset(., "\\~\\$", negate = TRUE)

# Iterate over the files
for (i in seq_along(wrn_excel_report_files)) {
  
  # Filename
  filename <- wrn_excel_report_files[i]
  report_year <- str_extract(string = basename(filename), pattern = "[0-9]{4}")
  
  # List sheets
  file_sheets <- excel_sheets(path = filename)
  
  ### Entry list data ###
  pedigree_sheet <- str_subset(string = file_sheets, pattern = "Entry")
  pedigree_raw_dat <- read_excel(path = filename, sheet = pedigree_sheet)
  
  # Use the pedigree renaming column to determine the greatest number
  # of the row that matches the pattern
  ped_colname_last_row <- pedigree_raw_dat %>%
    map(~str_which(string = ., pattern = paste0(names(ped_col_rename), collapse = "|"))) %>% 
    map_dbl(~ifelse(length(.) == 0, NA, max(.)))
  
  # Subset rows and column; rename
  pedigree_raw_dat1 <- pedigree_raw_dat[-seq(1, max(ped_colname_last_row, na.rm = TRUE)), !is.na(ped_colname_last_row), drop = FALSE] %>%
    rename_all(~map_chr(pedigree_raw_dat, max(ped_colname_last_row, na.rm = TRUE))[!is.na(ped_colname_last_row)]) %>%
    # Remove rows with all NA
    filter_all(any_vars(!is.na(.))) %>%
    # Remove punctuation
    rename_all(~str_remove_all(., "[:punct:]")) %>%
    # Rename
    rename_all(~str_replace_all(string = ., ped_col_rename)) %>% 
    ## Parse
    mutate_at(vars(which(names(.) != "Number")), parse_character)
  
  ## Edit the pedigree data
  ## Edit pedigree data
  ## Combine into a single DF
  pedigree_info_df <- pedigree_raw_dat1 %>%
    rename_all(tolower) %>%
    mutate(entry = str_remove(entry, "\\*") %>% str_trim(),
           alias = str_extract(entry, "\\(.*\\)"),
           alias = str_remove_all(alias, "\\(|\\)"),
           entry = str_trim(str_remove(entry, "\\(.*\\)"))) %>%
    # Capitalize all
    mutate_if(is.character, toupper) %>%
    # Remove completely empty rows
    filter_all(any_vars(!. == "")) %>%
    mutate(nursery = toupper(names(pedigree_raw_dat)[1])) %>%
    filter(!is.na(entry))
  
  
  ### Trait data ###
  # Sheets with data have a location, state pattern to the name of the sheet
  trait_data_sheets <- str_subset(string = file_sheets, ", ")
  trait_data_list <- meta_data_list <- list()
  
  # Iterate over these sheets
  for (p in seq_along(trait_data_sheets)) {
    
    # Read in the data
    trait_data_raw_p <- read_excel(path = filename, sheet = trait_data_sheets[p])
    
    # If the data.frame is empty, move on
    if (is_empty(trait_data_raw_p)) next
    
    ## Metadata
    metadata_patterns <- c(location = "Location:", planting_date = "Seed Date:|Planting Date:",
                           harvest_date = "Harvest Date:")
    
    metadata_list <- map(metadata_patterns, ~{
      pattern <- .x
      map(trait_data_raw_p, ~str_subset(., pattern)) %>% 
        unlist() %>% 
        str_remove(string = ., pattern = pattern) %>%
        str_trim() })
    
    # Convert to df
    metadata_df <- as_tibble(metadata_list) %>%
      # Add the sheet location name if it is shorter than what is recorded or if 
      # what is recorded is blank
      mutate(raw = ifelse(location == "" | nchar(location) > max(nchar(trait_data_sheets)), 
                               trait_data_sheets[p], location)) %>%
      mutate(state = str_extract(raw, ", .*") %>% str_extract(., "[A-Z]{2}"),
             location = str_split(raw, ", ")[[1]][1]) %>%
      slice(1)
    
    meta_data_list[[p]] <- metadata_df
    
    
    # Search for the column/row containing "cultivar/designation"
    cult_search_text <- map(trait_data_raw_p, ~str_extract(., "DESIGNATION|Cultivar")) %>%
      map(~str_subset(., "DESIGNATION|Cultivar")) %>%
      unlist()
    
    cult_search <- map(trait_data_raw_p, ~str_which(., cult_search_text))
    
    col <- min(which(map_lgl(cult_search, ~length(.) > 0)))
    row <- min(cult_search[[col]])
    row <- ifelse(cult_search_text == "Cultivar", row - 1, row)
    
    trait_data_raw_p1 <- trait_data_raw_p[-seq(1, row), -seq(1, col - 1), drop = FALSE]
    # Remove heavily NA rows
    trait_data_raw_p2 <- trait_data_raw_p1[!apply(X = trait_data_raw_p1, MARGIN = 1, FUN = function(r) mean(is.na(r))) >= 0.90, , drop = FALSE]
      
    # Search for trait units
    unit_pattern <- c(GrainYield = "Mg ha-1|Mg/Ha", GrainYield = "bu a-1|Bu/A", TestWeight = "kg m-3|kg/m3", 
                      TestWeight = "lb bu", HeadingDate = "from 1/1|Julian", PlantHeight = "cm", PlumpGrain = ">2.4mm|>2.2mm", 
                      ThinGrain = "<2.2mm", GrainProtein = "%")
    
    # Select columns that match these patterns
    trait_cols <- unit_pattern %>%
      map(~{y <- .; which(map_lgl(trait_data_raw_p2[-1], ~any(str_detect(., y), na.rm = TRUE))) }) %>%
      # If that column is mostly NA, skip
      map(function(y) map(y, ~ifelse( mean(is.na(trait_data_raw_p2[-1][[.]])) > 0.8, NA, .) )) %>%
      map_dbl(~ifelse(any(is.na(.x)), NA, first(.))) %>%
      subset(., !is.na(.))
      
    # Make sure one column per trait is selected
    trait_cols <- map(unique(names(trait_cols)), ~trait_cols[which(names(trait_cols) %in% .x)]) %>% 
      map_dbl(max) %>%
      set_names(., unique(names(trait_cols)))
    
    trait_data_raw_p3 <- trait_data_raw_p2 %>%
      select(c(line_name = 1, trait_cols + 1)) %>%
      filter(!is.na(line_name)) %>%
      # Parse
      mutate_all(str_trim) %>%
      mutate_if(is.character, parse_guess) %>%
      select_if(~!all(is.na(.)))
    
    # Extract only entry data
    trait_data_df <- trait_data_raw_p3 %>%
      filter(toupper(line_name) %in% toupper(pedigree_info_df$entry)) %>%
      # Add metadata
      mutate(location = metadata_df[["location"]],
             management = "irrigated") %>%
      gather(trait, value, -line_name, -location, -management) %>%
      mutate_if(is.character, parse_guess)
    
    # Add to data list
    trait_data_list[[p]] <- trait_data_df
    
  }
    
  
  ## Bind rows for metadata, traits, etc.
  metadata_combine <- bind_rows(meta_data_list) %>%
    mutate_all(str_trim) %>%
    # Reformat dates
    mutate_at(vars(contains("date")), ~parse_date_time(., orders = date_formats) %>% `year<-`(., as.numeric(report_year))) %>%
    mutate_if(is.POSIXct, as.character)

  trait_data_combine <- bind_rows(trait_data_list)
    
  ## Add data to the extracted data list
  extracted_data_list[[basename(filename)]] <- list(
    filename = basename(filename),
    year = report_year,
    station_info = metadata_combine,
    pedigree = pedigree_info_df,
    trait_data = trait_data_combine
  )
  
}


wrn_extracted_data_list <- extracted_data_list


save("wrn_extracted_data_list", file = file.path(extr_dir, "wrn_extracted_metadata.RData"))


#############################
## Clean up the results
#############################

## Read in PDF text to aid in cleaning
all_pdf_text <- map(pdf_files, ~pdf_text(pdf = .) %>% str_to_upper())

# Load the data
load(file.path(extr_dir, "wrn_extracted_metadata.RData"))



#### Pedigree

## 2013 has a very odd pedigree table
## Manually curate it

## Import this manual table
wrn2013pedigree <- read_excel(path = file.path(extr_dir, "wrn2013_pedigree_manual.xlsx"), sheet = "pedigree")


## Pull the pedigree data and assemble into a single df
pedigree_df <- map(wrn_extracted_data_list, "pedigree") %>% 
  map(~mutate_if(., is.character, parse_guess)) %>%
  imap_dfr(~mutate(.x, file = .y)) %>%
  rename(entry = number, name = entry) %>%
  mutate(report = str_remove(file, ".pdf") %>% str_replace(., " ", "_")) %>%
  select(-file) %>%
  as_tibble() %>%
  # Replace 2013
  filter(str_detect(report, "2013", negate = T)) %>%
  bind_rows(., wrn2013pedigree) %>%
  # parse all characters and remove strings
  mutate_if(is.character, parse_guess) %>%
  mutate_all(str_trim)

## Output this raw metadata information
write_csv(x = pedigree_df, path = file.path(extr_dir, "wrn_extracted_entry_metadata.csv"))



# Remove spaces and capitalize entry names
pedigree_df_distinct_entry <- pedigree_df %>% 
  mutate_at(vars(name, pedigree, alias), ~str_remove_all(., " ") %>% toupper()) %>%
  rename(line_name = name) %>%
  # distinct line name
  distinct(line_name)

## Find distinct pedigrees, aliases, program, row type, grade
cols_distinct <- c("pedigree", "alias", "type", "grade")

pedigree_distinct_cols <- pedigree_df %>%
  mutate(type = parse_number(type)) %>%
  select(line_name = name, cols_distinct)

pedigree_distinct_cols1 <- cols_distinct %>%
  map(~select(pedigree_distinct_cols, line_name, .)) %>%
  map(~filter_all(., all_vars(!is.na(.)))) %>%
  map(distinct) %>%
  c(list(pedigree_df_distinct_entry), .) %>%
  reduce(left_join, by = "line_name")


## Check for duplicate line names with type and grade; fix problems
pedigree_distinct_cols1 %>%
  group_by(line_name) %>% 
  filter(n_distinct(type) > 1) %>% 
  arrange(line_name)

pedigree_distinct_cols1 %>%
  group_by(line_name) %>% 
  filter(n_distinct(grade) > 1) %>% 
  arrange(line_name)

## Clean up some columns
pedigree_distinct_cols2 <- pedigree_distinct_cols1 %>%
  mutate(type = ifelse(line_name == "STEPTOE", 2, type)) %>%
  # Extract types from grade
  mutate(grade = str_extract(toupper(grade), "MALT|FEED|FOOD")) %>%
  distinct()
  
# Grade priority
grade_priority <- c("MALT", "FEED", "FOOD")

# For each line name, select non-NA pedigree, alias, feed, and grade
# Prioritize grade and type
pedigree_distinct_cols3 <- pedigree_distinct_cols2 %>%
  split(.$line_name) %>%
  # For each line, remove some rows:
  # If a column is incompletely NA, filter out NA
  # Priortize malt, feed, food
  .[order(map_dbl(., nrow), decreasing = TRUE)]  %>%
  map_df(~{
    # Are there any columns with missing data > 0 & < 1?
    filter_cols <- map_dbl(.x, ~mean(is.na(.))) %>% 
      subset(., . < 1 & . > 0) %>%
      names()
    
    if (!is_empty(filter_cols)) {
      .x2 <- filter_at(.x, .vars = vars(filter_cols), .vars_predicate = all_vars(!is.na(.)))
    } else {
      .x2 <- .x
    }
    
    # Find unique grades, if all are not NA
    if (!all(is.na(.x2$grade))) {
      mutate(.x2, grade = factor(grade, levels = grade_priority), grade_no = as.numeric(grade)) %>%
        group_by_at(vars(-grade, -grade_no)) %>% 
        top_n(x = ., n = 1, wt = -grade_no) %>% 
        ungroup() %>% 
        select(-grade_no)
      
    } else {
      .x2
    }
    
  })
    



## Check for duplicate name with unmatching pedigrees
pedigree_distinct_cols3 %>% 
  group_by(line_name) %>% 
  filter(n_distinct(pedigree) > 1) %>% 
  arrange(line_name) %>%
  as.data.frame()

## Leave these duplicates in; the pedtools program will find them

## Check for other duplicates
pedigree_distinct_cols3 %>%
  group_by(line_name) %>%
  filter_at(vars(alias, type, grade), any_vars(n_distinct(.) > 1))

# none - good



# Re-write the CSV
write_csv(x = pedigree_distinct_cols3, path = file.path(extr_dir, "wrn_extracted_entry_metadata_edited.csv"))



#### Station information

## Load location information
location_info_df <- wrn_extracted_data_list %>%
  map(`[`, c("year", "station_info")) %>%
  transpose() %>%
  # Add year to the df
  pmap(~mutate(.y, year = .x)) %>%
  imap_dfr(~mutate(.x, file = .y)) %>%
  select(location, state, year, contains("date"), nursery, raw, file)

# Save this
write_csv(x = location_info_df, path = file.path(extr_dir, "wrn_extracted_location_metadata.csv"))


## Select relevant columns, make some edits, assign management
location_info_df1 <- location_info_df %>% 
  mutate(report = str_remove(file, ".pdf") %>% str_replace(., " ", "_"),
         location = str_replace_all(location, "_", " "),
         location = str_remove_all(location, "And|and"),
         location = str_to_title(str_trim(location))) %>%
  select(location, state, year, contains("date"), report)

## Output distinct location/state for corrections
location_info_df %>% 
  distinct(location, state) %>%
  arrange(location) %>%
  write_csv(x = ., path = file.path(extr_dir, "wrn_extracted_location_rename.csv"))


## Custom location renaming df
location_rename <- read_excel(path = file.path(extr_dir, "wrn_extracted_location_rename.xlsx"))

# Rename
location_info_df2 <- location_info_df1 %>% 
  inner_join(., location_rename) %>%
  select(location = location_rename, state = state_rename, year, contains("date"), report) %>%
  filter_at(vars(location, state), all_vars(!is.na(.))) %>%
  distinct()

## Save csv
write_csv(x = location_info_df2, path = file.path(extr_dir, "wrn_extracted_station_info.csv"))













