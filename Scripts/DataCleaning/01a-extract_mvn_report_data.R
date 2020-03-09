## Code to extract metadata from the MVN reports
## 
## Jeff Neyhart
## 

# Packages to use
packages <- c("tidyverse", "lubridate", "readxl", "pdftools")
invisible(lapply(packages, library, character.only = TRUE))

## Source a script for relevant functions
source(file.path(getwd(), "source_functions.R"))


# Trait metadata
trait_metadata <- read_csv("C:/GoogleDrive/BarleyLab/Breeding/PhenotypicData/Metadata/trait_metadata.csv")


# Directories containing data
shared_drive_dir <- c(paste0(LETTERS, ":/AGRO/BARLEY_LAB"), paste0(LETTERS, ":/BARLEY_LAB")) %>%
  subset(., dir.exists(.)) %>%
  head(1)

nursery_dir <- file.path(shared_drive_dir, "/Breeding/BARLEY/HistoricalNurseryData/")
raw_dir <- file.path(nursery_dir, "RawData")
agro_dir <- file.path(raw_dir, "Agronomic")
maltq_dir <- file.path(raw_dir, "Quality")
tidy_dir <- file.path(nursery_dir, "TidyData")
extr_dir <- file.path(nursery_dir, "ExtractedData")


# Change working directory to BARLEY LAB directory
mvn_dir <- file.path(nursery_dir, "PDF_Reports/MississippiValleyNursery")
report_dir <- mvn_dir



#################################
# Read and scrape PDF reports
#################################

## List all report files
pdf_files <- sort(list.files(report_dir, pattern = "MVN [0-9]{4}.pdf", full.names = TRUE))


# Create an empty list to store output
extracted_data_list <- list()

## Iterate over files
for (p in seq_along(pdf_files)) {
  
  pdf_file <- pdf_files[p]
  
  ## Extract text
  pdf_txt <- pdf_text(pdf = pdf_file) %>%
    # convert everything to upper-case
    str_to_upper()
  
  # Get coordinates of text
  pdf_dat <- pdf_data(pdf = pdf_file) %>%
    map(~mutate(., text = toupper(text)))
  
  # Which is the contents page?
  toc_page <- which(str_detect(pdf_txt, "CONTENTS"))

  
  
  ## Get the pages of all tables
  table_names <- str_extract(pdf_txt, "TABLE [0-9].*\r") %>%
    str_remove_all(string = ., pattern = "\n|\r")
  
  table_names[toc_page] <- NA
  
  
  #############################
  ## Year
  #############################
  
  report_year <- str_extract(pdf_txt[[1]], "[0-9]{4}")
  
  
  
  #############################
  ## Station Names
  #############################
  
  # Station names are organized by a table 
  station_page <- str_which(pdf_txt, "LIST OF STATIONS")
  
  # Get the data from this page
  station_data <- pdf_dat[[station_page]]

  ### Extract the table
  
  ## Search for column headers
  # The first column is station; the second column is cooperator.
  # This might vary depending on the report 
  pattern <- "STATIONS|COOPERATORS"
  col_coord <- station_data %>% 
    filter(str_detect(text, pattern)) %>% 
    filter(y == max(y))  %>%
    # Add rightmost x
    mutate(xend = c(tail(x, -1) - 1, max(station_data$x)))
  
  # Filter the rows beyond the columns
  page_dat_post_col <- station_data %>% 
    # Move past the columns
    filter(y > unique(col_coord$y))

  
  ## Elements in a column are left justified; so elements 
  ## within a column can be extracted by common x. Merge by common y
  
  ## Create a list to store the column cut coordinates
  col_data_list <- list()
  ## List over columns
  for (j in seq(nrow(col_coord))) {
    
    # Get x
    x1 <- col_coord$x[j]
    xend <- col_coord$xend[j]
    # Get the leftmost elements of the column
    left_x1 <- page_dat_post_col %>%
      filter(between(x, x1, xend)) %>%
      # Merge by common y
      split(.$y) %>%
      map_df(~tibble(y = unique(.x$y), text = paste0(.x$text, collapse = " ")))
    
    ## Add to list
    col_data_list[[j]] <- left_x1
  }
    
    
  ## Convert to df
  station_df <- reduce(col_data_list, coord_join, by = "y", tol = 0) %>% 
    select(-y) %>% 
    rename_all(~col_coord$text) %>%
    filter_all(all_vars(!is.na(.))) %>%
    rename_all(tolower)
  
  # Separate into location/state
  station_df1 <- station_df %>%
    filter(stations != "") %>%
    separate(stations, c("location", "state_province"), sep = ", ") %>%
    ## Arrange by state_province, then location
    arrange(state_province, location) %>%
    ## Assign this as the order
    mutate(location = factor(location, levels = .$location)) %>%
    # Remove NA
    filter_all(all_vars(!is.na(.)))
  
  
  
  
  
  ## We need to assign numbers to the stations. To do this, we need to determine whether 
  ## data from those stations was submitted. We can parse Table 5 (information on
  ## plot type) for this.
  
  # Get the page
  plot_type_page <- str_which(string = table_names, pattern = "PLOT TYPE")
  
  # Get column coordinates
  page_dat <- pdf_dat[[plot_type_page]]
  col_coord <- filter(page_dat, str_detect(text, "LOCATION|\\#|REPS|SOWN|KG|METERS|HARVESTED|ROWS"))
  
  # If more than one harvested, choose second
  col_coord1 <- anti_join(
    x = col_coord, 
    y = col_coord %>% 
      filter(text == "HARVESTED") %>% 
      filter(x != max(x)), by = c("width", "height", "x", "y", "space", "text")) %>%
    arrange(x)
  
  # If # and REPS are both present, choose reps
  x_keep <- col_coord1 %>%
    filter(str_detect(text, "\\#|REPS"))
  x_remove <- if (length(x_keep) > 1) subset(x_keep, str_detect(text, "REPS", negate = T), x, drop = T) else numeric(0)

  if (length(x_remove > 0)) {
    col_coord1 <- subset(col_coord1, x != x_remove)
  }
  
  
  # Filter the rows beyond the columns
  page_dat_post_col <- page_dat %>% 
    # Move past the columns
    filter(y > max(unique(col_coord1$y)))
  
  
  ## Iterate over columns
  ## Create a list to store the column cut coordinates
  col_data_list <- list()
  
  ## List over columns
  for (j in seq(nrow(col_coord1))) {
    
    # Determine what text overlaps with that column
    # x1 is first column
    # x2 is the next column
    
    ## If j == 1, x1 is the min x
    if (j == 1) {
      x1 <- min(page_dat$x)
      
    } else {
      x1 <- col_coord1$x[j]
      
    }
    
    # Width of the column name
    width1 <- col_coord1$width[j]
    
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
      map_df(~tibble(y = unique(.x$y), text = paste0(.x$text, collapse = " ")))
    
    col_data_list[[j]] <- col_data
    
  }
    
  # Filter out nrow = 0
  to_keep <- which(map_lgl(col_data_list, ~nrow(.) > 0))
  
  # Keep this in the data list and the column name list
  col_data_list <- col_data_list[to_keep]
  col_data_names <- col_coord1$text[to_keep]
    
  ## Re-order the list so the longest column is first
  col_orders <- order(map_dbl(col_data_list, nrow), decreasing = TRUE)
  col_data_list <- col_data_list[col_orders]
  # Reorder the names, too
  col_data_names <- col_data_names[col_orders]
  
  ## Merge by y
  ## Use coordinate join with a tolerance for the y value
  plot_type_df <- reduce(col_data_list, coord_join, by = "y", tol = 0) %>%
    select(-y) %>%
    # Add column names
    rename_all(~col_data_names) %>%
    # Remove punctuation
    rename_all(~str_remove_all(., "[:punct:]")) %>%
    # Lowercase and select important
    rename_all(tolower) %>%
    select(location, sown, harvested) %>%
    ## Parse dates
    mutate_at(vars(sown, harvested), ~parse_date_time(x = ., orders = c("%m %d")) %>% 
                `year<-`(., as.numeric(report_year)) %>% as.character()) %>%
    # Rename
    rename(planting_date = sown, harvest_date = harvested)
  
  
  ## Order based on levels in the station_df, assign numbers
  station_plot_df <- plot_type_df %>%
    filter(!is.na(planting_date)) %>%
    mutate(location = factor(location, levels = levels(station_df1$location))) %>%
    arrange(location) %>%
    filter(!is.na(location)) %>%
    mutate(number = seq(nrow(.))) %>%
    droplevels()
    
  
  
  
  #############################
  ## Determine the location names
  #############################
  
  ## Yield tends to be a trait with universal coverage, so use those pages
  
  # Find the yield pages
  yield_pages <- str_which(table_names, "AVE\\. STA\\. YIELD \\(KG/HA\\)")
  
  # Get the data from those pages
  yield_pages_table_data <- pdf_dat[yield_pages]
  location_name_list <- list()
  
  ## Iterate over the pages
  for (i in seq_along(yield_pages)) {
    
    ## Extract the data
    data_i <- yield_pages_table_data[[i]]
    # Get the y values of "AVERAGE", since this line will contain the locations
    data_i_sub <- subset(data_i, text == "AVERAGE")
    data_i_sub1 <- subset(data_i, str_detect(text, "STA\\."))
    # Get the coordinates
    dim_use <- unlist(subset(data_i_sub, y == min(y), c(x, y)))
    dim_use1 <- unlist(subset(data_i_sub1, x >= dim_use[1], c(x, y)))
    
    
    ## Data from the location line
    loc_data_i <- subset(data_i, x > dim_use[1] & y == dim_use[2])
    ## Data from the state line
    state_data_i <- subset(data_i, x > dim_use1[1] & y == dim_use1[2]) %>%
      mutate(i = seq(nrow(.)))
    
    ### Use the x values from the state df to determine if location names need to be combined ###
    # First, determine the average width of locations
    avg_log_width <- mean(loc_data_i$width)
    # Are there outliers? (if widths are less than 10 and are adjacent, combine)
    loc_data_i1 <- loc_data_i %>%
      mutate(outlier = width < avg_log_width - 10, i = seq(nrow(.)))
    
    # Combine if the number of outliers is > 1
    
    if (sum(loc_data_i1$outlier) > 1) {
      
      ## Determine groups to combine
      loc_outliers <- filter(loc_data_i1, outlier) %>%
        mutate(diff = c(1, diff(i)), group = i)
      
      for (q in seq(2, length(loc_outliers$i))) {
        group_q <- if (loc_outliers$diff[q] == 1) loc_outliers$group[q - 1] else loc_outliers$i[q]
        loc_outliers$group[q] <- group_q
      }
      
      ## Combine back to the loc_data
      loc_data_i2 <- left_join(loc_data_i1, loc_outliers, by = c("width", "height", "x", "y", "space", "text", "outlier", "i")) %>% 
        mutate(group = ifelse(is.na(group), i, group))
      
      ## Split on group and add back to the df
      location_name_list[[i]] <- sapply(split(loc_data_i2$text, loc_data_i2$group), paste0, collapse = "_")
      
      
    } else {
      # Subset data for those coordinates and get the location names
      location_name_list[[i]] <- subset(data_i, x > dim_use[1] & y == dim_use[2], text, drop = TRUE)
      
    }
    
    
  }
  
  ## Unlist the location names and assign numbers
  location_df <- tibble(location = unlist(location_name_list), environment_number = seq_along(location))
  
  
  
  
  
  #############################
  ## Scrape trait data
  #############################
  
  # First get the pages containing relevant tables
  trait_tables <- str_subset(table_names, "AVE\\. STA\\. [A-Z ]* \\(")
  trait_pages <- str_which(table_names, "AVE\\. STA\\. [A-Z ]* \\(")
  
  # A list to store df
  trait_data_list <- list()
  
  ## A pattern vector to match
  col_pattern <- c("NO\\.", "SELECTION", "RANK", "AVERAGE", location_df$location) %>%
    str_replace_all(string = ., "_", " ") %>%
    str_split(string = ., pattern = " ") %>% # Split on spaces
    map_chr(., 1) %>%
    paste0(collapse = "|")
  
  
  ## Iterate over the content df
  for (i in seq_along(trait_pages)) {
    
    # Get the page
    page <- trait_pages[i]
    
    # Get the text and coordinates at this page
    page_text <- pdf_txt[[page]]
    page_dat <- pdf_dat[[page]]
    
    
    ## If the page contains "NO DATA", pass
    if (str_detect(page_text, "NO DATA")) next
    
    
    # What is the trait name?
    trait <- str_extract(page_text, "AVE\\. STA\\. [A-Z ]*") %>% 
      str_remove(., "AVE\\. STA\\.") %>% 
      str_trim()
    
    ## Create a table name
    table_name <- paste0(trait, "_page", page)
    
    ## Find the column names and coordinates
    col_coord <- page_dat %>%
      # Remove special characters
      filter(str_detect(text, col_pattern)) %>%
      filter(y < 200) %>%
      arrange(x)
      
    # Determine the max line based on summary values
    last_y <- page_dat %>%
      filter(y > max(col_coord$y),
             str_detect(text, "AVERAGE"),
             y > mean(y)) %>%
      pull(y)
    last_y <- ifelse(length(last_y) == 0, max(page_dat$y), min(last_y))
    
    ## Subset the data to parse
    page_dat_post_col <- page_dat %>%
      filter(y > max(col_coord$y)) %>%
      # Remove the "Average" rows
      filter(y < last_y)
    
    
    col_data_list <- list()
    
    ## List over columns
    for (j in seq(nrow(col_coord))) {
      
      ## Determine the x value for the left hand cut
      ## x1 is left, x2 is right
      # If j == 1, the cut x is the minimum x
      if (j == 1) {
        x1 <- min(col_coord$x)
        
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
    
    ## Add the data to the list
    trait_data_list[[table_name]] <- col_data_df
    
  } # End loop
    
  
  # Combine into a single df
  trait_data_df <- trait_data_list %>%
    imap_dfr(~{
      # Determine character columns
      char_col <- summarize_all(.x, is.character) %>% unlist() %>% which()
      
      # Convert data from character to numeric
      mutate_at(.x, vars(setdiff(char_col, 2)), parse_number) %>%
        mutate(table_name = .y, trait = str_remove(.y, "_[a-z0-9]*$")) %>% 
        ## Rename "entry" as "number"
        rename_all(~str_replace_all(string = ., pattern = "Entry", replacement = "Number"))
      
    }) %>%
    rename(number = `NO.`, line_name = SELECTION) %>%
    select(-AVERAGE, -RANK) %>%
    gather(location, value, -1, -2, -table_name, -trait) %>%
    filter(!is.na(value))
  
  

  #############################
  ## Variety pedigree information
  #############################
  
  # Get the page
  pedigree_page <- str_which(string = table_names, pattern = "PARENTAGE|ORIGIN")
  # Text and data
  page_text <- pdf_txt[[pedigree_page]]
  page_dat <- pdf_dat[[pedigree_page]]
  
  
  # Vector of column names and renamings
  ped_col_rename <- c("^ENTRY|NTRY|^ENT" = "Number", "CONTRIBUTOR|^CI" = "Program", "NAME" = "Entry",
                      "PEDIGREE|PARENTAGE" = "Pedigree")
    
  ## Find the column names
  col_coord <- map_df(names(ped_col_rename), ~filter(page_dat, str_detect(text, .x)) %>%
                        group_by(y) %>% filter(width == max(width)) %>% ungroup()) %>%
    filter(str_detect(text, paste0(names(ped_col_rename), collapse = "|"))) %>%
    # subset for the most common y
    filter(y %in% names(which.max(table(y)))) %>%
    arrange(x)
  
  ## Find the max y value (where additional pedigree information is stored)
  ymax <- map(c("ADDITIONAL", "\\*E", "NEW"), ~min(subset(page_dat, str_detect(text, .x), y, drop = T))) %>%
    unlist() %>%
    subset(., is.finite(.)) %>%
    first()
    
    
  # Filter the rows beyond the columns
  page_dat_post_col <- page_dat %>% 
    # Move past the columns
    filter(y > unique(col_coord$y)) %>%
    # Stop at new entries
    filter(y < ymax)
  
  ## Data is inconsistently justified
  ## Create a list to store the column cut coordinates
  col_data_list <- list()
  
  ## List over columns
  for (j in seq(nrow(col_coord))) {
    
    # Determine what text overlaps with that column
    # x1 is first column
    # x2 is the next column
    
    ## If j == 1, x1 is the min x
    x1 <- if (j == 1) min(page_dat_post_col$x, min(col_coord$x)) else col_coord$x[j]
    
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

    # Merge within a column (sometimes rows can be off)
    col_data_list[[j]] <- col_data %>%
      mutate(ydiff = c(Inf, diff(y))) %>% 
      mutate(y = imap_dbl(ydiff, ~ifelse(.x < 5, col_data$y[.y - 1], col_data$y[.y]))) %>%
      split(.$y) %>%
      map_df(~tibble(y = unique(.x$y), text = paste0(.x$text, collapse = " ")))
    
  }
  
  # Filter out nrow = 0
  to_keep <- which(map_lgl(col_data_list, ~nrow(.) > 0))
  
  # Keep this in the data list and the column name list
  col_data_list <- col_data_list[to_keep]
  col_data_names <- col_coord$text[to_keep]
  
  ## Re-order the list so the longest column is first
  col_orders <- order(map_dbl(col_data_list, nrow), decreasing = TRUE)
  col_data_list <- col_data_list[col_orders]
  # Reorder the names, too
  col_data_names <- col_data_names[col_orders]

  
  ## Merge by y
  ## Use coordinate join with a tolerance for the y value
  pedigree_info_df <- reduce(col_data_list, coord_join, by = "y", tol = 5) %>%
    select(-y) %>%
    # Add column names
    rename_all(~col_data_names) %>%
    # Remove punctuation
    rename_all(~str_remove_all(., "[:punct:]")) %>%
    # Rename
    rename_all(~str_replace_all(string = ., ped_col_rename)) %>% 
    ## Parse
    mutate_at(vars(which(names(.) != "Number")), parse_character)

  
  
  
  ## Look for additional pedigree information beyond ymax
  
  # First determine if there is additional text beyond ymax
  ymax <- max(subset(page_dat, str_detect(text, "ADDITIONAL|\\*E|NEW"), y, drop = TRUE))
  addition_ped_dat <- subset(page_dat, y > ymax) %>%
    arrange(y, x)
  
  if (nrow(addition_ped_dat)) {
    
    ## Look for equals signs - these are pedigree strings
    eqls_i <- str_which(addition_ped_dat$text, "=")
    # For each equals, find the text that preceds it (this is the additional entry name)
    additional_ped_name <- addition_ped_dat[eqls_i - 1,]
    
    ## Find clusters of x
    x_list <- table(additional_ped_name$x) %>% 
      subset(. > 1) %>% 
      names() %>%
      as.numeric()
    # Add a final x
    x_list <- c(x_list, max(addition_ped_dat$x + 2))
    
    additional_ped_dat_list <- list()
    
    ## Iterate over x
    for (r in seq(1, length(x_list) - 1)) {
      xleft <- x_list[r]
      xright <- x_list[r + 1] - 1
      
      # Subset additional ped data by x
      additional_ped_dat_list[[r]] <- addition_ped_dat %>% 
        filter(x >= xleft, x < xright) %>%
        # split by y and concatenate
        split(.$y) %>%
        map_df(~tibble(y = unique(.x$y), text = paste0(.x$text, collapse = " "))) 
      
      
    }
    
    ## Merge additional pedigree data, then separate by the equals sign
    additional_ped_df <- additional_ped_dat_list %>%
      bind_rows() %>%
      # Remove any line with a colon
      filter(str_detect(text, ":", negate = TRUE)) %>%
      separate(text, c("Entry", "Pedigree"), sep = "=") %>%
      select(-y) %>%
      mutate_if(is.character, str_trim)
    
    
  } else {
    additional_ped_df <- NULL
    
  }
  
  
  
  
  
  ## Merge pedigree information
  ## Parse entry number
  pedigree_table_df1 <- bind_rows(pedigree_info_df, additional_ped_df) %>%
    rename_all(tolower) %>%
    ## Extract aliases
    mutate(alias = str_extract(entry, "\\(.*\\)"),
           alias = str_remove_all(alias, "\\(|\\)"),
           original = entry, # Keep original for
           entry = str_trim(str_remove(entry, "\\(.*\\)"))) %>%
    # Capitalize all
    mutate_if(is.character, toupper)
    
  
  
  #############################
  ## Assemble information in a list
  #############################
  
  extracted_data_list[[basename(pdf_file)]] <- list(
    filename = basename(pdf_file), 
    year = report_year, 
    station_info = station_plot_df,
    location_info = location_df, 
    pedigree = pedigree_table_df1,
    trait_data = trait_data_df
  )
  
} # Close the loop



## Read in manually curated 2018 pedigree information
ped_mvn_2018 <- read_excel(file.path(extr_dir, "mvn2018_pedigree_manual.xlsx"), col_types = "text") %>%
  rename_all(tolower) %>%
  ## Extract aliases
  mutate(alias = str_extract(entry, "\\(.*\\)"),
         alias = str_remove_all(alias, "\\(|\\)"),
         entry = str_trim(str_remove(entry, "\\(.*\\)"))) %>%
  # Capitalize all
  mutate_if(is.character, toupper)

## Modify list
extracted_data_list$`MVN 2018.pdf`$pedigree <- ped_mvn_2018



## Read in manually curated 2019 pedigree information
ped_mvn_2019 <- read_excel(file.path(raw_dir, "Agronomic/MISSBAR2019.xlsx"), sheet = "pedigree",
                           col_types = "text") %>%
  rename_all(tolower) %>%
  ## Extract aliases
  mutate(alias = str_extract(entry, "\\(.*\\)"),
         alias = str_remove_all(alias, "\\(|\\)"),
         entry = str_trim(str_remove(entry, "\\(.*\\)"))) %>%
  # Capitalize all
  mutate_if(is.character, toupper)

## Modify list
extracted_data_list$`MVN 2019.pdf`$pedigree <- ped_mvn_2019



### Extract trait data from the 2009 abridged report

pdf_file <- file.path(report_dir, "core tables MVBN2009.pdf")

## Extract text
pdf_txt <- pdf_text(pdf = pdf_file) %>%
  # convert everything to upper-case
  str_to_upper()

# Get coordinates of text
pdf_dat <- pdf_data(pdf = pdf_file) %>%
  map(~mutate(., text = toupper(text)))

# Which is the contents page?
toc_page <- which(str_detect(pdf_txt, "CONTENTS"))



## Get the pages of all tables
table_names <- str_extract(pdf_txt, "TABLE [0-9].*\r") %>%
  str_remove_all(string = ., pattern = "\n|\r")

table_names[toc_page] <- NA


#############################
## Year
#############################

report_year <- "2009"



#############################
## Determine the location names
#############################

## Yield tends to be a trait with universal coverage, so use those pages

# Find the yield pages
yield_pages <- str_which(table_names, "AVE\\. STA\\. YIELD \\(KG/HA\\)")

# Get the data from those pages
yield_pages_table_data <- pdf_dat[yield_pages]
location_name_list <- list()

## Iterate over the pages
for (i in seq_along(yield_pages)) {
  
  ## Extract the data
  data_i <- yield_pages_table_data[[i]]
  # Get the y values of "AVERAGE", since this line will contain the locations
  data_i_sub <- subset(data_i, text == "AVERAGE")
  data_i_sub1 <- subset(data_i, str_detect(text, "STA\\."))
  # Get the coordinates
  dim_use <- unlist(subset(data_i_sub, y == min(y), c(x, y)))
  dim_use1 <- unlist(subset(data_i_sub1, x >= dim_use[1], c(x, y)))
  
  
  ## Data from the location line
  loc_data_i <- subset(data_i, x > dim_use[1] & y == dim_use[2])
  ## Data from the state line
  state_data_i <- subset(data_i, x > dim_use1[1] & y == dim_use1[2]) %>%
    mutate(i = seq(nrow(.)))
  
  ### Use the x values from the state df to determine if location names need to be combined ###
  # First, determine the average width of locations
  avg_log_width <- mean(loc_data_i$width)
  # Are there outliers? (if widths are less than 10 and are adjacent, combine)
  loc_data_i1 <- loc_data_i %>%
    mutate(outlier = width < avg_log_width - 10, i = seq(nrow(.)))
  
  # Combine if the number of outliers is > 1
  
  if (sum(loc_data_i1$outlier) > 1) {
    
    ## Determine groups to combine
    loc_outliers <- filter(loc_data_i1, outlier) %>%
      mutate(diff = c(1, diff(i)), group = i)
    
    for (q in seq(2, length(loc_outliers$i))) {
      group_q <- if (loc_outliers$diff[q] == 1) loc_outliers$group[q - 1] else loc_outliers$i[q]
      loc_outliers$group[q] <- group_q
    }
    
    ## Combine back to the loc_data
    loc_data_i2 <- left_join(loc_data_i1, loc_outliers, by = c("width", "height", "x", "y", "space", "text", "outlier", "i")) %>% 
      mutate(group = ifelse(is.na(group), i, group))
    
    ## Split on group and add back to the df
    location_name_list[[i]] <- sapply(split(loc_data_i2$text, loc_data_i2$group), paste0, collapse = "_")
    
    
  } else {
    # Subset data for those coordinates and get the location names
    location_name_list[[i]] <- subset(data_i, x > dim_use[1] & y == dim_use[2], text, drop = TRUE)
    
  }
  
  
}

## Unlist the location names and assign numbers
location_df <- tibble(location = unlist(location_name_list), environment_number = seq_along(location))





#############################
## Scrape trait data
#############################

# First get the pages containing relevant tables
trait_tables <- str_subset(table_names, "AVE\\. STA\\. [A-Z ]* \\(")
trait_pages <- str_which(table_names, "AVE\\. STA\\. [A-Z ]* \\(")

# A list to store df
trait_data_list <- list()

## A pattern vector to match
col_pattern <- c("NO\\.", "SELECTION", "RANK", "AVERAGE", location_df$location) %>%
  str_replace_all(string = ., "_", " ") %>%
  str_split(string = ., pattern = " ") %>% # Split on spaces
  map_chr(., 1) %>%
  paste0(collapse = "|")


## Iterate over the content df
for (i in seq_along(trait_pages)) {
  
  # Get the page
  page <- trait_pages[i]
  
  # Get the text and coordinates at this page
  page_text <- pdf_txt[[page]]
  page_dat <- pdf_dat[[page]]
  
  
  ## If the page contains "NO DATA", pass
  if (str_detect(page_text, "NO DATA")) next
  
  
  # What is the trait name?
  trait <- str_extract(page_text, "AVE\\. STA\\. [A-Z ]*") %>% 
    str_remove(., "AVE\\. STA\\.") %>% 
    str_trim()
  
  ## Create a table name
  table_name <- paste0(trait, "_page", page)
  
  ## Find the column names and coordinates
  col_coord <- page_dat %>%
    # Remove special characters
    filter(str_detect(text, col_pattern)) %>%
    filter(y < 200) %>%
    arrange(x)
  
  # Determine the max line based on summary values
  last_y <- page_dat %>%
    filter(y > max(col_coord$y),
           str_detect(text, "AVERAGE"),
           y > mean(y)) %>%
    pull(y)
  last_y <- ifelse(length(last_y) == 0, max(page_dat$y), min(last_y))
  
  ## Subset the data to parse
  page_dat_post_col <- page_dat %>%
    filter(y > max(col_coord$y)) %>%
    # Remove the "Average" rows
    filter(y < last_y)
  
  
  col_data_list <- list()
  
  ## List over columns
  for (j in seq(nrow(col_coord))) {
    
    ## Determine the x value for the left hand cut
    ## x1 is left, x2 is right
    # If j == 1, the cut x is the minimum x
    if (j == 1) {
      x1 <- min(col_coord$x)
      
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
  
  ## Add the data to the list
  trait_data_list[[table_name]] <- col_data_df
  
} # End loop


# Combine into a single df
trait_data_df <- trait_data_list %>%
  imap_dfr(~{
    # Determine character columns
    char_col <- summarize_all(.x, is.character) %>% unlist() %>% which()
    
    # Convert data from character to numeric
    mutate_at(.x, vars(setdiff(char_col, 2)), parse_number) %>%
      mutate(table_name = .y, trait = str_remove(.y, "_[a-z0-9]*$")) %>% 
      ## Rename "entry" as "number"
      rename_all(~str_replace_all(string = ., pattern = "Entry", replacement = "Number"))
    
  }) %>%
  rename(number = `NO.`, line_name = SELECTION) %>%
  select(-AVERAGE, -RANK) %>%
  gather(location, value, -1, -2, -table_name, -trait) %>%
  filter(!is.na(value))


## Pedigree information will simply be taken from a trait df
pedigree_table_df1 <- trait_data_df %>%
  # Remove trailing numbers as product of poor column separation
  mutate(line_name = str_replace(line_name, " [0-9]{1,2}$", "") %>% str_replace("(\\))([0-9]{1,2}$)", "\\1")) %>%
  distinct(line_name) %>%
  select(entry = line_name) %>%
  ## Extract aliases
  mutate(alias = str_extract(entry, "\\(.*\\)"),
         alias = str_remove_all(alias, "\\(|\\)"),
         entry = str_trim(str_remove(entry, "\\(.*\\)"))) %>%
  # Capitalize all
  mutate_if(is.character, toupper)


## Add this information to the list
extracted_data_list[[basename(pdf_file)]] <- list(
  filename = basename(pdf_file), 
  year = report_year, 
  location_info = location_df, 
  pedigree = pedigree_table_df1,
  trait_data = trait_data_df
)


### Remove trailing numbers from line name, if present
extracted_data_list1 <- extracted_data_list %>%
  map(., ~modify_at(., c("trait_data"), ~mutate(.x, line_name = str_replace(line_name, " [0-9]{1,2}$", "") %>% 
                                                  str_replace("(\\))([0-9]{1,2}$)", "\\1")) ) )




## 2017 cannot be resolved. I cannot match locations to the environment
## number in the excel file, so this year for this nursery is skipped.



## Save this data
mvn_extracted_data_list <- extracted_data_list1

save("mvn_extracted_data_list", file = file.path(extr_dir, "mvn_extracted_metadata.RData"))









#############################
## Clean up the results
#############################


# Load the data
load(file.path(extr_dir, "mvn_extracted_metadata.RData"))



#### Pedigree


## Pull the pedigree data and assemble into a single df
pedigree_df <- map(mvn_extracted_data_list, "pedigree") %>% 
  imap_dfr(~mutate(.x, file = .y)) %>% 
  select(entry, program, pedigree, alias, file) %>% 
  mutate(report = str_remove(file, ".pdf") %>% str_replace(., " ", "_")) %>%
  select(-file) %>%
  as_tibble() %>%
  mutate_if(is.character, parse_guess) %>%
  filter(!is.na(entry))

## Save csv
write_csv(x = pedigree_df, path = file.path(extr_dir, "mvn_extracted_entry_metadata.csv"), na = "")





#### Station information

## Load location information
location_info_df <- map(mvn_extracted_data_list, "location_info") %>%
  imap_dfr(~mutate(.x, file = .y)) %>% 
  select(environment_number, location, file) %>% 
  mutate(report = str_remove(file, ".pdf") %>% str_replace(., " ", "_"),
         location = str_replace_all(location, "_", " ") %>% str_to_title()) %>%
  select(-file) %>%
  as_tibble()

## Output distinct location/state for corrections
location_info_df %>% 
  distinct(location) %>%
  arrange(location) %>%
  write_csv(x = ., path = file.path(extr_dir, "mvn_extracted_location_rename.csv"))

## Custom location renaming df
location_rename <- read_excel(path = file.path(extr_dir, "mvn_extracted_location_rename.xlsx"))


# Rename
location_info_df1 <- location_info_df %>% 
  left_join(., location_rename) %>%
  select(environment_number, location = rename, report)


## Pull the station info data and assemble into a single df
station_info_df <- map(mvn_extracted_data_list, "station_info") %>% 
  imap_dfr(~mutate(.x, file = .y)) %>% 
  select(number, location, planting_date, harvest_date, file) %>% 
  mutate(report = str_remove(file, ".pdf") %>% str_replace(., " ", "_"),
         location = str_to_title(location)) %>%
  select(-file) %>%
  as_tibble() %>%
  mutate_at(vars(contains("date")), ymd)

## Add the planting/harvest dates (where applicable) to the location information df
location_info_df2 <- location_info_df1 %>%
  left_join(., station_info_df) %>%
  select(report, location, environment_number, contains("date"))

## Save csv
write_csv(x = location_info_df2, path = file.path(extr_dir, "mvn_extracted_station_info.csv"))










