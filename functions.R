## Relevant functions for this analysis
## 
## 

# A function to join two df using a numeric given some tolerance
coord_join <- function(x, y, by, tol, remove.NA = FALSE) {
  
  # Make sure x and y are dfs
  stopifnot(is.data.frame(x))
  stopifnot(is.data.frame(y))
  
  # By must be a character that is a column shared by x and y
  stopifnot(is.character(by))
  # By must be length 1 and positive
  stopifnot(length(by) == 1)
  stopifnot(by >= 0)
  
  if (! by %in% intersect(names(x), names(y))) stop (paste0("Cannot join on '", by, "' because it is not a shared column."))
  
  # By must target a numeric column
  if ( any(!is.numeric(x[[by]]), !is.numeric(y[[by]])) ) stop(paste0("Column ", by, " must be numeric in both x and y."))
  
  ## Remove rows in x or y that contain NA in the by column
  x1 <- x[!is.na(x[[by]]),,drop = FALSE]
  y1 <- y[!is.na(y[[by]]),,drop = FALSE]
  
  # For all elements of by in y, find those matching in x within the specified tolerance
  by1_list <- sapply(X = y1[[by]], FUN = function(r) subset(x1[[by]], x1[[by]] <= r + tol & x1[[by]] >= r - tol), simplify = FALSE)
  
  # Replace empty elements with NA
  by1_list[sapply(X = by1_list, is_empty)] <- NA
  y1[[by]] <- unlist(by1_list)
  
  ## Remove the NAs - if called
  if (remove.NA) {
    y1 <- y1[!is.na(y1[[by]]),,drop = FALSE]
  }
  
  
  # Join and return
  full_join(x = x1, y = y1, by = by)
  
}










## A function to read in malting quality data from the formatted excel sheet
read_maltq <- function (file, sheets, tidy = FALSE) {
  
  if (!all(endsWith(x = file, suffix = c(".xlsx")) | endsWith(x = file, 
                                                              suffix = c(".xls")))) 
    stop("The filenames in 'file' do not have the extension '.xls' or '.xlsx'.")
  
  maltq.compiled <- list()
  
  ## If sheets is missing, just take the first sheet
  sheets <- if (missing(sheets)) 1L else sheets
  sheet_names <- if (missing(sheets)) { "sheet1" } else { sheets }
  
  for (s in sheets) {
    print(s)
    
    data.f <- read_excel(path = file, sheet = s, col_names = FALSE) %>% 
      filter_all(any_vars(!is.na(.))) %>%
      select_if(~!all(is.na(.)))
    
    # Table number
    table_col <- map_lgl(data.f, ~any(str_detect(., "Table"), na.rm = TRUE))
    table_row <- str_which(string = data.f[,table_col,drop = T], pattern = "^Table[ 0-9]*")[1]
    table.no <- data.f[table_row,table_col,drop = T]
    data.f1 <-  slice(data.f, -seq(1, table_row))
    var.names <- slice(data.f1, 1:3) %>% apply(X = ., MARGIN = 2, 
                                               FUN = function(var) paste0(na.omit(var), collapse = "")) %>% 
      str_replace_all(pattern = " ", replacement = "") %>% 
      str_replace_all(pattern = "\"", replacement = "") %>% 
      as.character()
    data.f2 <- slice(data.f1, -c(1:3))
    names(data.f2) <- make.names(var.names, unique = TRUE)
    lab_no_col <- which(map_lgl(data.f, ~any(str_detect(., "LabNo|Lab No|Lab \\#"), na.rm = TRUE)))
    
    table.row <- which(data.f2[[lab_no_col]] == table.no)
    if (length(table.row) != 0) {
      to.remove <- as.numeric(sapply(X = table.row, FUN = function(row) seq(row, row + 3)))
      data.f2 <- slice(data.f2, -to.remove)
    }
    # The first column name is always LabNo.
    names(data.f2)[1] <- "LabNo."
    
    coef.row <- which(select(data.f2, LabNo.) == "Coefficients of Variation")
    data.f3 <- slice(data.f2, -((coef.row + 1):nrow(data.f2)))
    data.f4 <- suppressWarnings({
      data.f3 %>% mutate_at(vars(-LabNo., -VarietyorSelection), 
                            .funs = parse_number)
    }) %>% slice(which(rowMeans(is.na(.)) != 1)) %>% mutate(Batch = table.no, sheet = s)
    data.check <- data.f4 %>% filter(str_detect(VarietyorSelection, 
                                                "MALT CHECK"))
    data.stats <- data.f4 %>% filter(str_detect(LabNo., "^[^0-9]{1,}"))
    data.f5 <- dplyr::setdiff(data.f4, bind_rows(data.stats, 
                                                 data.check))
    maltq.compiled[[s]] <- list(data = data.f5, checks = data.check, 
                                stats = data.stats)
  }
  maltq.data <- maltq.compiled %>% map(function(table_list) table_list$data) %>% 
    bind_rows() %>% arrange(VarietyorSelection)
  maltq.stats <- maltq.compiled %>% map(function(table_list) table_list$stats) %>% 
    bind_rows()
  maltq.checks <- maltq.compiled %>% map(function(table_list) table_list$checks) %>% 
    bind_rows()
  df_names <- c("maltq.data", "maltq.stats", "maltq.checks")
  df_list <- list(maltq.data, maltq.stats, maltq.checks)
  names(df_list) <- df_names
  if (tidy) {
    df_list <- df_list %>% map(gather, key = "Parameter", 
                               value = "value", -VarietyorSelection, -Batch, -LabNo.)
    return(df_list)
  }
  else {
    return(df_list)
  }
}
  

## Internal function to reconstruct a pedigree using math notation
reconstruct_pedigree <- function(x) {
  
  # rename
  .x <- x
  
  # detect bc
  detect_bc <- str_detect(string = .x, pattern = "[0-9]*\\*|\\*[0-9]*")
  
  # If the length of the vector is != 1, do a check for backcrossing
  if (length(.x) == 1) {
    return(.x[[1]])
    
  } else if (any(detect_bc)) {
    ent <- .x[[which(detect_bc)]]
    
    ## If length is == 1, proceed normally
    if (length(ent) == 1) {
      
      ## prepare the backcross string
      bc_detect_pattern <- "[0-9]{1,2}\\*$|\\*[0-9]{1,2}$|^[0-9]{1,2}\\*|^\\*[0-9]{1,2}"
      ## Find the backcross parent
      bc_parent <- str_remove(string = ent, pattern = bc_detect_pattern)
      # number of backcrosses
      n_bc <- as.numeric(str_remove(string = str_extract(string = ent, pattern = bc_detect_pattern), "\\*"))
      
      # Left or right parens
      if (which(detect_bc) == 1) {
        ent_new <- paste0("(", rep(bc_parent, n_bc), "/", collapse = "")
        # Replace
        .x[[which(detect_bc)]] <- ent_new
        ct1 <- paste0("(", paste0(.x, collapse = ""), paste0(rep(")", n_bc + 1), collapse = ""))
        
      } else {
        ent_new <- paste0("/", rep(bc_parent, n_bc), ")", collapse = "")
        # Replace
        .x[[which(detect_bc)]] <- ent_new
        ct1 <- paste0(paste0(rep("(", n_bc + 1), collapse = ""), paste0(.x, collapse = ""), ")")
      }
      
      return(ct1)
      
      # If length > 1, proceed differently 
    } else {
      # Detect again
      which_ent <- str_which(string = ent, pattern = bc_detect_pattern)
      ent1 <- ent[[which_ent]]
      
      bc_parent <- str_remove(string = ent1, pattern = bc_detect_pattern)
      # number of backcrosses
      n_bc <- as.numeric(str_remove(string = str_extract(string = ent1, pattern = bc_detect_pattern), "\\*"))
      
      ## Combine ent, then backcross with n_bc - 1
      ent2 <- paste0("(", setdiff(ent, ent1), "/", bc_parent, ")", collapse = "/")
      ent_new <- paste0(paste0(rep("(", n_bc - 1)), ent2, paste0("/", rep(bc_parent, n_bc - 1), ")"))
      
      # Replace in .x
      .x[[which(detect_bc)]] <- ent_new
      return(.x)
      
    }
    
  } else {
    return(paste0("(", paste0(.x, collapse = "/"), ")"))
  }
}





## Write a function to parse a pedigree and convert it to mathematical notation
purdy2math <- function(x) {

  purdy <- x
  
  # if NA, return NA
  if (is.na(purdy)) return(NA)
  
  # Stop if the length is > 1
  stopifnot(length(purdy) == 1)

  
  ## Start the parsing process
  
  ## Extract the unique cross notations
  cross_notation <- unique(str_extract_all(string = purdy, pattern = "/[0-9]*/|/")[[1]])
  # What is the number of the most basal cross?
  # For / and //, simply count the number of /
  cross_number <- as.numeric(str_extract(string = cross_notation, pattern = "[0-9]{1,}"))
  cross_number[is.na(cross_number)] <- str_count(string = cross_notation, pattern = "/")[is.na(cross_number)]
  
  ## Iterate over the sorted unique cross numbers
  cross_number_sorted <- sort(x = cross_number, decreasing = TRUE)
  purdy_n <- purdy
  
  for (i in seq_along(cross_number_sorted)) {
    
    n <- cross_number_sorted[i]
    # Create the cross notation
    cross_notation_i <- ifelse(n > 2, paste0("/", n, "/"), paste0(rep("/", n), collapse = ""))
    
    ## Split the pedigree on this notation
    ped_split_n <- map_depth(.x = purdy_n, .depth = vec_depth(purdy_n), .f = ~str_split(string = ., pattern = cross_notation_i))
    ped_split_n <- modify_depth(.x = ped_split_n, .depth = vec_depth(ped_split_n) - 2, .f = unlist)
    
    ## Purdy n is the next level in the list
    purdy_n <- ped_split_n
    
  }
  
  ## Find any backcrosses and replace
  ped_list1 <- purdy_n[[1]]
  
  ## Reconstruct the pedigree from the bottom-up
  for (d in rev(seq_len(vec_depth(ped_list1) - 1))) {
    ped_list1 <- modify_depth(.x = ped_list1, .depth = d, .f = reconstruct_pedigree)
    
  }
  
  ## Collapse the last level to finalize
  new_purdy <- reconstruct_pedigree(ped_list1) %>%
    structure(., class = c(class(.), "pedigree.arithmetic"))
  
  ## Return this pedigree
  return(new_purdy)
  
}

## A function to convert an arithmetic pedigree string into 
## a nested list
as.list.pedigree.arithmetic <- function(x, individual = NULL) {
  
  # Replace "(" with "list("
  x1 <- str_replace_all(x, "\\(", "list\\(")
  # Replace / with ,
  x2 <- str_replace_all(string = x1, pattern = "/", replacement = ", ")

  # Get the individuals from the arithmetic string
  pedigree_indv <- str_replace_all(string = x, pattern = "\\(|\\)|/", replacement = "  ") %>%
    str_trim() %>%
    str_split(string = ., pattern = " ") %>% 
    first() %>% 
    subset(x = ., . != "") %>%
    unique()

  # Add quotes
  quoted_indiv <- paste0("'", pedigree_indv, "'") %>%
    # Add literal indicator for punctuation
    set_names(., paste0(str_replace_all(string = pedigree_indv, pattern = "([:punct:])", "\\\\\\1"), "\\b"))
  # Rename if punctuation is present
  contains_punct <- str_detect(string = names(quoted_indiv), pattern = "\\.")
  names(quoted_indiv)[contains_punct] <- str_remove(string = names(quoted_indiv)[contains_punct], 
                                                    pattern = "\\\\b")
    
  
  
  # Replace in x2
  x3 <- str_replace_all(string = x2, quoted_indiv)
  # Replace inner crosses with c instead of list
  x4 <- str_replace_all(string = x3, 
                        pattern = "(list)(\\(\\'[A-Za-z0-9.-]*\\'\\, \\'[A-Za-z0-9.-]*\\'\\))", 
                        replacement = "c\\2")
  
  # Parse and evaluatate the list
  x4_parse <- eval(parse(text = x4))
  # Wrap in list if individual name is provided
  if (!is.null(individual)) {
    x4_list <- list(x4_parse)
    names(x4_list) <- individual
  } else {
    x4_list <- x4_parse
  }
  class(x4_list) <- union(class(x4_list), "pedigree.list")

  return(x4_list)
}


## A function to collapse a pedigree list into a pedigree string; either
## arithmetic or purdy
collapse.pedigree.list <- function(x) {
  
  # Working from the deepest part of the list to the shallowest, collapse the pedigree
  ped_list1 <- x

  # What is the pedigree list depth; subtract 1
  ped_depth <- vec_depth(ped_list1) - 1
  ped_depth_seq <- rev(seq_len(ped_depth))
  
  ## Reconstruct the pedigree from the bottom-up
  for (d in ped_depth_seq) {
    
    # Recursively collapse the pedigree
    ped_list1 <- map_depth(.x = ped_list1, .depth = d, .ragged = TRUE, .f = ~{
      if (length(.x) == 1) { .x[[1]] } else { paste0("(", paste0(.x, collapse = "/"), ")") } })

  }

  out <- ped_list1[[1]]
  
  ## Assign class
  class(out) <- union(class(out), "pedigree.arithmetic")
  return(out)
    
  # paste0("(", paste0(ped_list1[[1]], collapse = "/"), ")")
  # ## Collapse the last level without parens
  # paste0(, collapse = "/")
  
}





## A function to expand the genealogy of a pedigree by referencing a
## database
## 
## It requires the following:
## 1. sample line_name
## 2. sample pedigree
## 3. sample aliases
## 4. reference (line_name, pedigree, aliases)
## 

expand_pedigree <- function(line_name, pedigree, aliases, reference) {
  
  ## Check input classes
  stopifnot(is.character(line_name))
  stopifnot(is.character(pedigree))
  stopifnot(is.character(aliases))
  stopifnot(is.data.frame(reference))
  
  # Ensure column names and types in reference are correct
  correct_col_names <- c("line_name", "pedigree")
  correct_col_classes <- setNames(c("character", "character"), correct_col_names)
  if (!all(correct_col_names %in% names(reference)))
    stop("The 'reference' input must have the columns 'line_name', 'pedigree', and 'aliases.'")
  
  if (!all(sapply(reference[correct_col_names], class) == correct_col_classes))
    stop(paste0("Columns in reference must have the following classes: ", 
                paste(names(correct_col_classes), correct_col_classes, sep = ": ", collapse = ", ")))
  
  ## There should be no NA line names or NA pedigrees
  if (any(map_dbl(reference[correct_col_names], ~sum(is.na(.))) > 0)) {
    stop ("There should be no missing data in the 'reference' input.")
  }
  
  
  
  
  ## Here are the steps of the algorithm:
  ## 1. Convert the pedigree to arithmetic
  ## 2. Extract individuals in the pedigree
  ## 3. For each individual; perform a look up by line name (first) or alias (second) in
  ## the reference input
  ## 4. Extract the pedigree
  ## 5. Repeat 1-4 until the pedigree cannot be expanded further.
  ## 6. Repeat 1-5 for each individual in the original pedigree
  
  
  ##### All this below will be a new function ####
  
  pedigree_to_check <- list(pedigree)
  names(pedigree_to_check) <- line_name
  # Is it a pedigree?
  any_ped <- any(grepl(pattern = "/", x = unlist(x = pedigree_to_check, recursive = TRUE)), na.rm = TRUE)
  
  ## A list to store extra individual names as the loop dives into the pedigree
  accessory_individuals <- list(character(0))
  r = 1
  
  # While loop
  while (any_ped) {
    
    # Convert pedigree to arithmetic
    # pedigree_arith <- rapply(object = pedigree_to_check, f = purdy2math, how = "replace")
    pedigree_arith <- pedigree_to_check %>%
      modify_depth(.x = ., .depth = vec_depth(.), .ragged = TRUE, .f = ~{
        modify_if(.x = ., .p = ~any(grepl(pattern = "/", x = .), na.rm = TRUE), purdy2math)
      })

    # Convert to list
    # pedigree_list <- rapply(object = pedigree_arith, how = "replace", f = function(x) {
    #   x1 <- as.list.pedigree.arithmetic(x)
    #   structure(x1, class = NULL)
    # })
    pedigree_list <- pedigree_arith %>%
      map_depth(.x = ., .depth = vec_depth(.) - 1, .ragged = TRUE, .f = ~{
        map(.x = .x, as.list.pedigree.arithmetic)
      })
    
    ## Unlist the pedigree and add all names to the accessory list
    accessory_individuals[[r]] <- unlist(pedigree_list, recursive = TRUE)
    
    ## Modify only elements that are not pedigree lists
    pedigree_list1 <- pedigree_list %>%
      modify_depth(.x = ., .depth = vec_depth(.) - 1, .ragged = TRUE, .f = ~{
        # Only modify if not a pedigree list
        l1 <- .x
        modify_if(.x = l1, .p = ~!inherits(., "pedigree.list"), .f = ~{
          
          # Look up this individual in "reference"
          indiv_ref <- head(subset(reference, line_name == .x), 1)
          # Should simply the individual be returned?
          flag <- any((nrow(indiv_ref) == 0), (grepl(pattern = "UNKNOWN|=", x = indiv_ref$pedigree)))

          if (flag) {
            .x
          } else {
            indiv_ref$pedigree
          }
          
        })
      })
    
    # Reassign
    pedigree_to_check <- pedigree_list1
    
    # Are any of the elements a pedigree?
    any_ped <- any(grepl(pattern = "/", x = unlist(x = pedigree_to_check, recursive = TRUE)), na.rm = TRUE)
    r = r + 1
    
  } # CLose the loop
  
  
  ## Collapse the large pedigree string
  ped_collapse <- collapse.pedigree.list(pedigree_to_check)
  # Reduce all of the accessory individuals
  accessory_individuals1 <- reduce(accessory_individuals, union)
  
  # Return a data.frame
  tibble(line_name = line_name, pedigree = pedigree, expanded_pedigree = ped_collapse,
         accessory = list(accessory_individuals1))
  
}



## A different version of the expand pedigree function
## 
## 
expand_pedigree2 <- function(line_name, pedigree, reference) {
  
  ## Check input classes
  stopifnot(is.character(line_name))
  stopifnot(is.character(pedigree))
  stopifnot(is.data.frame(reference))
  
  # Ensure column names and types in reference are correct
  correct_col_names <- c("line_name", "pedigree")
  correct_col_classes <- setNames(c("character", "character"), correct_col_names)
  if (!all(correct_col_names %in% names(reference)))
    stop("The 'reference' input must have the columns 'line_name' and 'pedigree'")
  
  if (!all(sapply(reference[correct_col_names], class) == correct_col_classes))
    stop(paste0("Columns in reference must have the following classes: ", 
                paste(names(correct_col_classes), correct_col_classes, sep = ": ", collapse = ", ")))
  
  ## There should be no NA line names or NA pedigrees
  if (any(map_dbl(reference[correct_col_names], ~sum(is.na(.))) > 0)) {
    stop ("There should be no missing data in the 'reference' input.")
  }
  
  
  
  
  ## Here are the steps of the algorithm:
  ## 1. Convert the pedigree to arithmetic
  ## 2. Extract individuals in the pedigree
  ## 3. For each individual; perform a look up by line name (first) or alias (second) in
  ## the reference input
  ## 4. Extract the pedigree
  ## 5. Repeat 1-4 until the pedigree cannot be expanded further.
  ## 6. Repeat 1-5 for each individual in the original pedigree
  
  
  pedigree_to_check <- list(pedigree)
  names(pedigree_to_check) <- line_name
  # Is it a pedigree?
  any_ped <- any(grepl(pattern = "/", x = unlist(x = pedigree_to_check, recursive = TRUE)), na.rm = TRUE)
  
  ## A list to store extra individual names as the loop dives into the pedigree
  accessory_individuals <- list(character(0))
  r = 1
  
  # While loop
  while (any_ped) {
    
    # Convert pedigree to arithmetic
    pedigree_arith <- map_if(.x = pedigree_to_check, .p = ~any(grepl(pattern = "/", x = .)), .f = purdy2math)
  
    # Convert to a list; identify individuals in the pedigree string
    pedigree_list <- pedigree_arith %>%
      map(as.list.pedigree.arithmetic)
    
    ## Unlist the pedigree and add all names to the accessory list
    accessory_individuals[[r]] <- accessory1 <- unique(unlist(pedigree_list, recursive = TRUE))
    
    # For each accessory individual, look them up in the reference set
    accessory1_lookup <- accessory1 %>%
      map(~{
        # Look up this individual in "reference"
        indiv_ref <- head(subset(reference, line_name == .x), 1)
        # Should simply the individual be returned?
        flag <- any((nrow(indiv_ref) == 0), (grepl(pattern = "UNKNOWN|=", x = indiv_ref$pedigree)))
        
        if (flag) {
          .x
        } else {
          indiv_ref$pedigree
        }
      })
    
    # Reassign
    pedigree_to_check <- accessory1_lookup
    
    # Are any of the elements a pedigree?
    any_ped <- any(grepl(pattern = "/", x = unlist(x = pedigree_to_check, recursive = TRUE)), na.rm = TRUE)
    r = r + 1
    
  } # CLose the loop
  
  
  # Reduce all of the accessory individuals
  accessory_individuals1 <- reduce(accessory_individuals, union)
  
  # Return a data.frame
  tibble(line_name = line_name, pedigree = pedigree,
         accessory = list(accessory_individuals1))
  
}


# Create a function to merge entry information iterate over line_name
merge_entry_info <- function(x) {
  data <- x
  
  # Pull out alias, reduce
  alias_reduced <- reduce(data$aliases, union) %>%
    na.omit() %>%
    as.character()
  .x <- distinct(select(data, -aliases))
  
  # For each line, remove some rows:
  # If a column is incompletely NA, filter out NA
  # Priortize malt, feed, food
  
  ## 
  # Are there any columns with missing data > 0 & < 1?
  filter_cols <- map_dbl(.x, ~mean(is.na(.))) %>% 
    subset(., . < 1 & . > 0) %>%
    names()
  
  if (!is_empty(filter_cols)) {
    # Coalesce these columns
    .x2 <- .x %>%
      mutate(.id = paste0("id", seq(nrow(.)))) %>% 
      gather(col, value, -.id) %>% 
      spread(.id, value) %>% 
      mutate(final = coalesce(!!!select(., -col))) %>%
      select(col, final) %>%
      spread(col, final)
    
  } else {
    .x2 <- .x
  }
  
  # Find unique grades, if all are not NA
  if (!all(is.na(.x2$grade))) {
    .x3 <- mutate(.x2, grade = factor(grade, levels = grade_priority), grade_no = as.numeric(grade)) %>%
      group_by_at(vars(-grade, -grade_no)) %>% 
      top_n(x = ., n = 1, wt = -grade_no) %>% 
      ungroup() %>% 
      select(-grade_no)
    
  } else {
    .x3 <- .x2
  }
  
  # Add alias back in
  mutate(.x3, aliases = list(alias_reduced), row_type = as.numeric(row_type))
}



## Function to check a pedigree
check_pedigree <- function(x) {
  
  if (is.na(x)) return(NA)
  
  ## Are there any slashes?
  f1 <- grepl(pattern = "/", x = x)
  if (!f1) return(f1)
  
  # Check if the correct order of slashes is present
  slashes <- str_extract_all(string = x, pattern = "/[0-9]{1,}/|//|/")[[1]]
  # Highest cross number
  cross_numbers <- as.numeric(suppressWarnings(parse_number(slashes)))
  cross_numbers[is.na(cross_numbers)] <- str_count(string = slashes[is.na(cross_numbers)], pattern = "/")
  
  # Highest
  max_cross <- max(cross_numbers)
  # Are all preceding crosses present?
  f2 <- all(seq_len(max_cross) %in% cross_numbers)
    
  if (!f2) return(f2)
  
  # You can't have more of a higher cross number than a lower cross number
  if (max_cross >= 2) {
    f3 <- map_lgl(seq(max_cross, 2), ~sum(cross_numbers == .x) > sum(cross_numbers == (.x - 1)))
    if (any(f3)) return(FALSE)
  }
  
  # A pedigree cannot contain solely >= 2 of a cross number
  f4 <- length(unique(cross_numbers)) == 1 & length(cross_numbers) >= 2
  if (f4) return(!f4)
  
  # You cannot have identical sequential cross numbers
  f5 <- any(diff(cross_numbers) == 0)
  if (f5) return(!f5)
      
  ## Finally print TRUE
  return(TRUE)
  
}







## Calculate heritability based on sommer mixed model output
herit.mmer <- function(x, gen.var = "u:line_name", method = c("Cullis")) {
  
  # Match method
  method <- match.arg(method)
  
  # Extract common components
  varG <- as.numeric(x$sigma[[gen.var]]) # Genetic variance component
  n_g <- n_distinct(x$data$line_name) # Number of genotypes
  C22 <- x$PevU[[gen.var]]$value
  
  # Direct by method
  if (method == "Cullis") {
    trC22 <- sum(diag(C22))
    # Mean variance of a difference between BLUPs
    vdBLUP <- 2 / n_g * (trC22 - (sum(C22) - trC22)/(n_g - 1))
    
    H2 <- 1 - (vdBLUP / 2 / varG)
    
  } else {
    stop("!")
    
  }
  
  # Return heritability
  return(H2)
  
}





## Function to perform the logit transformation
logitTransform <- function(x, tol = 0.001) { 
  # Values to add to 0 or 1
  x[x == 0] <- x[x == 0] + tol
  x[x == 1] <- x[x == 1] - tol
  log(x/(1-x))
}

# Exponentiate from logit transformation
logitExp <- function(x) exp(x) / (1 + exp(x))


# function for calculating angle
angle <- function(v1, v2, degrees = TRUE) {
  rad <- acos( sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2))) )
  ifelse(degrees, rad * (180/pi), rad)
}






  


