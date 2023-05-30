#____________________________ 
# ┬  ┬┌┐ ┬─┐┌─┐┬─┐┬┌─┐┌─┐
# │  │├┴┐├┬┘├─┤├┬┘│├┤ └─┐
# ┴─┘┴└─┘┴└─┴ ┴┴└─┴└─┘└─┘
#__________________________

library(here)
library(dplyr)
library(stringr)
library(gt)

#______________________________________________________________________
# ╔═╗╦ ╦╔╗╔╔═╗╔╦╗╦╔═╗╔╗╔  ╔═╗╔╗╔╔═╗       ┬─┐┌─┐┌─┐┌┬┐  ┌─┐┬┬  ┌─┐┌─┐
# ╠╣ ║ ║║║║║   ║ ║║ ║║║║  ║ ║║║║║╣   ───  ├┬┘├┤ ├─┤ ││  ├┤ ││  ├┤ └─┐
# ╚  ╚═╝╝╚╝╚═╝ ╩ ╩╚═╝╝╚╝  ╚═╝╝╚╝╚═╝       ┴└─└─┘┴ ┴─┴┘  └  ┴┴─┘└─┘└─┘
#______________________________________________________________________

files=list.files(here("data","barbac_testdata"), pattern = "*barcodes.txt$", full.names = TRUE)   # select all the files that have isolate barcodes 
files         # all files

files[1]      # only one file

raw_counts <- function(file_path) {
  # Read the table from the file
  data <- tryCatch(
    {
      read.table(file_path, sep = ",")
    },
    error = function(e) {
      stop("Error: Failed to read the table from the file. Check if the file exists and the path is correct.")
    }
  )
  
  # Remove dashes from the 'barcode' column
  data <- tryCatch(
    {
      data %>% mutate(barcode = str_replace_all(barcode, "-", ""))
    },
    error = function(e) {
      stop("Error: Failed to remove dashes from the 'barcode' column. Check if the column exists and contains valid data.")
    }
  )
  
  # Add a new column 'length_barcode' with the length of 'barcode'
  data <- tryCatch(
    {
      data %>% mutate(length_barcode = str_length(barcode))
    },
    error = function(e) {
      stop("Error: Failed to calculate the length of 'barcode'. Check if the 'barcode' column exists and contains valid data.")
    }
  )
  
  # Filter out rows with a 'length_barcode' of "0"
  data <- tryCatch(
    {
      data %>% filter(length_barcode != "0") %>% as_tibble()
    },
    error = function(e) {
      stop("Error: Failed to filter rows with a 'length_barcode' of '0'. Check if the 'length_barcode' column exists and contains valid data.")
    }
  )
  
  # Return the resulting data
  return(data)
}
test <- raw_counts(files[1])

#_______________________________________________________________________________
# ╔═╗╦ ╦╔╗╔╔═╗╔╦╗╦╔═╗╔╗╔  ╔═╗╔╗╔╔═╗ ╔═╗╔╗╔╔═╗       ┬─┐┌─┐┌─┐┌┬┐  ┌─┐┬┬  ┌─┐┌─┐
# ╠╣ ║ ║║║║║   ║ ║║ ║║║║  ║ ║║║║║╣  ║ ║║║║║╣   ───  ├┬┘├┤ ├─┤ ││  ├┤ ││  ├┤ └─┐
# ╚  ╚═╝╝╚╝╚═╝ ╩ ╩╚═╝╝╚╝  ╚═╝╝╚╝╚═╝o╚═╝╝╚╝╚═╝       ┴└─└─┘┴ ┴─┴┘  └  ┴┴─┘└─┘└─┘
#_______________________________________________________________________________

library(dplyr)
library(gt)
library(gtExtras)

table <- function(file, barcode_length) {
  # Check if the 'barcode_length' parameter is a numeric vector of length 2
  if (!is.numeric(barcode_length) || length(barcode_length) != 2) {
    stop("Error: 'barcode_length' should be a numeric vector of length 2.")
  }
  
  # Check if 'file' is a data frame or tibble
  if (!is.data.frame(file) && !is.tbl(file)) {
    stop("Error: 'file' should be a data frame or tibble.")
  }
  
  file %>%
    mutate(
      bin = case_when(
        (length_barcode < barcode_length[1] & length_barcode > 0) ~ paste0("(", 0, ",", barcode_length[1], ")"),
        (length_barcode >= barcode_length[1] & length_barcode <= barcode_length[2]) ~ paste0("[", barcode_length[1], ",", barcode_length[2], "]"),
        TRUE ~ paste0("[", barcode_length[2] + 1, ",∞)")
      )
    ) %>%
    group_by(bin) %>%
    count() %>%
    tableGrob()
    # gt() %>%
    # cols_label(n = "No of barcodes") %>%
    # tab_source_note(
    #   source_note = "Barbac: A versatile tool to quantify barcodes."
    # ) %>%
    # gt_theme_espn()
}

# Usage example
test
file <- test  # Assuming 'data' is your input dataframe
barcode_length <- c(22, 28)  # Specify the bin ranges here

table(file, barcode_length)

test %>% 
  ggplot(., aes(x=length_barcode)) +
  geom_histogram(fill="#FF5349", color=NA, alpha=0.8) +
  labs(x="Barcode length (bp)", y="Number of barcodes", title="Barcode length distribution") + 
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))


test %>%   
ggplot(., aes(x=log(counts))) +
geom_histogram(fill="#FF5349", color=NA, alpha=0.8) +
  labs(x="log(Barcode counts)", y="Number of barcodes", title="Barcode counts distribution") + 
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

# Function to calculate Shannon entropy
shannon_entropy <- function(seq) {
  n <- nchar(seq)
  nucleotides <- unique(strsplit(seq, "")[[1]])
  counts <- sapply(nucleotides, function(nuc) sum(strsplit(seq, "")[[1]] == nuc))
  probabilities <- counts / n
  entropy <- -sum(probabilities * log2(probabilities))
  return(entropy)
}

test %>% 
  mutate(entropy = sapply(barcode, shannon_entropy)) %>% 
  ggplot() +
  geom_histogram(aes(x=entropy)) +
  theme_b

#__________________________________________________________________________ 
# ╔═╗╦ ╦╔╗╔╔═╗╔╦╗╦╔═╗╔╗╔  ╔╦╗╦ ╦╔═╗       ┌─┐┬  ┬ ┬┌─┐┌┬┐┌─┐┬─┐┬┌┐┌┌─┐
# ╠╣ ║ ║║║║║   ║ ║║ ║║║║   ║ ║║║║ ║  ───  │  │  │ │└─┐ │ ├┤ ├┬┘│││││ ┬
# ╚  ╚═╝╝╚╝╚═╝ ╩ ╩╚═╝╝╚╝   ╩ ╚╩╝╚═╝       └─┘┴─┘└─┘└─┘ ┴ └─┘┴└─┴┘└┘└─┘
#__________________________________________________________________________

library(dplyr)
library(stringdist)
library(igraph)

df <- test
group.compmat <- function(df, method = "lv", distance) {
  
  # Check if the input dataframe has the required columns
  # Check if 'df' is a data frame or tibble
  if (!is.data.frame(df) && !is(tbl_df(df))) {
    stop("Input 'df' must be a data frame or tibble.")
  }
  
  # Check if 'barcode' column exists in 'df'
  if (!"barcode" %in% colnames(df)) {
    stop("Column 'barcode' not found in the input data frame.")
  }
  
  # Check if 'counts' column exists in 'df'
  if (!"counts" %in% colnames(df)) {
    stop("Column 'counts' not found in the input data frame.")
  }
  
  # Check if 'distance' is numeric
  if (!is.numeric(distance)) {
    stop("Parameter 'distance' must be numeric.")
  }
  
  df <- df %>% mutate(id = row_number())
  
  df1 <- df %>% dplyr::select(1)
  
  my_dist <- function(x, y, ...) {
    1 - stringdist(x, y, method = method)
  }
  
  p <- pair_minsim(df1, on = "barcode", deduplication = TRUE,
                   default_comparator = my_dist, minsim = 1 - distance)
  
  g <- graph.data.frame(p, directed = TRUE)
  g1 <- p %>% mutate(col3 = str_c("group", clusters(g)$membership[as.character(.x)]))
  
  deduplicate_equivalence(p, "group")
  
  to_add <- data.frame(.x = seq_len(nrow(df1)), .y = seq_len(nrow(df1)),
                       simsum = 1)
  to_add <- to_add[!(to_add$.x %in% g1$.x), , drop = FALSE]
  to_add <- to_add[!(to_add$.y %in% g1$.y), , drop = FALSE]
  
  p2 <- bind_rows(g1, to_add)
  attributes(p2) <- attributes(g1)
  
  p3 <- p2 %>%
    as.data.frame() %>%
    mutate(
      maxgrp = max(as.integer(gsub("[^0-9]", "", col3)), na.rm = TRUE),
      col3 = if_else(is.na(col3), paste0("group", maxgrp + cumsum(is.na(col3))), col3)
    ) %>%
    select(-maxgrp)
  
  p4 <- p3 %>%
    mutate(across(1:2, ~ df$barcode[match(., df$id)], .names = "barcodes_{.col}")) %>%
    mutate(across(5:6, .names = "{col}.counts", ~ df$counts[match(., df$barcode)])) %>%
    filter(simsum <= 0 ) %>% 
    mutate(simsum=abs(simsum)) %>% 
    mutate(lv_distane=simsum + 1) %>% 
    relocate(lv_distance, .after = barcodes_.y) %>%
    dplyr::rename(groups = col3) %>%
    dplyr::rename(id = groups, barcode_col1 = barcodes_.x,
                  barcode_col2 = barcodes_.y,
                  counts_col1 = barcodes_.x.counts,
                  counts_col2 = barcodes_.y.counts) %>%
    select(-c(1:3, 7)) %>%
    pivot_longer(cols = -id, names_to = c(".value", "grp"),
                 names_pattern = "(.*)_(col\\d+)") %>%
    ungroup() %>%
    group_by(id) %>%
    distinct(barcode, .keep_all = TRUE) %>%
    mutate(sum_counts = sum(counts)) %>%
    mutate(barcodes = paste(barcode, collapse = ",")) %>%
    select(id, central_barcode = barcode, barcodes, sum_counts, counts) %>%
    ungroup()
  
  return(p4)
}

