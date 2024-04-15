

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

