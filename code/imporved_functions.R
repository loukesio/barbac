library(dplyr)
library(stringr)
library(gt)

process_data <- function(file_path) {
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
files[1]
table <- function(file, barcode_length) {
  # Check if the 'barcode_length' parameter is a numeric vector of length 2
  if (!is.numeric(barcode_length) || length(barcode_length) != 2) {
    stop("Error: 'barcode_length' should be a numeric vector of length 2.")
  }
  
  # Check if 'file' is a data frame or tibble
  if (!is.data.frame(file) && !is_tibble(file)) {
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
    gt() %>%
    cols_label(n = "No of barcodes") %>%
    tab_source_note(
      source_note = "Barbac: A versatile tool to quantify barcodes."
    ) %>%
    gt_theme_espn()
}
plot_rawbar <- function(raw_barcodes) {
  # Check if 'raw_barcodes' is a data frame or tibble
  if (!is.data.frame(raw_barcodes) && !is_tibble(raw_barcodes)) {
    stop("Error: 'raw_barcodes' should be a data frame or tibble.")
  }
  
  # Check if 'length_barcode' column exists in 'raw_barcodes'
  if (!"length_barcode" %in% colnames(raw_barcodes)) {
    stop("Error: 'raw_barcodes' does not contain the 'length_barcode' column.")
  }
  
  ggplot(raw_barcodes, aes(x = length_barcode)) +
    geom_histogram(fill = "#FF5349", color = NA, alpha = 0.8) +
    labs(
      x = "Barcode length (bp)",
      y = "Number of barcodes"
    ) +
    theme_bw()
}


