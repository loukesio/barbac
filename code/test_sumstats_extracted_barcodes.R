#____________________________ 
# ┬  ┬┌┐ ┬─┐┌─┐┬─┐┬┌─┐┌─┐
# │  │├┴┐├┬┘├─┤├┬┘│├┤ └─┐
# ┴─┘┴└─┘┴└─┴ ┴┴└─┴└─┘└─┘
#__________________________

library(here)
library(dplyr)
library(stringr)
library(gt)
library(gridExtra)


#______________________________________________________________________
# ╔═╗╦ ╦╔╗╔╔═╗╔╦╗╦╔═╗╔╗╔  ╔═╗╔╗╔╔═╗       ┬─┐┌─┐┌─┐┌┬┐  ┌─┐┬┬  ┌─┐┌─┐
# ╠╣ ║ ║║║║║   ║ ║║ ║║║║  ║ ║║║║║╣   ───  ├┬┘├┤ ├─┤ ││  ├┤ ││  ├┤ └─┐
# ╚  ╚═╝╝╚╝╚═╝ ╩ ╩╚═╝╝╚╝  ╚═╝╝╚╝╚═╝       ┴└─└─┘┴ ┴─┴┘  └  ┴┴─┘└─┘└─┘
#______________________________________________________________________

files=list.files(here("data","barbac_testdata"), pattern = "*barcodes.txt$", full.names = TRUE)   # select all the files that have isolate barcodes 
files         # all files

read_outbam <- function(file_path) {
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

test <- read_outbam(files[1])

###########################
# table 1 barcode - bins
##########################

file <- test  # Assuming 'data' is your input dataframe
barcode_length <- c(22, 28)  # Specify the bin ranges here

sumstats(file,barcode_length, fill_color = "blue")

