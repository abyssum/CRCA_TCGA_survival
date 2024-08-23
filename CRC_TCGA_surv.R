# Load libraries
library("survival")
library("survminer")
library(dplyr)
library(GSVA)

# Function to filter the DFs and extract top DE genes
filter_and_trim <- function(data, column_name, lFC_threshold, n_rows) {
  # Filter rows based on the lFC_threshold
  filtered_data <- data[data[[column_name]] >= lFC_threshold | data[[column_name]] <= -lFC_threshold, ]
  
  # Keep only the top n rows
  trimmed_data <- head(filtered_data, n_rows)
  
  return(trimmed_data)
}

