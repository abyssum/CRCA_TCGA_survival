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

# Load in the data
clinical <- read.table("/home/fotakis/myScratch/PHD_Project/TCGA_survival/clinicalCRC.txt")

expression.clean <- read.csv("/home/fotakis/myScratch/PHD_Project/TCGA_survival/expression_COADREAD_clean.tsv",
                             header = T,
                             sep = ",")

mappings <- read.csv("/home/fotakis/myScratch/PHD_Project/TCGA_survival/GDCidMapping_20190712_full.tsv",
                     header = T,
                     sep = "\t")

crca.de.genes <- read.csv2("/home/fotakis/myScratch/CRCA_survival/immune_infiltration-DESeq2_result.tsv", 
                           header = T, 
                           sep = "\t")