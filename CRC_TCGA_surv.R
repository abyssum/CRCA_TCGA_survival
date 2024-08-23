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

# Edit the colnames to match with the IDs (will be used for mapping)
colnames(expression.clean) <- gsub("\\.", "-", colnames(expression.clean))

# Clean up the expression DF and use te gene symbols as index
gene.names <- expression.clean$GeneName
expression.clean <- expression.clean[,-c(1)]
# Convert to matrix - required by GSVA
expression.clean <- as.matrix(expression.clean)
rownames(expression.clean) <- gene.names

# Define the DE comparisons
comparisons <- c("B", "desert", "M", "T")

# Create a list to store the results
de.genes.list <- lapply(comparisons, function(comp) {
  genes <- filter(crca.de.genes, comparison == comp)
  genes$log2FoldChange <- as.numeric(genes$log2FoldChange)
  return(genes)
})

# Assign names to the list elements for easy access
names(de.genes.list) <- comparisons

# Access specific comparison from the initial DE gene DF
b.de.genes <- de.genes.list[["B"]]
desert.de.genes <- de.genes.list[["desert"]]
m.de.genes <- de.genes.list[["M"]]
t.de.genes <- de.genes.list[["T"]]

# Define the filtering params and filter the DE gene DFs
n_rows <- 15
lFC_trheshold <- 0.5
column <- "log2FoldChange"
# B group
b.de.genes.top <- filter_and_trim(b.de.genes, column, lFC_trheshold, n_rows)
# Desert group
desert.de.genes.top <- filter_and_trim(desert.de.genes, column, lFC_trheshold, n_rows)
# M group
m.de.genes.top <- filter_and_trim(m.de.genes, column, lFC_trheshold, n_rows)
# T group
t.de.genes.top <- filter_and_trim(t.de.genes, column, lFC_trheshold, n_rows)

# Define the named list of top DE genes per group
goi <- list(c(b.de.genes.top$symbol),
            c(desert.de.genes.top$symbol),
            c(m.de.genes.top$symbol),
            c(t.de.genes.top$symbol))
names(goi) <- c("B", "desert", "M", "T")

# Run the GSVA and transpose the results for easier processing
gsva.res <- gsva(expression.clean, goi, verbose = FALSE)
test_df <- as.data.frame(t(gsva.res), stringsAsFactors = FALSE)
test_df$case_ID <- rownames(test_df)

# Create a mapping dataframe with unique rows, only for the common case_IDs
df.map <- unique(data.frame(
  bcr_patient_barcode = mappings$cases.0.submitter_id,
  case_ID = mappings$cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id,
  stringsAsFactors = FALSE
))
df.map <- df.map[df.map$case_ID %in% test_df$case_ID, ]

# Merge the dataframes: GSVA results, mappings, and clinical data
test.final_df <- merge(clinical, merge(test_df, df.map, by = "case_ID"), by = "bcr_patient_barcode")

# Convert OS time to years
test.final_df$OS.time.years <- test.final_df$OS.time / 365

# Define the cutpoint - to be used for the ranking: high/low
surv_cut <- surv_cutpoint(
  test.final_df,
  time = "OS.time.years",
  event = "OS",
  variables = c("B", "desert", "M", "T")
)
summary(surv_cut)


