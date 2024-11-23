# This script analyzes gene expression data focusing on Microbial Metabolite Receptors (MMR) and Cancer Hallmark genes.
# It begins by loading libraries such as tidyr, dplyr, ggplot2, and others.
# The working directory is set, and a list of target genes (MMR and Hallmark) is defined.
# A directory named "filtered" is created to store filtered data if it doesn't already exist.
# The script processes all CSV files starting with "FULL_" to filter rows with significant genes (adjusted p-value < 0.05).
# Filtered data is saved in the "filtered" directory, ensuring each file has unique entries and reordered columns.
# Next, the script consolidates data by combining all filtered files into a single dataframe, using gene symbols as row names.
# Rows corresponding to MMR and Hallmark genes are extracted for correlation analysis.
# Pairwise Spearman correlations are calculated between MMR and Hallmark genes.
# The correlation matrix is saved as "correlation_full.csv".
# A summary of the strongest correlations for each MMR gene is saved as "highest_absolute_correlationa.csv".
# Correlations above 0.8 or below -0.8 are identified and visualized using a heatmap.
# The heatmap depicts correlations between MMR and Hallmark genes across 24 cancer types.

library(tidyr)
library(dplyr)
library(purrr)
library(readr)
library(ggplot2)
library(GGally)
library(Hmisc)
setwd("/YOUR/WORKING/DIRECTORY")

MMR <-c("FFAR1", "FFAR2", "FFAR3", "FFAR4", "HCAR1", "HCAR2", "HCAR3", "GPR84", "PPARA", "PPARD", "PPARG", "PPARGC1A", "PPARGC1B", "LPAR1", "LPAR2", "LPAR3", "LPAR4", "LPAR5", "LPAR6", "AHR", "AHRR", "GPR35", "HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F", "HTR2A", "HTR2B", "HTR2C", "HTR3A", "HTR3B", "HTR3C", "HTR4", "HTR5A", "HTR6", "HTR7", "GPBAR1", "CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5", "NR1H4", "VDR", "NR1I2", "NR1I3", "GABRA1", "GABRB2", "GABRG2", "AR", "ESR1", "ESR2", "S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5", "RARA", "RARB", "RARG", "RXRA", "RXRB", "RXRG", "TRPV1", "TRPV2", "TRPV3", "TRPV4", "TRPA1", "TRPM8", "CNR1", "CNR2", "GPR55", "GPR119", "HRH1", "HRH2", "HRH3", "HRH4", "DRD1", "DRD2", "DRD3", "DRD4", "ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C", "ADRB1", "ADRB2", "ADRB3", "ADORA1", "ADORA2A", "ADORA3", "P2RY1", "P2RY2", "P2RY4", "P2RY6", "P2RY8", "P2RY10", "P2RY11", "P2RY12", "P2RY13", "P2RY14", "P2RX1", "P2RX4", "P2RX7")


HALLMARK<- read.csv("big_hallmark.csv", header= F)
HALLMARK<- HALLMARK$V1

# Create the filtered subdirectory if it doesn't exist
if (!dir.exists("filtered")) {
  dir.create("filtered")
}

# List all CSV files that start with "FULL_"
files <- list.files(pattern = "^FULL_.*\\.csv$")

# Loop through each file
for (file in files) {
  # Read the CSV file
  data <- read.csv(file)
  
  # Filter the rows where the symbol column matches any value in MMR or HALLMARK
  filtered_data <- data %>%
    filter((symbol %in% MMR | symbol %in% HALLMARK) & padj < 0.05)%>%
    select(log2FoldChange, symbol) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    distinct(symbol, .keep_all = TRUE)
  
  # Reorder columns to set symbol as the first column
  filtered_data <- filtered_data %>%
    select(symbol, log2FoldChange)
  
  # Create the new file name for the filtered data
  new_file <- file.path("filtered", file)
  
  # Save the filtered data to the new file
  write.csv(filtered_data, new_file, row.names = FALSE)
}


###############################SECOND PART

setwd("/YOUR/WORKING/DIRECTORY/filtered")


# List all CSV files in the directory
files <- list.files(pattern = "FULL_.*\\.csv")

# Function to read each CSV and extract relevant columns
read_and_extract <- function(file) {
  df <- read.csv(file)
  df <- df %>% select(symbol, log2FoldChange)
  colnames(df)[2] <- sub("FULL_", "", tools::file_path_sans_ext(file))
  return(df)
}

# Read all files and extract data
data_list <- map(files, read_and_extract)
data_list <- data_list[lapply(data_list, nrow) > 0]
# Merge all dataframes by the symbol column
combined_df <- reduce(data_list, full_join, by = "symbol")
combined_df[is.na(combined_df)] <- 0

# Set symbol as row names and remove the symbol column
rownames(combined_df) <- combined_df$symbol
combined_df <- combined_df %>% select(-symbol)


# Extract rows for MMR and HALLMARK symbols from combined_df
MMR_df <- combined_df[rownames(combined_df) %in% MMR, ]
HALLMARK_df <- combined_df[rownames(combined_df) %in% HALLMARK, ]


# Initialize matrices to store results
correlation_matrix <- matrix(NA, nrow = nrow(MMR_df), ncol = nrow(HALLMARK_df))
p_value_matrix <- matrix(NA, nrow = nrow(MMR_df), ncol = nrow(HALLMARK_df))

# Calculate pairwise Spearman correlations and p-values
for (i in 1:nrow(MMR_df)) {
  for (j in 1:nrow(HALLMARK_df)) {
    # Perform rcorr for each pair of rows
    result <- rcorr(as.numeric(MMR_df[i, ]), as.numeric(HALLMARK_df[j, ]), type = "spearman")
    
    # Extract the correlation and p-value
    correlation_matrix[i, j] <- result$r[1, 2]
    p_value_matrix[i, j] <- result$P[1, 2]
  }
}

# Set row and column names for the matrices
rownames(correlation_matrix) <- rownames(MMR_df)
colnames(correlation_matrix) <- rownames(HALLMARK_df)

rownames(p_value_matrix) <- rownames(MMR_df)
colnames(p_value_matrix) <- rownames(HALLMARK_df)






# Convert correlation matrix to data frame for ggplot2
correlation_df <- as.data.frame(as.table(correlation_matrix))
correlation_df$P_value <- as.vector(p_value_matrix)

write.csv (correlation_df,"correlation_full.csv")


result <- correlation_df %>%
  group_by(Var1) %>%
  slice(which.max(abs(Freq)))

write.csv (result,"highest_absolute_correlationa.csv")

# Filter for correlations > 0.8 and < -0.8
filtered_df <- correlation_df %>% 
  filter(Freq > 0.8 | Freq < -0.8)


# Visualize using ggplot2
ggplot(filtered_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() +
  labs(x = "Microbial Metabolite Receptors", y = "Cancer Hallmark Genes", title = "Correlation Heatmap across the 24 cancer types")





