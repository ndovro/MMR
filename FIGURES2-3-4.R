# This script generates Figures 2, 3, and 4 for a pan-cancer analysis by processing and visualizing gene expression data.
# It begins by loading required libraries such as ggplot2, dplyr, pheatmap, and others and setting the working directory.
# All CSV files in the directory are loaded, with cancer type information extracted from filenames.
# Significant genes (adjusted p-value < 0.05 and log2 fold change > 1 or < -1) are filtered, and mean log2 fold change values are calculated for each gene per cancer type.

# **Figure 3**:
# The top 2 upregulated and bottom 2 downregulated genes for each cancer type are identified.
# A bar plot is created showing these genes' mean log2 fold changes, categorized by cancer type, with dashed lines indicating thresholds at y = -1 and y = 1.
# The results are saved as "PINAKAS.csv".

# **Figure 2A**:
# A heatmap is generated to visualize the expression patterns of genes across cancer types.
# Hierarchical clustering is applied to rows (cancer types) and columns (genes).
# The heatmap is styled with a red-white-blue color gradient and saved as "heatmap.csv".

# **Figure 2B**:
# Ligand-specific gene expressions are computed using a list of genes associated with each ligand.
# Mean expression levels for these ligands are calculated across samples and visualized in a heatmap using hierarchical clustering.

# **Figure 4**:
# A correlation matrix is created to explore relationships between genes across cancer types using Spearman correlation.
# The correlation matrix is visualized with a color gradient heatmap.

# The script produces heatmaps, bar plots, and correlation visuals to analyze and present key patterns in gene expression data across various cancer types.


# Step 1: Install and load necessary packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(Hmisc) 


setwd("/YOUR/WORKING/DIRECTORY")
# Get all CSV files in the current directory
file_paths <- list.files(pattern = "\\.csv$")

# Function to load CSV files and add a cancer type column extracted from file name
load_csv <- function(file_path) {
  # Extract cancer type from file name
  cancer_type <- str_extract(basename(file_path), "(?<=_)[^.]+")
  df <- read.csv(file_path)
  df$cancer_type <- cancer_type
  return(df)
}

# Load and combine all CSV files
data_list <- lapply(file_paths, load_csv)
combined_data <- bind_rows(data_list)

# Select relevant columns
combined_data <- combined_data %>% dplyr::select(symbol, log2FoldChange, cancer_type,padj)

# Filter out rows with missing log2FoldChange values
combined_data <- combined_data %>% filter(!is.na(log2FoldChange))
combined_data <- combined_data %>% filter(log2FoldChange < -1 | log2FoldChange > 1 )
combined_data <- combined_data %>% filter(padj<0.05)
# Group by cancer type and symbol, calculate mean log2FoldChange
combined_data <- combined_data %>%
  group_by(cancer_type, symbol) %>%
  summarize(mean_log2FoldChange = mean(log2FoldChange))



#=====================================4 genes per cancer type - FIGURE 3 ==================================


# Get top 2 and bottom 2 mean log2FoldChange values for each cancer type
top_bottom_data <- combined_data %>%
  
  group_by(cancer_type) %>%
  arrange(mean_log2FoldChange) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 2 | rank > n() - 2) # Change this line to get top 2 and bottom 2 for each group

# Reorder cancer_type factor levels by mean log2FoldChange within each cancer type
top_bottom_data <- top_bottom_data %>%
  arrange(cancer_type, rank) %>%
  mutate(cancer_type = factor(cancer_type, levels = unique(cancer_type)))

top_bottom_data<-as.data.frame(top_bottom_data)
top_bottom_data$FC<- ifelse(top_bottom_data$mean_log2FoldChange > 0, 2^top_bottom_data$mean_log2FoldChange, -1/2^top_bottom_data$mean_log2FoldChange)
top_bottom_data$FC




# Create a color palette for cancer types
cancer_type_colors <- scales::hue_pal()(length(unique(top_bottom_data$cancer_type)))
# Make the plot
p <- ggplot(top_bottom_data, aes(x = factor(interaction(cancer_type, rank), levels = unique(interaction(top_bottom_data$cancer_type, top_bottom_data$rank))), y = mean_log2FoldChange, fill = cancer_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = cancer_type), position = position_stack(vjust = 0.5), angle = 0, size = 3, fontface = "bold") +  # Add text labels for cancer type inside each bar with bold font
  scale_x_discrete(labels = top_bottom_data$symbol) +  # Change x-axis labels to top_bottom_data$symbol
  scale_fill_manual(values = cancer_type_colors) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "red") +  # Add red dashed lines at y = -1 and y = 1
  theme_minimal() +
  labs(title = "Top 2 UP- and DOWN-Regulated Genes Across Cancer Types",
       x = "Perturbed Gene",
       y = "Mean Log2 Fold Change",
       fill = "Cancer Type") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  theme(legend.position = "none", legend.title = element_blank())+ coord_flip()

p


###############################################HEATMAP AND corrplot - FIGURE 2 A
library(gplots)
library(reshape2)
# Reshape data for heatmap plotting


heatmap_data <- dcast(combined_data, cancer_type ~ symbol, value.var = "mean_log2FoldChange")
# Reshape data for heatmap plotting
  rownames(heatmap_data) <- heatmap_data$cancer_type
heatmap_data <- heatmap_data[, -1]
heatmap_data[is.na(heatmap_data)] <- 0

write.csv(heatmap_data, "summary/heatmap.csv")


library(pheatmap)
heatmap_data <- as.matrix(heatmap_data)



color_palette <- colorRampPalette(c("blue", "white", "red"))(50)

# Define breaks for the color gradient
breaks <- seq(min(heatmap_data, na.rm = TRUE), max(heatmap_data, na.rm = TRUE), length.out = 51)

# Create the heatmap with custom colors, breaks, and legend
pheatmap(heatmap_data, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols =  "euclidean", 
         scale = "none", 
         color = color_palette,
         fontsize_col = 8,
         breaks = breaks,legend = FALSE)

##############################FIGURE 4####################################


# Transpose Data for Correlation
heatmap_data_t <- t(heatmap_data)

# Calculate Spearman Correlation and P-Values
corr_results <- rcorr(heatmap_data_t, type = "spearman")
corr_matrix <- corr_results$r # Correlation matrix
p_matrix <- corr_results$P   # P-value matrix

# Create Significance Levels (Asterisks)
sig_levels <- matrix("", nrow = nrow(p_matrix), ncol = ncol(p_matrix))
sig_levels[p_matrix <= 0.001] <- "***"
sig_levels[p_matrix > 0.001 & p_matrix <= 0.01] <- "**"
sig_levels[p_matrix > 0.01 & p_matrix <= 0.05] <- "*"



long_corr <- as.data.frame(as.table(corr_matrix))
long_corr$Significance <- as.vector(sig_levels)




break_points <- c(-1,  0, 1)
break_colors <- c("blue", "white",  "red")

# Plot Correlation Matrix
ggplot(long_corr, aes(Var1, Var2)) +
  geom_tile(aes(fill = Freq)) +
  scale_fill_gradientn(colors = break_colors, values = scales::rescale(break_points), limits = c(-1, 1)) +
  geom_text(aes(label = Significance), color = "black", size = 4) +
  theme_minimal() +
  labs(x = "", y = "", fill = "Correlation", title = "Spearman correlation of cancer types based on their MMR transcriptomic profile" ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#####################HEATMAP LIGANDS - FIGURE 2B


# Load necessary libraries
library(dplyr)
library(reshape2)
setwd("/YOUR/WORKING/DIRECTORY")

# Read in the files
heatmap_data <- read.csv("heatmap.csv", row.names = 1)
ligands_data <- read.csv("table1.csv")

# Convert ligands_data from wide to long format for easier manipulation
ligands_long <- ligands_data %>%
  pivot_longer(cols = everything(), names_to = "Ligand", values_to = "Gene") %>%
  filter(Gene != "")

# Initialize a data frame to store mean expressions
mean_expressions <- data.frame()

# Loop through each unique ligand
for (ligand in unique(ligands_long$Ligand)) {
  # Get the genes corresponding to the current ligand
  genes <- ligands_long %>%
    filter(Ligand == ligand) %>%
    pull(Gene)
  
  # Subset the heatmap data to include only these genes
  subset_heatmap <- heatmap_data[, colnames(heatmap_data) %in% genes]
  
  # Calculate the mean expression for each row (sample) across these genes
  mean_expression <- rowMeans(subset_heatmap, na.rm = TRUE)
  
  # Store the result in the mean_expressions data frame
  mean_expressions <- rbind(mean_expressions, data.frame(Sample = rownames(heatmap_data), Ligand = ligand, MeanExpression = mean_expression))
}


library(tibble)
library(pheatmap)

# Reshape the mean_expressions data frame to a wide format suitable for pheatmap
mean_expressions_wide <- mean_expressions %>%
  pivot_wider(names_from = Ligand, values_from = MeanExpression) %>%
  column_to_rownames(var = "Sample")




mean_expressions_transposed <- t(mean_expressions_wide)

# Convert the transposed matrix back to a data frame
mean_expressions_transposed_df <- as.data.frame(mean_expressions_transposed)


# Define the breaks and corresponding colors

color_palette <- colorRampPalette(c("blue", "white", "red"))(50)

# Define breaks for the color gradient
breaks <- seq(min(mean_expressions_wide, na.rm = TRUE), max(mean_expressions_wide, na.rm = TRUE), length.out = 51)


# Create the heatmap without a legend
pheatmap(mean_expressions_wide, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         scale = "none", 
         color = color_palette, 
         breaks=breaks,
         legend = F)



# Display the result
print(mean_expressions)








