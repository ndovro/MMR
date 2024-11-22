# This script analyzes correlations between Microbial Metabolite Receptors (MMR) and Cancer Hallmark Pathways.
# It starts by loading libraries such as dplyr, tidyr, and ggplot2 and setting the working directory.
# The script reads two input CSV files: "correlation_full.csv" for correlation data and "categories.csv" for gene categories.
# The categories are transformed into a long format and merged with the correlation data.
# Rows with missing correlation values are removed, and the median correlation frequency is calculated for each gene pair.
# The top 25% and bottom 25% correlations are filtered, and the filtered data is saved as "network.csv".
# A heatmap visualizes correlations using ggplot2, showing relationships between MMR and Cancer Hallmark categories.
# To manage visualization size, the data is split into two parts and displayed as separate heatmaps.
# The script identifies the top 3 and bottom 3 correlations for each category, saving the results in "forpaperpathways.csv".




library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)


setwd("/YOUR/WORKING/DIRECTORY")


df_correlation <- read.csv('correlation_full.csv')
df_categories <- read.csv('categories.csv')

df_long <- df_categories %>%
  pivot_longer(cols = everything(), 
               names_to = "Category", 
               values_to = "Var2")

df_correlation$Var2 <- trimws(df_correlation$Var2)
df_long$Var2 <- trimws(df_long$Var2)


df_long<-as.data.frame(df_long)

df_merged <- merge(df_correlation, df_long, by = "Var2")

df_merged <- df_merged %>% select(-X)
df_merged <- df_merged %>% select(-Var2)




clean_df <- df_merged %>%
  filter(!is.na(Freq))

result_df <- clean_df %>%
  group_by(Var1, Category) %>%
  summarize(Freq = median(Freq, na.rm = TRUE)) %>%
  ungroup()


filtered_df<-result_df

# Calculate the quantiles for the top 25% and bottom 25%
quantiles <- quantile(filtered_df$Freq, probs = c(0.25, 0.75))

# Filter the dataframe to keep only the top 25% and bottom 25%
filtered_df <- filtered_df %>%
  filter(Freq <= quantiles[1] | Freq >= quantiles[2])

write.csv(filtered_df, "network.csv", row.names = FALSE)

# Visualize using ggplot2
ni<-ggplot(filtered_df, aes(Category, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, space = "Lab",
                       name="Mean Correlations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1),
        axis.text.y = element_text(size = 6),legend.position = "none") +  # Reduce y-axis text size
  coord_fixed(ratio = 0.5) +  # Adjust aspect ratio
  labs(x = "Cancer Hallmarks", y = "Microbial Metabolite Receptors", title = "Correlation Heatmap across the 24 cancer types")



library(gridExtra)  # For arranging multiple plots

# Divide the data into two parts
half_index <- ceiling(nrow(filtered_df) / 2)
filtered_df_part1 <- filtered_df[1:half_index, ]
filtered_df_part2 <- filtered_df[(half_index + 1):nrow(filtered_df), ]

# First plot
plot1 <- ggplot(filtered_df_part1, aes(Category, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, space = "Lab",
                       name="Mean Correlations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 6, hjust = 1),
        axis.text.y = element_text(size = 6),legend.position = "none") +
  coord_fixed(ratio = 0.5) +
  labs(x = "Cancer Hallmarks", y = "Microbial Metabolite Receptors", title = "Pan-cancer Correlation between MMRs and Cancer Hallmarks")

# Second plot
plot2 <- ggplot(filtered_df_part2, aes(Category, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, space = "Lab",
                       name="Mean Correlations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 6, hjust = 1),
        axis.text.y = element_text(size = 6),legend.position = "none") +
  coord_fixed(ratio = 0.5) +
  labs(x = "", y = "Microbial Metabolite Receptors", title = "")

# Arrange the plots side by side
grid.arrange(plot1, plot2, ncol = 2, widths = unit(c(0.2, 0.2), "npc"), padding = unit(0.5, "line"))




# Function to get top and bottom 5 based on frequency
get_top_bottom_5 <- function(df) {
  top5 <- df %>% arrange(desc(Freq)) %>% head(3)
  bottom5 <- df %>% arrange(Freq) %>% head(3)
  bind_rows(top5, bottom5)
}

# Apply the function to each category
result <- filtered_df %>%
  group_by(Category) %>%
  do(get_top_bottom_5(.))%>%
  arrange(Category, desc(Freq))

write.csv(result, "forpaperpathways.csv", row.names = FALSE)




