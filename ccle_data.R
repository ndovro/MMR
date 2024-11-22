# This script analyzes CCLE (Cancer Cell Line Encyclopedia) RNA-seq data to assess the expression of Microbial Metabolite Receptors (MMR) across tissue types and ligand classes.

# **Data Preparation**:
# - RNA-seq TPM data is loaded and preprocessed, with gene symbols mapped from Ensembl IDs using `org.Hs.eg.db`.
# - The dataset is filtered for specific genes of interest (MMR-related genes).
# - Quantile normalization is applied to standardize expression values, and the normalized data is saved as "ccle_tpm_quantiles.csv".
# - Tissue site information is extracted from sample identifiers and formatted for analysis.

# **Figure 1A**:
# - The tissue type with the highest expression of each MMR gene is identified.
# - A bar plot displays the log2-transformed mean expression values of the most expressed MMR gene per tissue type, highlighting gene and tissue type associations.

# **Supplementary Figure 1A**:
# - The MMR gene with the highest mean expression for each tissue type is identified.
# - A bar plot visualizes log2-transformed mean TPM values for these genes across tissues.

# **Figure 1B**:
# - Ligand categories are mapped to genes, and the site with the highest expression of each ligand category is identified.
# - A bar plot shows log2-transformed mean TPM values for each ligand category and the corresponding tissue type.



library(tidyverse)
library(preprocessCore)
library(org.Hs.eg.db)

setwd("/YOUR/WORKING/DIR")
ccle<-read.table("CCLE_RNAseq_rsem_genes_tpm_20180929.txt")


colnames(ccle)<- ccle[1,]
ccle<-ccle[-1,]


rownames(ccle)<- ccle[,1]
ccle<-ccle[,-1]



ccle$symbol <- mapIds(org.Hs.eg.db, 
                            keys=gsub("\\..*","", rownames(ccle)), 
                            column="SYMBOL", 
                            keytype="ENSEMBL",
                            multiVals="first")

genes_of_interest <-c("FFAR1", "FFAR2", "FFAR3", "FFAR4", "HCAR1", "HCAR2", "HCAR3", "GPR84", "PPARA", "PPARD", "PPARG", "PPARGC1A", "PPARGC1B", "LPAR1", "LPAR2", "LPAR3", "LPAR4", "LPAR5", "LPAR6", "AHR", "AHRR", "GPR35", "HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F", "HTR2A", "HTR2B", "HTR2C", "HTR3A", "HTR3B", "HTR3C", "HTR4", "HTR5A", "HTR6", "HTR7", "GPBAR1", "CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5", "NR1H4", "VDR", "NR1I2", "NR1I3", "GABRA1", "GABRB2", "GABRG2", "AR", "ESR1", "ESR2", "S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5", "RARA", "RARB", "RARG", "RXRA", "RXRB", "RXRG", "TRPV1", "TRPV2", "TRPV3", "TRPV4", "TRPA1", "TRPM8", "CNR1", "CNR2", "GPR55", "GPR119", "HRH1", "HRH2", "HRH3", "HRH4", "DRD1", "DRD2", "DRD3", "DRD4", "ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C", "ADRB1", "ADRB2", "ADRB3", "ADORA1", "ADORA2A", "ADORA3", "P2RY1", "P2RY2", "P2RY4", "P2RY6", "P2RY8", "P2RY10", "P2RY11", "P2RY12", "P2RY13", "P2RY14", "P2RX1", "P2RX4", "P2RX7")


cclefilt <- subset(ccle, symbol %in% genes_of_interest)

rownames(cclefilt)<-cclefilt$symbol

cclefilt <- cclefilt[, colnames(cclefilt) != "symbol"]

nik<-cclefilt

# Remove row names
rownames(nik) <- NULL

# Remove column names
colnames(nik) <- NULL

# Convert the data frame to numeric
nik <- as.data.frame(lapply(nik, as.numeric))

nik <- as.matrix(nik)

###PERFORM QUANTILES NORMALIZATION
normalized_data <- normalize.quantiles(nik)

nik<-as.data.frame (normalized_data)
rownames(nik)<-rownames(cclefilt)
colnames(nik)<-colnames(cclefilt)

write.csv(nik, file="ccle_tpm_quantiles.csv")



sites<-read.csv("sites.csv", header = T)
  
tpm<-read.csv("ccle_tpm_quantiles.csv")

rownames(tpm)<-tpm$Gene
tpm<-tpm[,-1]

tpm<-as.data.frame(t(tpm))       
         

tpm$site <- as.factor(sub("^[^_]+_", "", rownames(tpm)))


to_sentence_case <- function(text) {
  text <- str_to_lower(text)
  text <- gsub("_", " ", text)
  text <- str_to_sentence(text)
  return(text)
}

# Apply function to the site column
tpm$site <- as.factor(sapply(tpm$site, to_sentence_case))

write.csv(tpm, file="ccle_tpm_site.csv")


library(reshape2)
molten_data <- melt(tpm, id.vars = "site")
molten_data<-na.omit(molten_data)

means <- molten_data %>%
  group_by(site, variable) %>%
  summarise(mean = mean(value))




####CCLE by site - FIGURE 1 A

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
# Get the variable with the highest mean value for each site
max_means <- means %>%
  group_by(site) %>%
  filter(mean == max(mean)) %>%
  ungroup()

# Convert mean values to log2 scale
max_means <- max_means %>%
  mutate(mean_log2 = log2(mean + 1)) # Adding 1 to avoid log2(0)


library(stringr)

# Create the plot
ggplot(max_means, aes(x = mean_log2, y = reorder(site, mean_log2), fill = mean_log2)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = variable), hjust = -0.1, size = 3) +
  labs(title = "Highest Expressed MMR for each Tissue Type", 
       x = "Log2 Mean normalized TPM", 
       y = "Site") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y= element_text(size = 10) )+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))

#########CCLE BY MMR - SUPPL FIGURE 1 A

library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# Get the site with the highest mean value for each variable
max_means <- means %>%
  group_by(variable) %>%
  filter(mean == max(mean)) %>%
  ungroup()

# Convert mean values to log2 scale
max_means <- max_means %>%
  mutate(mean_log2 = log2(mean + 1)) # Adding 1 to avoid log2(0)


# Create the plot
ggplot(max_means, aes(x = mean_log2, y = reorder(variable, mean_log2), fill = mean_log2)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = site), hjust = -0.1, size = 2) +
  labs(title = "Highest Mean Normalized TPM by Site for Each Gene (log2)", 
       x = "Log2 Mean normalized TPM", 
       y = "Gene") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))




######################LIGAND TYPES - FIGURE 1B
# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(readr)


# Load the CSV file containing gene categories
gene_categories <- read_csv("ligand_class.csv")
# Reshape the data from wide to long format
gene_mapping <- gene_categories %>%
  pivot_longer(cols = everything(), names_to = "category", values_to = "gene") %>%
  filter(!is.na(gene))

# Assuming 'means' dataframe has columns: 'site', 'variable' (gene names), and 'mean'
# Join the means dataframe with the gene mapping to get the categories
means <- means %>%
  left_join(gene_mapping, by = c("variable" = "gene"))

# Group by category and site to find the highest mean for each category
max_means <- means %>%
  group_by(category, site) %>%
  summarise(max_mean = max(mean, na.rm = TRUE)) %>%
  ungroup()

# Find the site with the highest mean for each category
max_means <- max_means %>%
  group_by(category) %>%
  filter(max_mean == max(max_mean)) %>%
  ungroup()

# Convert mean values to log2 scale
max_means <- max_means %>%
  mutate(mean_log2 = log2(max_mean + 1)) # Adding 1 to avoid log2(0)

ggplot(max_means, aes(x = mean_log2, y = reorder(category, mean_log2), fill = mean_log2)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = site), position = position_stack(vjust = 0.5), size = 3, color = "white") +
  labs(title = "Tissue Type containing the Highest Expression of ligands", 
       x = "Log2 Mean normalized TPM", 
       y = "Ligands") +
  scale_fill_gradient(low = "black", high = "darkblue") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y= element_text(size = 12))