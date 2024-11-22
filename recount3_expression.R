# This script processes RNA-seq data for multiple cancer projects using the recount3 and DESeq2 libraries.
# It begins by setting up the environment, loading necessary libraries, and configuring parallel processing for efficiency.
# A list of MMRs and cancer projects is defined.
# The script retrieves project metadata and extracts count data along with sample information such as tissue types and batch details.
# It filters out projects that lack "Solid Tissue Normal" samples or sufficient group diversity for comparison.
# For each valid project, it creates a DESeq2 dataset, incorporating batch effect correction if needed.
# Differential expression analysis is performed to compare tumor and normal samples.
# The results are annotated with gene symbols, filtered for genes of interest, and exported to CSV files.
# These CSV files include full results and targeted results for the specified genes.
# Finally, a summary of sample types across projects is logged to a text file for review.

library(recount3)
library(DESeq2)
library(org.Hs.eg.db)
library(BiocParallel)
library(dplyr)
register(MulticoreParam(10))

setwd("/YOUR/WORKING/DIRECTORY")

genes_of_interest <- c("FFAR2", "FFAR4", "GPR84", "HCAR2", "HCAR3", "PPARG", "NR1I2", "NR1H4", "S1PR1", "S1PR3", "SPR4", "S1PR4", "CNR1", "TRPV1", "TRPV8")


options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")


TD <- c("BRCA", "KIRC", "LUAD", "UCEC", "THCA", "PRAD", "LUSC", "HNSC", "COAD", "LGG", "SKCM", "STAD", "BLCA", "OV", "LIHC", "KIRP", "CESC", "SARC", "ESCA", "PCPG", "PAAD", "LAML", "READ", "GBM", "TGCT", "THYM", "KICH", "MESO", "UVM", "ACC", "UCS", "DLBC", "CHOL")
options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")
available_projects <- available_projects()

for (project in TD) {
  proj_info <- available_projects[available_projects$project_type == "data_sources" & available_projects$project == project, ]
   rse_gene <- create_rse(proj_info)
  
  
  counts <- transform_counts(rse_gene)
  counts <- as.data.frame(counts)
  
  types <- rse_gene$tcga.cgc_sample_sample_type
  ids <- rse_gene$external_id
  batch <- rse_gene$tcga.cgc_case_batch_number
  
  metadata <- data.frame(ids, types, batch)
  
  summary(metadata$batch)
  
  metadata <- metadata[complete.cases(metadata), ] 
  
  counts <- counts[, (colnames(counts) %in% metadata$ids)]
  
  
  if (!("Solid Tissue Normal" %in% levels(as.factor(metadata$types))) || length(levels(as.factor(metadata$types))) < 2) {
    # Skip this iteration and move to the next loop
    print(paste0("SKIPPING",project))
    next
  }
  
  if (length(levels(as.factor(metadata$batch))) >= 2) {
    metadata$batch<-as.factor(metadata$batch)
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = metadata, 
                                design = ~types + batch) }
  else {
    dds <- DESeqDataSetFromMatrix(countData = counts, 
                                  colData = metadata, 
                                  design = ~types )
   
  }
  
  dds$types <- relevel(dds$types, ref = "Solid Tissue Normal")
  
  register(MulticoreParam(10))
  dds <- DESeq(dds, parallel = TRUE)
  
  resultsNames(dds) #LOCATION
  
  res <- results(dds, name = "types_Primary.Tumor_vs_Solid.Tissue.Normal", parallel = TRUE)
  resOrdered <- res[order(res$padj), ]
  
  resOrdered$symbol <- mapIds(org.Hs.eg.db, 
                              keys = gsub("\\..*", "", rownames(resOrdered)), 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL",
                              multiVals = "first")
  write.csv(as.data.frame(resOrdered), file = paste0("FULL_", project, ".csv"))
  resFiltered <- subset(resOrdered, symbol %in% genes_of_interest)
  resFiltered$FC<- ifelse(resFiltered$log2FoldChange > 0, 2^resFiltered$log2FoldChange, -1/2^resFiltered$log2FoldChange)
  
  write.csv(as.data.frame(resFiltered), file = paste0("metabolites2_", project, ".csv"))
}


sink(nik.txt)
for (project in TD) {
  proj_info <- available_projects[available_projects$project_type == "data_sources" & available_projects$project == project, ]
  rse_gene <- create_rse(proj_info)
  
  print(project)
  print(table(as.factor(rse_gene$tcga.cgc_sample_sample_type)))

  }
sink()
