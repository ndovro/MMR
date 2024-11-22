```markdown
# Pan-cancer insights: A Study of Microbial Metabolite Receptors in Malignancy Dynamics.


---

## Features

1. **Preprocessing**:
   - Quantile normalization of CCLE RNA-seq data.
   - Mapping gene identifiers (Ensembl to gene symbols).
   - Filtering genes of interest (e.g., MMRs and hallmark genes).

2. **Differential Expression**:
   - Analysis across cancer types from TCGA datasets using recount3.
   - Identification of significantly up- and down-regulated genes.

3. **Ligand Analysis**:
   - Mapping ligands to genes and assessing tissue-specific expression patterns.

4. **Visualizations**:
   - Heatmaps for expression and correlation analysis.
   - Bar plots for tissue-specific and ligand-specific expression.

5. **Outputs**:
   - Normalized expression data files.
   - Figures depicting expression trends across tissue types and pathways.

---

## Setup

### Required Libraries
- **Base Libraries**: `tidyverse`, `dplyr`, `ggplot2`, `reshape2`, `readr`
- **Bioconductor Libraries**: `preprocessCore`, `org.Hs.eg.db`, `pheatmap`, `recount3`
- **Visualization Libraries**: `RColorBrewer`, `ggcorrplot`, `gplots`


Install other required libraries:
```r
install.packages(c("tidyverse", "ggplot2", "reshape2", "readr", "RColorBrewer", "ggcorrplot", "gplots"))
```

---
## Overview
Research increasingly shows that bacteria, the largest part of the gut microbiome, may play a significant role in cancer. These microorganisms produce metabolites that can travel through the body, and interact with host cells via microbial metabolite receptors, potentially affecting cancer development. This study investigates the involvement of these receptors in human cells across twenty-three types of cancer. By analyzing data from both cancer cell lines and human tumor samples, we examined how these interactions may impact key cancer-related processes, such as immune response, tumor growth, and spread. Notably, we identified several receptors that are consistently altered in cancer, which might serve as helpful biomarkers for diagnosis or treatment. This research highlights the potential of targeting the gut microbiome in cancer therapy and provides valuable insights for developing new cancer treatments based on microbiome interactions.

## Usage

### Preprocess CCLE Data
- Use `ccle_data.R` to normalize CCLE RNA-seq data, map genes, and identify tissue-specific expression patterns.
- **Outputs**: `ccle_tpm_quantiles.csv`, `ccle_tpm_site.csv`.
- This File is required for Figure 1 and Supplementary Figure 1
- Figure 1. Microbial Metabolite Receptor (MMR) expression in the Cancer Cell Line Encyclopedia (CCLE) dataset. A) MMR with the highest expression per tissue-specific cell line B) MMRs summarized by their ligand showing which tissue-specific cell line those are highly expressed in.

### Pan-Cancer Analysis
- Use `recount3_expression.R` to analyze TCGA data, identify differentially expressed genes, and annotate results.
- **Outputs**: Differential expression results for genes.

### Ligand and Hallmark Pathway Analysis
- Use `hallmark_genes.R` and `hallmark_pathways.R` to assess ligand and hallmark pathway correlations to MMRs.
- **Outputs**: Ligand-specific and pathway-specific heatmaps and bar plots. These create Figure 5 components.
- Figure 5. A) Pairwise Spearman correlation of Cancer Hallmark Genes (CHGs) and Microbial Metabolite Receptors (MMRs) based on their expression B) Pairwise Spearman correlation of MMRs when the CHGs are aggregated into Cancer Hallmark Pathways (CHPs) C) Summary of the strongest correlations between MMRs and CHPs

### Generate Figures
- Use `FIGURES2-3-4.R` to generate visualizations:
  - Figure 2. Differential Expression of Microbial Metabolite Receptors (MMRs) in a pancancer setting. A) Each MMR’s dysregulation versus control samples per studied cancer type  B) MMR expression dysregulation summarized by ligand type per studied cancer type. Red color signifies Upregulation and blue Downregulation.
  - Figure 3. Top dysregulated Microbial Metabolite Receptors (MMRs) per studied cancer type. X-axis represents log2FoldChange values while the dashed red lines signify log2FoldChange of ±1.
  - Figure 4. Spearman Correlation Heatmap showcasing the relationship between cancer types based on MMR expression profiles. Red indicates a positive correlation while blue highlights inverse correlations.

---

### Important!!
- Change the directories where necessary
- All necessary input files are provided.

---

## License
This project is licensed under the MIT License. 

---

## Acknowledgments
- **Data Sources**: TCGA and CCLE for providing publicly available RNA-seq data.

