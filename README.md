```markdown
# Pan-Cancer Gene Expression Analysis and Visualization

## Overview
This repository provides a comprehensive pipeline for analyzing and visualizing pan-cancer RNA-seq data, including CCLE RNA-seq data and TCGA datasets. The analysis includes quantile normalization, differential expression, ligand and receptor mapping, heatmap generation, and bar plot visualizations. It is designed to explore tissue-specific and ligand-specific expression patterns, focusing on microbial metabolite receptors (MMRs) and hallmark genes/pathways.

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
- **Bioconductor Libraries**: `preprocessCore`, `org.Hs.eg.db`, `pheatmap`
- **Visualization Libraries**: `RColorBrewer`, `ggcorrplot`, `gplots`

### Installation
Install Bioconductor libraries if not already installed:
```r
BiocManager::install("preprocessCore")
BiocManager::install("org.Hs.eg.db")
```

Install other required libraries:
```r
install.packages(c("tidyverse", "ggplot2", "reshape2", "readr", "RColorBrewer", "ggcorrplot", "gplots"))
```

---

## Usage

### Preprocess CCLE Data
- Use `ccle_data.R` to normalize CCLE RNA-seq data, map genes, and identify tissue-specific expression patterns.
- **Outputs**: `ccle_tpm_quantiles.csv`, `ccle_tpm_site.csv`.

### Pan-Cancer Analysis
- Use `recount3_expression.R` to analyze TCGA data, identify differentially expressed genes, and annotate results.
- **Outputs**: Differential expression results for genes and pathways.

### Ligand and Hallmark Pathway Analysis
- Use `hallmark_genes.R` and `hallmark_pathways.R` to assess ligand and hallmark pathway expression.
- **Outputs**: Ligand-specific and pathway-specific heatmaps and bar plots.

### Generate Figures
- Use `FIGURES2-3-4.R` to generate visualizations:
  - **Figure 2**: Heatmap of gene expression and ligand pathways.
  - **Figure 3**: Bar plot of up- and down-regulated genes per cancer type.
  - **Figure 4**: Correlation analysis heatmap.

---

## Output Files
- **Normalized Data**: `ccle_tpm_quantiles.csv`, `ccle_tpm_site.csv`
- **Heatmaps**: `heatmap.csv`
- **Bar Plots**: Plots saved in the specified directory for each figure.
- **Differential Expression**: CSV files for full and filtered gene sets.

---

## Directory Structure
```
/data          # Input and output CSV files
/scripts       # R scripts for each step of the analysis
/figures       # Generated plots and visualizations
```

---

## Contributing
Contributions to improve the pipeline, fix bugs, or add new features are welcome. Please create a pull request or open an issue for discussion.

---

## License
This project is licensed under the MIT License. 

---

## Acknowledgments
- **Data Sources**: TCGA and CCLE for providing publicly available RNA-seq data.
- **Libraries**: Developers of the `recount3` and `preprocessCore` libraries for enabling robust data handling.
```
