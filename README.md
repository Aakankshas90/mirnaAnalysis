# miRNA_seq meta-analysis

This repository contains the **miRNA Analysis Pipeline**, developed to analyze publicly available NGS datasets (miRNA-seq) for Gastric Cancer, collected from public databases, like TCGA. The workflow is implemented using automated scripts and R packages for preprocessing, differential expression analysis, meta-analysis, and functional analysis.

## Key Features

1. **Preprocessing and Alignment**
   - The script `mirna.sh` performs the following tasks:
     - Quality check of raw reads using **FastQC**.
     - Trimming of adapters and low-quality bases using **fastp** or similar tools.
     - Alignment of reads to the genome and read counts using **miRDeep2**.

2. **Differential Expression Analysis**
   - Differential expression analysis is conducted on case vs. control datasets from six bioprojects.
   - The results include:
     - Log Fold Change (LFC) values.
     - Standard error estimates for each miRNA.

3. **Meta-Analysis**
   - Meta-analysis of miRNAs is conducted using the **metafor** package in R.
   - Filtering criteria:
     - miRNAs must be detected in multiple studies.
     - Log Fold Change (LFC) and standard error ratio (SE) are used to refine results.
   - Outputs include forest plots for selected miRNAs.

4. **Functional Analysis**
   - Target genes of significant miRNAs are identified using tools such as **miRDB**.
   - Gene ontology (GO) term enrichment is performed using **clusterProfiler**.
   - GO terms are visualized with **ggplot2**.

## Tools and Packages Used

- **Shell Scripting:** Automated preprocessing and alignment steps.
- **R Packages:**
  - `DESeq2` for differential expression analysis.
  - `metafor` for meta-analysis.
  - `clusterProfiler` for functional enrichment analysis.
  - `ggplot2` for data visualization.

## Workflow Overview

1. **Preprocessing:** Quality check, trimming, and alignment.
2. **Differential Expression Analysis:** Case vs. control comparisons for individual datasets.
3. **Meta-Analysis:** Integration of results across multiple studies.
4. **Functional Analysis:** GO term enrichment and visualization.
