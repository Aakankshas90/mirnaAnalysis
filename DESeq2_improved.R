# Differential expression analysis using DESeq2
# Multiple independent miRNA-seq datasets

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)

setwd("D:/setolabo/DESeq2")

# List of datasets
datasets <- c(
  "PRJNA956559",
  "PRJNA946800",
  "PRJNA540915",
  "PRJNA350003",
  "PRJNA347596",
  "PRJNA866134"
)

for (dataset in datasets) {

  cat("Processing:", dataset, "\n")

  # Load data
  count_data <- read.csv(paste0(dataset, "counts.csv"), row.names = 1)
  sample_data <- read.csv(paste0(dataset, ".csv"), row.names = 1)

  # Create DESeq object
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_data,
    design = ~ condition
  )

  # Run DESeq2
  dds <- DESeq(dds)

  # Save results
  res <- results(dds)
  write.csv(as.data.frame(res), file = paste0(dataset, "_results.csv"))

  # rlog transformation
  rld <- rlog(dds, blind = FALSE)

  # PCA
  pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
    geom_point(size = 3) +
    geom_text(vjust = -1, size = 2) +
    labs(
      title = paste("PCA:", dataset),
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance")
    ) +
    theme_minimal()

  ggsave(paste0(dataset, "_PCA.png"), plot = p, width = 6, height = 5)

}
