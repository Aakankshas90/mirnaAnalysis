if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)
setwd("D:\\setolabo\\DESeq2")

#111111111111111111111111111111111111111111111111111111111111111111111111111

count_data1 <- read.csv("PRJNA956559counts.csv", row.names=1)
sample_data1 <- read.csv("PRJNA956559.csv", row.names=1)

# Create DESeqDataSet object
dds1 <- DESeqDataSetFromMatrix(countData = count_data1,
                               colData = sample_data1,
                               design = ~ condition)

# Run DESeq2 pipeline
dds1 <- DESeq(dds1)

# Extract results
results1 <- results(dds1)

# Save results to a file
write.csv(as.data.frame(results1), file = "PRJNA956559_results.csv")

# use rlog transformation (slower, better for smaller datasets)
rld1 <- rlog(dds1, blind=FALSE)

# Extract PCA data with sample labels
pca_data1 <- plotPCA(rld1, intgroup="condition", returnData=TRUE)

# Plot PCA
ggplot(pca_data1, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -1, hjust = 0.5, size = 1.5) +  # Add sample labels
  labs(title = "PCA Plot", x = paste0("PC1: ", round(100 * attr(pca_data, "percentVar")[1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * attr(pca_data, "percentVar")[2], 1), "% variance")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),  # Adjust x-axis title font size
    axis.title.y = element_text(size = 16),  # Adjust y-axis title font size
    axis.text = element_text(size = 14),     # Adjust axis tick labels font size
    plot.title = element_text(size = 10, face = "bold")  # Adjust title font size
  )

#22222222222222222222222222222222222222222222222222222222222222222222

count_data2 <- read.csv("PRJNA946800counts.csv", row.names=1)
sample_data2 <- read.csv("PRJNA946800.csv", row.names=1)

# Create DESeqDataSet object
dds2 <- DESeqDataSetFromMatrix(countData = count_data2,
                               colData = sample_data2,
                               design = ~ condition)

# Run DESeq2 pipeline
dds2 <- DESeq(dds2)

# Extract results
results2 <- results(dds2)

# Save results to a file
write.csv(as.data.frame(results2), file = "PRJNA946800_results.csv")

# use rlog transformation (slower, better for smaller datasets)
rld2 <- rlog(dds2, blind=FALSE)

# Extract PCA data with sample labels
pca_data2 <- plotPCA(rld2, intgroup="condition", returnData=TRUE)

# Plot PCA
ggplot(pca_data2, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -1, hjust = 0.5, size = 1.5) +  # Add sample labels
  labs(title = "PCA Plot", x = paste0("PC1: ", round(100 * attr(pca_data, "percentVar")[1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * attr(pca_data, "percentVar")[2], 1), "% variance")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),  # Adjust x-axis title font size
    axis.title.y = element_text(size = 16),  # Adjust y-axis title font size
    axis.text = element_text(size = 14),     # Adjust axis tick labels font size
    plot.title = element_text(size = 10, face = "bold")  # Adjust title font size
  )

#3333333333333333333333333333333333333333333333333333333333333333333333333333

count_data3 <- read.csv("PRJNA540915counts.csv", row.names=1)
sample_data3 <- read.csv("PRJNA540915.csv", row.names=1)

# Create DESeqDataSet object
dds3 <- DESeqDataSetFromMatrix(countData = count_data3,
                               colData = sample_data3,
                               design = ~ condition)

# Run DESeq2 pipeline
dds3 <- DESeq(dds3)

# Extract results
results3 <- results(dds3)

# Save results to a file
write.csv(as.data.frame(results3), file = "PRJNA540915_results.csv")

# use rlog transformation (slower, better for smaller datasets)
rld3 <- rlog(dds3, blind=FALSE)

# Extract PCA data with sample labels
pca_data3 <- plotPCA(rld3, intgroup="condition", returnData=TRUE)

# Plot PCA
ggplot(pca_data3, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -1, hjust = 0.5, size = 1.5) +  # Add sample labels
  labs(title = "PCA Plot", x = paste0("PC1: ", round(100 * attr(pca_data, "percentVar")[1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * attr(pca_data, "percentVar")[2], 1), "% variance")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),  # Adjust x-axis title font size
    axis.title.y = element_text(size = 16),  # Adjust y-axis title font size
    axis.text = element_text(size = 14),     # Adjust axis tick labels font size
    plot.title = element_text(size = 10, face = "bold")  # Adjust title font size
  )

#44444444444444444444444444444444444444444444444444444444444444444444444

count_data4 <- read.csv("PRJNA350003counts.csv", row.names=1)
sample_data4 <- read.csv("PRJNA350003.csv", row.names=1)

# Create DESeqDataSet object
dds4 <- DESeqDataSetFromMatrix(countData = count_data4,
                               colData = sample_data4,
                               design = ~ condition)

# Run DESeq2 pipeline
dds4 <- DESeq(dds4)

# Extract results
results4 <- results(dds4)

# Save results to a file
write.csv(as.data.frame(results4), file = "PRJNA350003_results.csv")

# use rlog transformation (slower, better for smaller datasets)
rld4 <- rlog(dds4, blind=FALSE)

# Extract PCA data with sample labels
pca_data4 <- plotPCA(rld4, intgroup="condition", returnData=TRUE)

# Plot PCA
ggplot(pca_data4, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -1, hjust = 0.5, size = 1.5) +  # Add sample labels
  labs(title = "PCA Plot", x = paste0("PC1: ", round(100 * attr(pca_data, "percentVar")[1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * attr(pca_data, "percentVar")[2], 1), "% variance")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),  # Adjust x-axis title font size
    axis.title.y = element_text(size = 16),  # Adjust y-axis title font size
    axis.text = element_text(size = 14),     # Adjust axis tick labels font size
    plot.title = element_text(size = 10, face = "bold")  # Adjust title font size
  )

#55555555555555555555555555555555555555555555555555555555555555555555555

count_data5 <- read.csv("PRJNA347596counts.csv", row.names=1)
sample_data5 <- read.csv("PRJNA347596.csv", row.names=1)

# Create DESeqDataSet object
dds5 <- DESeqDataSetFromMatrix(countData = count_data5,
                               colData = sample_data5,
                               design = ~ condition)

# Run DESeq2 pipeline
dds5 <- DESeq(dds5)

# Extract results
results5 <- results(dds5)

# Save results to a file
write.csv(as.data.frame(results5), file = "PRJNA347596_results.csv")

# use rlog transformation (slower, better for smaller datasets)
rld5 <- rlog(dds5, blind=FALSE)

# Extract PCA data with sample labels
pca_data5 <- plotPCA(rld5, intgroup="condition", returnData=TRUE)

# Plot PCA
ggplot(pca_data5, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -1, hjust = 0.5, size = 1.5) +  # Add sample labels
  labs(title = "PCA Plot", x = paste0("PC1: ", round(100 * attr(pca_data, "percentVar")[1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * attr(pca_data, "percentVar")[2], 1), "% variance")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),  # Adjust x-axis title font size
    axis.title.y = element_text(size = 16),  # Adjust y-axis title font size
    axis.text = element_text(size = 14),     # Adjust axis tick labels font size
    plot.title = element_text(size = 10, face = "bold")  # Adjust title font size
  )

#6666666666666666666666666666666666666666666666666666666666666666666666666

count_data6 <- read.csv("PRJNA866134counts.csv", row.names=1)
sample_data6 <- read.csv("PRJNA866134.csv", row.names=1)

# Create DESeqDataSet object
dds6 <- DESeqDataSetFromMatrix(countData = count_data6,
                               colData = sample_data6,
                               design = ~ condition)

# Run DESeq2 pipeline
dds6 <- DESeq(dds6)

# Extract results
results6 <- results(dds6)

# Save results to a file
write.csv(as.data.frame(results6), file = "PRJNA866134_results.csv")

# use rlog transformation (slower, better for smaller datasets)
rld6 <- rlog(dds6, blind=FALSE)

# Extract PCA data with sample labels
pca_data6 <- plotPCA(rld6, intgroup="condition", returnData=TRUE)

# Plot PCA
ggplot(pca_data6, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -1, hjust = 0.5, size = 1.5) +  # Add sample labels
  labs(title = "PCA Plot", x = paste0("PC1: ", round(100 * attr(pca_data, "percentVar")[1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * attr(pca_data, "percentVar")[2], 1), "% variance")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),  # Adjust x-axis title font size
    axis.title.y = element_text(size = 16),  # Adjust y-axis title font size
    axis.text = element_text(size = 14),     # Adjust axis tick labels font size
    plot.title = element_text(size = 10, face = "bold")  # Adjust title font size
  )

