# Install necessary Bioconductor packages
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

setwd("D:/setolabo/Metafor/4 datasets meta-analysis")

# Step 1: Read the CSV file
gene_data <- read.csv("target_genes_from_miRNAs.csv", header = TRUE)

# Step 2: Extract the gene symbols from the file
genes <- gene_data$Gene.Symbol

# Step 3: Convert gene symbols to Entrez IDs
gene_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Step 4: Perform GO enrichment analysis (Biological Process ontology)
go_enrich <- enrichGO(gene = gene_ids$ENTREZID, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID", 
                      ont = "BP",  # Ontology: BP = Biological Process, CC = Cellular Component, MF = Molecular Function
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.05)

# Step 5: Display the top enriched GO terms
head(go_enrich)

# Step 6: Visualize the results with a dot plot
dotplot(go_enrich, showCategory = 20) + ggtitle("GO Enrichment Analysis (Biological Process)")
