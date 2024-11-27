# Install and load necessary packages
install.packages("metafor")
library(metafor)

# Set the working directory
setwd("D:/setolabo/Metafor")

# Load data from a CSV file
input_data <- read.csv("miRNA_meta.csv")

# Inspect the column names and the first few rows of data (optional)
print(colnames(input_data))
print(head(input_data))

# Remove rows with missing EffectSize or SE values
cleaned_data <- input_data %>%
  filter(!is.na(EffectSize) & !is.na(SE))

# Remove rows with zero or negative SE, and filter by EffectSize and SE
cleaned_data <- cleaned_data %>%
  filter(SE > 0 & SE < 1.0 & (EffectSize >= 0.5 | EffectSize <= -0.5))

# Filter miRNAs detected in at least 2 samples
filtered_data <- cleaned_data %>%
  group_by(miRNA) %>%
  filter(n() >= 2) %>%
  ungroup()

# Write the filtered data to a CSV file
write.csv(filtered_data, file = "filtered_data.csv", row.names = FALSE)

# Create a list to store the results
results_list <- list()

# Loop through each unique miRNA and perform meta-analysis
for (miRNA in unique(filtered_data$miRNA)) {
  
  # Filter the dataset for the current miRNA
  miRNA_data <- filtered_data %>% filter(miRNA == !!miRNA)
  
  # Check if the miRNA has enough studies for meta-analysis (at least 2 data points) (optional)
  if (nrow(miRNA_data) >= 2) {
    
    # Conduct a random-effects meta-analysis using REML
    meta_res <- rma(yi = EffectSize, sei = SE, data = miRNA_data, method = "REML")
    
    # Save the results into the results list
    results_list[[miRNA]] <- list(
      estimate = meta_res$b,
      ci_lb = meta_res$ci.lb,
      ci_ub = meta_res$ci.ub,
      pval = meta_res$pval,
      tau2 = meta_res$tau2,
      i2 = meta_res$I2,
      Q_test = meta_res$QE,
      meta_object = meta_res  # Store the meta-analysis object for further analysis
    )
    
    # Print a summary for each miRNA
    cat("\nMeta-Analysis Results for", miRNA, ":\n")
    print(summary(meta_res))
    
    # Generate and save the forest plot for each miRNA
    pdf(paste0("forest_plot_", miRNA, ".pdf"))
    forest(meta_res, slab = miRNA_data$Study, main = paste("Forest Plot for", miRNA))
    dev.off()
    
  } else {
    # Print message if insufficient studies for this miRNA
    cat("\nNot enough studies for meta-analysis of", miRNA, "\n")
  }
}

# Convert the results list to a data frame for saving
results_df <- do.call(rbind, lapply(names(results_list), function(miRNA) {
  res <- results_list[[miRNA]]
  data.frame(
    miRNA = miRNA,
    estimate = res$estimate,
    ci_lb = res$ci_lb,
    ci_ub = res$ci_ub,
    pval = res$pval,
    tau2 = res$tau2,
    i2 = res$i2,
    Q_test = res$Q_test,
    stringsAsFactors = FALSE
  )
}))

# Save the results to a CSV file
write.csv(results_df, file = "miRNA_meta_analysis_results_filtered.csv", row.names = FALSE)

# Display the final combined results
print("Final results for all miRNAs:")
print(results_df)
