############
setwd("/home/ceglab358/BUDDHA/BMCBIO-EoI/Segmental_deletion/STAP1_deletion/CDS_filtered/")

# Load required libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)

# Read the dataset
df <- read.table("final_output.no_gaps.sorted.csv", header = TRUE, sep = "\t", check.names = FALSE)
colnames(df)
colnames(df) <- gsub("-", "_", colnames(df))

cols_to_check <- c( "EPHA5_CENPC",  "EPHA5_UBA6",   "EPHA5_YTHDC1","CENPC_UBA6",   "CENPC_YTHDC1", "UBA6_YTHDC1")

# Function to check if a value is within IQR range
is_within_iqr <- function(value, q1, q3, iqr) {
  (value >= (q1 - 1.5 * iqr)) & (value <= (q3 + 1.5 * iqr))
}

# Compute IQR per Group and filter
df_filtered <- df %>%
  group_by(Group) %>%
  filter(across(all_of(cols_to_check), 
                ~ is_within_iqr(.x, quantile(.x, 0.25), quantile(.x, 0.75), IQR(.x))))
View(df_filtered)
# Display the first few rows of the filtered dataset
head(df_filtered, 10)
str(df)
head(df,10)
# Convert STAP1_status to factor
df_filtered$STAP1_status <- as.factor(df_filtered$STAP1_status)
str(df_filtered)
# Ensure column names are properly formatted (replace dashes with underscores)
colnames(df_filtered) <- gsub("-", "_", colnames(df_filtered))
hist(df_filtered$EPHA5_CENPC)
hist(df$EPHA5_CENPC)

# Define the column names for pairwise distances
pairwise_cols <- c( "EPHA5_CENPC",  "EPHA5_UBA6",   "EPHA5_YTHDC1","CENPC_UBA6",   "CENPC_YTHDC1", "UBA6_YTHDC1")


# Create a matrix to store p-values
p_values_matrix <- matrix(NA, nrow = length(pairwise_cols), ncol = length(pairwise_cols), 
                          dimnames = list(pairwise_cols, pairwise_cols))

# Compute Wilcoxon test p-values for ratios
for (i in 1:length(pairwise_cols)) {
  for (j in 1:length(pairwise_cols)) {
    if (i != j) {  # Avoid self-ratios
      ratio_values <- df_filtered[[pairwise_cols[i]]] / df_filtered[[pairwise_cols[j]]]  # Compute ratio
      ratio_values[is.infinite(ratio_values) | is.nan(ratio_values)] <- NA  # Handle division errors
      
      # Perform Wilcoxon test
      test <- wilcox.test(ratio_values ~ df_filtered$STAP1_status, na.action = na.omit, alternative = "greater")
      p_values_matrix[i, j] <- test$p.value
    }
  }
}
p_values_adj_matrix <- matrix(p.adjust(as.vector(p_values_matrix), method = "fdr"), 
                              nrow = nrow(p_values_matrix), 
                              ncol = ncol(p_values_matrix), 
                              dimnames = dimnames(p_values_matrix))

print(p_values_adj_matrix)

# Plot heatmap of p-values
pheatmap(p_values_adj_matrix, 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Heatmap of p-values (Ratios by IL34 Status)", 
         display_numbers = TRUE, 
         cluster_rows = T, 
         cluster_cols = T, 
         na_col = "grey")  # Mark NA values in grey
