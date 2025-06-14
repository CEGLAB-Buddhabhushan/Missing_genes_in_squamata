setwd("/home/ceglab358/BUDDHA/BMCBIO-EoI/Segmental_deletion/IL34_deletion/CDS_filtered//")


##########

library(ggplot2)
library(reshape2)
library(pheatmap)

# Read the dataset
df <- read.table("final_output.no_gaps.sorted.tsv", header = TRUE, sep = "\t", check.names = FALSE)
colnames(df)
colnames(df) <- gsub("-", "_", colnames(df))

cols_to_check <- c("FCSK_COG4", "FCSK_SF3B3", "FCSK_MTSS2", "FCSK_VAC14", 
                   "FCSK_HYDIN", "COG4_SF3B3", "COG4_MTSS2", "COG4_VAC14", 
                   "COG4_HYDIN", "SF3B3_MTSS2", "SF3B3_VAC14", "SF3B3_HYDIN", 
                   "MTSS2_VAC14", "MTSS2_HYDIN", "VAC14_HYDIN")

cols_to_check <- c("COG4_SF3B3",   "COG4_MTSS2",   "COG4_VAC14",  
                    "SF3B3_MTSS2",  "SF3B3_VAC14",  "MTSS2_VAC14")
cols_to_check <- c("FCSK_SF3B3",   "FCSK_MTSS2",   "FCSK_HYDIN",  
                    "SF3B3_MTSS2",  "SF3B3_HYDIN",  "MTSS2_HYDIN")

# Function to check if a value is within IQR range
is_within_iqr <- function(value, q1, q3, iqr) {
  (value >= (q1 - 0.5 * iqr)) & (value <= (q3 + 0.5 * iqr))
}
#df_filtered <- df
# Compute IQR per Group and filter
df_filtered <- df %>%
  group_by(Group) %>%
  filter(across(all_of(cols_to_check), 
                ~ is_within_iqr(.x, quantile(.x, 0.25), quantile(.x, 0.75), IQR(.x))))

# Convert STAP1_status to factor
df_filtered$IL34_status <- as.factor(df_filtered$IL34_status)
r1< df_filtered$
var(df_filtered)

  
#View(df)
#View(df_filtered)
# Ensure column names are properly formatted (replace dashes with underscores)
hist(df_filtered$COG4_MTSS2)
max(df_filtered$VAC14_HYDIN)
R1<- df_filtered$COG4_MTSS2/df_filtered$COG4_VAC14
var(R1)
#hist(df$FCSK_HYDIN)
head(df,2)
# Define the column names for pairwise distances
pairwise_cols <- c("FCSK_COG4", "FCSK_SF3B3", "FCSK_MTSS2", "FCSK_VAC14", 
                   "FCSK_HYDIN", "COG4_SF3B3", "COG4_MTSS2", "COG4_VAC14", 
                   "COG4_HYDIN", "SF3B3_MTSS2", "SF3B3_VAC14", "SF3B3_HYDIN", 
                   "MTSS2_VAC14", "MTSS2_HYDIN", "VAC14_HYDIN")
pairwise_cols <- c("COG4_SF3B3",   "COG4_MTSS2",   "COG4_VAC14",  
                   "SF3B3_MTSS2",  "SF3B3_VAC14",  "MTSS2_VAC14")
pairwise_cols <- c("FCSK_SF3B3",   "FCSK_MTSS2",   "FCSK_HYDIN",  
                   "SF3B3_MTSS2",  "SF3B3_HYDIN",  "MTSS2_HYDIN")
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
      test <- wilcox.test(ratio_values ~ df_filtered$IL34_status, na.action = na.omit, alternative = "less")
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
         color = colorRampPalette(c("yellow","skyblue","purple","blue", "white", "red", "pink"))(50), 
         main = "Heatmap of p-values (Ratios by IL34 Status)", 
         display_numbers = TRUE, 
         cluster_rows = T, 
         cluster_cols = T, 
         na_col = "grey")  # Mark NA values in grey


library(pheatmap)

# Define a custom matrix to display only values < 0.01
display_matrix <- ifelse(p_values_adj_matrix < 0.01, round(p_values_matrix, 4), "")

# Generate heatmap
pheatmap(p_values_adj_matrix, 
         color = colorRampPalette(c("yellow", "skyblue", "purple", "blue", "white", "red", "pink"))(100), 
         main = "Heatmap of p-values (Ratios by IL34 Status)", 
         display_numbers = display_matrix,  # Display only values < 0.01
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         na_col = "grey") 

