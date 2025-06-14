setwd("./../CDS_filtered/")
getwd()
library(ggplot2)
library(reshape2)
library(pheatmap)

# Read the dataset
df <- read.table("final_output_with_status.tsv", header = TRUE, sep = "\t", check.names = FALSE)
colnames(df)[7:21]
df<-df[df$Num_Gaps == 0, ]
################## HEATMAPS ####################
df$Gene_status <- as.factor(df$Gene_status)
pairwise_cols <- colnames(df)[7:21]



p_values_matrix <- matrix(NA, nrow = length(pairwise_cols), ncol = length(pairwise_cols), 
                          dimnames = list(pairwise_cols, pairwise_cols))

# Compute Wilcoxon test p-values for ratios
for (i in 1:length(pairwise_cols)) {
  for (j in 1:length(pairwise_cols)) {
    if (i != j) {  # Avoid self-ratios
      ratio_values <- df[[pairwise_cols[i]]] / df[[pairwise_cols[j]]]  # Compute ratio
      ratio_values[is.infinite(ratio_values) | is.nan(ratio_values)] <- NA  # Handle division errors
      
      # Perform Wilcoxon test
      test <- wilcox.test(ratio_values ~ df$Gene_status, na.action = na.omit, alternative = "less")
      p_values_matrix[i, j] <- test$p.value
    }
  }
}
p_values_adj_matrix <- matrix(p.adjust(as.vector(p_values_matrix), method = "fdr"), 
                              nrow = nrow(p_values_matrix), 
                              ncol = ncol(p_values_matrix), 
                              dimnames = dimnames(p_values_matrix))



# Plot heatmap of p-values
png("Heatmap_of_deletion_plot.png",height=6,width=10,units = "in",res=600)
pheatmap(p_values_adj_matrix, 
         color = colorRampPalette(c("yellow","skyblue","lightblue", "lightcyan", "white", "pink", "lightpink"))(50), 
         display_numbers = TRUE, 
         cluster_rows = F, 
         cluster_cols = F, 
         na_col = "grey") 
dev.off()
