setwd("/home/ceglab358/BUDDHA/BMCBIO-EoI/STAP1_deletion/CDS_filtered/")
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr) 
library(patchwork)
library(gridExtra)
library(cowplot)
# Load data
df1 <- read.table("final_output.no_gaps.sorted.tsv", header = TRUE, sep = "\t")
colnames(df1)
# Identify numeric columns (excluding "gene_status")
numeric_cols1 <- setdiff(names(df1)[sapply(df1, is.numeric)], c("Group", "Order", "Species_name", "STAP1_status"))
df1 <- df1[!(df1$Group %in% c("Chondrichthyes")), ]
head(df1,10)
df1$CENPC.UBA6_EPHA5.CENPC <- df1$CENPC.UBA6/df1$EPHA5.CENPC
df1$CENPC.UBA6_UBA6.YTHDC1 <- df1$CENPC.UBA6/df1$UBA6.YTHDC1
df1$EPHA5.CENPC_UBA6.YTHDC1 <- df1$EPHA5.CENPC/df1$UBA6.YTHDC1


df_melt_CENPC.UBA6_EPHA5.CENPC <- melt(df1[,c("Group", "Species_name", "STAP1_status", "CENPC.UBA6_EPHA5.CENPC")], id.vars = c("Group", "Species_name", "STAP1_status"))
df_melt_CENPC.UBA6_UBA6.YTHDC1 <- melt(df1[,c("Group", "Species_name", "STAP1_status", "CENPC.UBA6_UBA6.YTHDC1")], id.vars = c("Group", "Species_name", "STAP1_status"))
df_melt_EPHA5.CENPC_UBA6.YTHDC1 <- melt(df1[,c("Group", "Species_name", "STAP1_status", "EPHA5.CENPC_UBA6.YTHDC1")], id.vars = c("Group", "Species_name", "STAP1_status"))

# Create boxplot
# Create boxplot with significance annotations
##STAP1_gene_status
T1gene_p_df_melt_CENPC.UBA6_EPHA5.CENPC<-ggplot(df_melt_CENPC.UBA6_EPHA5.CENPC, aes(x = factor(STAP1_status), y = value)) +
  geom_boxplot(aes(fill = factor(STAP1_status))) + 
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "less"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Gene Status", y = "Gene distance ratio", fill = "Gene Status") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom")  # Move legend to bottom
wilcox.test(value ~ STAP1_status, data = df_melt_EPHA5.CENPC_UBA6.YTHDC1)
head (df_melt_CENPC.UBA6_EPHA5.CENPC,5)
T1gene_p_df_melt_CENPC.UBA6_UBA6.YTHDC1<-ggplot(df_melt_CENPC.UBA6_UBA6.YTHDC1, aes(x = factor(STAP1_status), y = value)) +
  geom_boxplot(aes(fill = factor(STAP1_status))) +
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "less"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Gene Status", y = "Gene distance ratio", fill = "Gene Status") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom")

T1gene_p_df_melt_EPHA5.CENPC_UBA6.YTHDC1<-ggplot(df_melt_EPHA5.CENPC_UBA6.YTHDC1, aes(x = factor(STAP1_status), y = value)) +
  geom_boxplot(aes(fill = factor(STAP1_status))) + 
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Gene Status", y = "Gene distance ratio", fill = "Gene Status") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom")

T1Gp1 <- T1gene_p_df_melt_CENPC.UBA6_EPHA5.CENPC + theme(legend.position = "none")
T1Gp2 <- T1gene_p_df_melt_CENPC.UBA6_UBA6.YTHDC1 + theme(legend.position = "none")
T1Gp3 <- T1gene_p_df_melt_EPHA5.CENPC_UBA6.YTHDC1 + theme(legend.position = "none")

legend <- get_legend(
  T1gene_p_df_melt_CENPC.UBA6_EPHA5.CENPC + 
    theme(legend.position = "bottom") + 
    guides(fill = guide_legend(nrow = 1))  # Force legend to be in one row
)
# Arrange plots with the common legend
Gene_T1_plot<- grid.arrange(
  arrangeGrob(T1Gp1, T1Gp2, T1Gp3, ncol = 3),  # Arrange plots in a column
  legend,  # Add common legend below
  ncol = 1,  # Stack vertically
  heights = c(3, 0.2)  # Adjust legend size
)

ggsave(filename = "STAP1_Gene_plot.png", plot = Gene_T1_plot, width = 16, height = 9,dpi = 600) 


############
setwd("/home/ceglab358/BUDDHA/BMCBIO-EoI/IL34_Rearrangement/Segmental_deletion/gene_distance/")

# Load required libraries
library(ggplot2)
library(reshape2)
library(pheatmap)

# Read the dataset
df <- read.table("final_output.filtered.heatmap.removed_rearrengement.tsv", header = TRUE, sep = "\t", check.names = FALSE)
colnames(df)
# Convert STAP1_status to factor
df$IL34_status <- as.factor(df$IL34_status)
str(df)
# Ensure column names are properly formatted (replace dashes with underscores)
colnames(df) <- gsub("-", "_", colnames(df))
hist(df$FCSK_HYDIN)
max(df$FCSK_HYDIN)
head(df,2)
# Define the column names for pairwise distances
pairwise_cols <- c("FCSK_COG4", "FCSK_SF3B3", "FCSK_MTSS2", "FCSK_VAC14", 
                   "FCSK_HYDIN", "COG4_SF3B3", "COG4_MTSS2", "COG4_VAC14", 
                   "COG4_HYDIN", "SF3B3_MTSS2", "SF3B3_VAC14", "SF3B3_HYDIN", 
                   "MTSS2_VAC14", "MTSS2_HYDIN", "VAC14_HYDIN")
pairwise_cols <- c("FCSK_SF3B3","COG4_SF3B3", "SF3B3_MTSS2", 
                   "MTSS2_VAC14", "MTSS2_HYDIN")
# Create a matrix to store p-values
p_values_matrix <- matrix(NA, nrow = length(pairwise_cols), ncol = length(pairwise_cols), 
                          dimnames = list(pairwise_cols, pairwise_cols))

# Compute Wilcoxon test p-values for ratios
for (i in 1:length(pairwise_cols)) {
  for (j in 1:length(pairwise_cols)) {
    if (i != j) {  # Avoid self-ratios
      ratio_values <- df[[pairwise_cols[i]]] / df[[pairwise_cols[j]]]  # Compute ratio
      ratio_values[is.infinite(ratio_values) | is.nan(ratio_values)] <- NA  # Handle division errors
      
      # Perform Wilcoxon test
      test <- wilcox.test(ratio_values ~ df$IL34_status, na.action = na.omit)
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


##########

library(ggplot2)
library(reshape2)
library(pheatmap)

# Read the dataset
df <- read.table("final_output.filtered.heatmap.tsv", header = TRUE, sep = "\t", check.names = FALSE)
colnames(df)
colnames(df) <- gsub("-", "_", colnames(df))

cols_to_check <- c("FCSK_COG4", "FCSK_SF3B3", "FCSK_MTSS2", "FCSK_VAC14", 
                   "FCSK_HYDIN", "COG4_SF3B3", "COG4_MTSS2", "COG4_VAC14", 
                   "COG4_HYDIN", "SF3B3_MTSS2", "SF3B3_VAC14", "SF3B3_HYDIN", 
                   "MTSS2_VAC14", "MTSS2_HYDIN", "VAC14_HYDIN")

# Function to check if a value is within IQR range
is_within_iqr <- function(value, q1, q3, iqr) {
  (value >= (q1 - 1.5 * iqr)) & (value <= (q3 + 1.5 * iqr))
}

# Compute IQR per Group and filter
df_filtered <- df %>%
  group_by(Group) %>%
  filter(across(all_of(cols_to_check), 
                ~ is_within_iqr(.x, quantile(.x, 0.25), quantile(.x, 0.75), IQR(.x))))

# Convert STAP1_status to factor
df_filtered$IL34_status <- as.factor(df_filtered$IL34_status)
str(df_filtered)
View(df)
View(df_filtered)
# Ensure column names are properly formatted (replace dashes with underscores)
hist(df_filtered$VAC14_HYDIN)
hist(df$FCSK_HYDIN)
head(df,2)
# Define the column names for pairwise distances
pairwise_cols <- c("FCSK_COG4", "FCSK_SF3B3", "FCSK_MTSS2", "FCSK_VAC14", 
                   "FCSK_HYDIN", "COG4_SF3B3", "COG4_MTSS2", "COG4_VAC14", 
                   "COG4_HYDIN", "SF3B3_MTSS2", "SF3B3_VAC14", "SF3B3_HYDIN", 
                   "MTSS2_VAC14", "MTSS2_HYDIN", "VAC14_HYDIN")

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
      test <- wilcox.test(ratio_values ~ df_filtered$IL34_status, na.action = na.omit)
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
pheatmap(p_values_matrix, 
         color = colorRampPalette(c("yellow","skyblue","purple","blue", "white", "red", "pink"))(50), 
         main = "Heatmap of p-values (Ratios by IL34 Status)", 
         display_numbers = TRUE, 
         cluster_rows = T, 
         cluster_cols = T, 
         na_col = "grey")  # Mark NA values in grey


library(pheatmap)

# Define a custom matrix to display only values < 0.01
display_matrix <- ifelse(p_values_adj_matrix < 0.05, round(p_values_matrix, 4), "")

# Generate heatmap
pheatmap(p_values_adj_matrix, 
         color = colorRampPalette(c("yellow", "skyblue", "purple", "blue", "white", "red", "pink"))(100), 
         main = "Heatmap of p-values (Ratios by IL34 Status)", 
         display_numbers = display_matrix,  # Display only values < 0.01
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         na_col = "grey") 

