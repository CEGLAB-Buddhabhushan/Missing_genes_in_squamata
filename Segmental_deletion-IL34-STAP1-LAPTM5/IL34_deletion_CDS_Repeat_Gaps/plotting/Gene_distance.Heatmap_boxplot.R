setwd("/home/ceglab358/BUDDHA/BMCBIO-EoI/Segmental_deletion/IL34_deletion_CDS_Repeat_Gaps/plotting//")
library(ggplot2)
library(reshape2)
library(pheatmap)

# Read the dataset
df <- read.table("final_output.sorted.tsv", header = TRUE, sep = "\t", check.names = FALSE)
colnames(df)
################## HEATMAPS ####################
##### Without IQR, with gaps and with repeats
df$IL34_status <- as.factor(df$IL34_status)
pairwise_cols <- c("COG4_SF3B3",   "COG4_MTSS2",   "COG4_VAC14",  
                   "SF3B3_MTSS2",  "SF3B3_VAC14",  "MTSS2_VAC14")

p_values_matrix <- matrix(NA, nrow = length(pairwise_cols), ncol = length(pairwise_cols), 
                          dimnames = list(pairwise_cols, pairwise_cols))

# Compute Wilcoxon test p-values for ratios
for (i in 1:length(pairwise_cols)) {
  for (j in 1:length(pairwise_cols)) {
    if (i != j) {  # Avoid self-ratios
      ratio_values <- df[[pairwise_cols[i]]] / df[[pairwise_cols[j]]]  # Compute ratio
      ratio_values[is.infinite(ratio_values) | is.nan(ratio_values)] <- NA  # Handle division errors
      
      # Perform Wilcoxon test
      test <- wilcox.test(ratio_values ~ df$IL34_status, na.action = na.omit, alternative = "less")
      p_values_matrix[i, j] <- test$p.value
    }
  }
}
p_values_adj_matrix <- matrix(p.adjust(as.vector(p_values_matrix), method = "fdr"), 
                              nrow = nrow(p_values_matrix), 
                              ncol = ncol(p_values_matrix), 
                              dimnames = dimnames(p_values_matrix))


# Plot heatmap of p-values
without_IQR_with_gaps_with_repeats <- pheatmap(p_values_adj_matrix, 
                                               color = colorRampPalette(c("yellow","skyblue","lightblue", "lightcyan", "white", "pink", "lightpink"))(50), 
                                               display_numbers = TRUE, 
         cluster_rows = F, 
         cluster_cols = F, 
         na_col = "grey") 


##### within IQR, with gaps, with repeat
# Function to check if a value is within IQR range
is_within_iqr <- function(value, q1, q3, iqr) {
  (value >= (q1 - 0.5 * iqr)) & (value <= (q3 + 0.5 * iqr))
}
cols_to_check <- c("COG4_SF3B3",   "COG4_MTSS2",   "COG4_VAC14",  
                   "SF3B3_MTSS2",  "SF3B3_VAC14",  "MTSS2_VAC14")
df_filtered_IQR <- df %>%
  group_by(Group) %>%
  filter(across(all_of(cols_to_check), 
                ~ is_within_iqr(.x, quantile(.x, 0.25), quantile(.x, 0.75), IQR(.x))))

df_filtered_IQR$IL34_status <- as.factor(df_filtered_IQR$IL34_status)
pairwise_cols <- c("COG4_SF3B3",   "COG4_MTSS2",   "COG4_VAC14",  
                   "SF3B3_MTSS2",  "SF3B3_VAC14",  "MTSS2_VAC14")

p_values_matrix <- matrix(NA, nrow = length(pairwise_cols), ncol = length(pairwise_cols), 
                          dimnames = list(pairwise_cols, pairwise_cols))

# Compute Wilcoxon test p-values for ratios
for (i in 1:length(pairwise_cols)) {
  for (j in 1:length(pairwise_cols)) {
    if (i != j) {  # Avoid self-ratios
      ratio_values <- df_filtered_IQR[[pairwise_cols[i]]] / df_filtered_IQR[[pairwise_cols[j]]]  # Compute ratio
      ratio_values[is.infinite(ratio_values) | is.nan(ratio_values)] <- NA  # Handle division errors
      
      # Perform Wilcoxon test
      test <- wilcox.test(ratio_values ~ df_filtered_IQR$IL34_status, na.action = na.omit, alternative = "less")
      p_values_matrix[i, j] <- test$p.value
    }
  }
}
p_values_adj_matrix <- matrix(p.adjust(as.vector(p_values_matrix), method = "fdr"), 
                              nrow = nrow(p_values_matrix), 
                              ncol = ncol(p_values_matrix), 
                              dimnames = dimnames(p_values_matrix))


# Plot heatmap of p-values
with_IQR_with_gaps_with_repeats <- pheatmap(p_values_adj_matrix, 
                                            color = colorRampPalette(c("yellow","skyblue","lightblue", "lightcyan", "white", "pink", "lightpink"))(50), 
                                            display_numbers = TRUE, 
                                            cluster_rows = F, 
                                            cluster_cols = F, 
                                            na_col = "grey") 
##### within IQR, No gaps, with repeat
df_filtered_IQR_no_gaps <- df_filtered_IQR[df_filtered_IQR$Num_Gaps == 0, ]
df_filtered_IQR_no_gaps$IL34_status <- as.factor(df_filtered_IQR_no_gaps$IL34_status)
pairwise_cols <- c("COG4_SF3B3",   "COG4_MTSS2",   "COG4_VAC14",  
                   "SF3B3_MTSS2",  "SF3B3_VAC14",  "MTSS2_VAC14")

p_values_matrix <- matrix(NA, nrow = length(pairwise_cols), ncol = length(pairwise_cols), 
                          dimnames = list(pairwise_cols, pairwise_cols))

# Compute Wilcoxon test p-values for ratios
for (i in 1:length(pairwise_cols)) {
  for (j in 1:length(pairwise_cols)) {
    if (i != j) {  # Avoid self-ratios
      ratio_values <- df_filtered_IQR_no_gaps[[pairwise_cols[i]]] / df_filtered_IQR_no_gaps[[pairwise_cols[j]]]  # Compute ratio
      ratio_values[is.infinite(ratio_values) | is.nan(ratio_values)] <- NA  # Handle division errors
      
      # Perform Wilcoxon test
      test <- wilcox.test(ratio_values ~ df_filtered_IQR_no_gaps$IL34_status, na.action = na.omit, alternative = "less")
      p_values_matrix[i, j] <- test$p.value
    }
  }
}
p_values_adj_matrix <- matrix(p.adjust(as.vector(p_values_matrix), method = "fdr"), 
                              nrow = nrow(p_values_matrix), 
                              ncol = ncol(p_values_matrix), 
                              dimnames = dimnames(p_values_matrix))


# Plot heatmap of p-values
With_IQR_No_gaps_with_repeats <- pheatmap(p_values_adj_matrix, 
                                          color = colorRampPalette(c("yellow","skyblue","lightblue", "lightcyan", "white", "pink", "lightpink"))(50), 
                                          display_numbers = TRUE, 
                                          cluster_rows = F, 
                                          cluster_cols = F, 
                                          na_col = "grey") 

#####  within IQR, No gaps, No repeat
df_filtered_IQR_no_gaps$COG4_SF3B3<- df_filtered_IQR_no_gaps$COG4_SF3B3-df_filtered_IQR_no_gaps$Repeat_COG4_SF3B3
df_filtered_IQR_no_gaps$COG4_MTSS2<- df_filtered_IQR_no_gaps$COG4_MTSS2-df_filtered_IQR_no_gaps$Repeat_COG4_MTSS2
df_filtered_IQR_no_gaps$COG4_VAC14<- df_filtered_IQR_no_gaps$COG4_VAC14-df_filtered_IQR_no_gaps$Repeat_COG4_VAC14
df_filtered_IQR_no_gaps$SF3B3_MTSS2<- df_filtered_IQR_no_gaps$SF3B3_MTSS2-df_filtered_IQR_no_gaps$Repeat_SF3B3_MTSS2
df_filtered_IQR_no_gaps$SF3B3_VAC14<- df_filtered_IQR_no_gaps$SF3B3_VAC14-df_filtered_IQR_no_gaps$Repeat_SF3B3_VAC14
df_filtered_IQR_no_gaps$MTSS2_VAC14<- df_filtered_IQR_no_gaps$MTSS2_VAC14-df_filtered_IQR_no_gaps$Repeat_MTSS2_VAC14

df_filtered_IQR_no_gaps$IL34_status <- as.factor(df_filtered_IQR_no_gaps$IL34_status)
colnames(df_filtered_IQR_no_gaps)
pairwise_cols <- c("COG4_SF3B3",   "COG4_MTSS2",   "COG4_VAC14",  
                   "SF3B3_MTSS2",  "SF3B3_VAC14",  "MTSS2_VAC14")

p_values_matrix <- matrix(NA, nrow = length(pairwise_cols), ncol = length(pairwise_cols), 
                          dimnames = list(pairwise_cols, pairwise_cols))

# Compute Wilcoxon test p-values for ratios
for (i in 1:length(pairwise_cols)) {
  for (j in 1:length(pairwise_cols)) {
    if (i != j) {  # Avoid self-ratios
      ratio_values <- df_filtered_IQR_no_gaps[[pairwise_cols[i]]] / df_filtered_IQR_no_gaps[[pairwise_cols[j]]]  # Compute ratio
      ratio_values[is.infinite(ratio_values) | is.nan(ratio_values)] <- NA  # Handle division errors
      
      # Perform Wilcoxon test
      test <- wilcox.test(ratio_values ~ df_filtered_IQR_no_gaps$IL34_status, na.action = na.omit, alternative = "less")
      p_values_matrix[i, j] <- test$p.value
    }
  }
}
p_values_adj_matrix <- matrix(p.adjust(as.vector(p_values_matrix), method = "fdr"), 
                              nrow = nrow(p_values_matrix), 
                              ncol = ncol(p_values_matrix), 
                              dimnames = dimnames(p_values_matrix))

# Plot heatmap of p-values
With_IQR_No_gaps_No_repeats <- pheatmap(p_values_adj_matrix, 
                                        color = colorRampPalette(c("yellow","skyblue","lightblue", "lightcyan", "white", "pink", "lightpink"))(50), 
                                        display_numbers = TRUE, 
                                        cluster_rows = F, 
                                        cluster_cols = F, 
                                        na_col = "grey")



##### All heatmaps save
png("without_IQR_with_gaps_with_repeats.png", width = 4 , height = 4, units = "in", res = 600)
print(without_IQR_with_gaps_with_repeats)
dev.off()

png("with_IQR_with_gaps_with_repeats.png", width = 4 , height = 4, units = "in", res = 600)
print(with_IQR_with_gaps_with_repeats)
dev.off()
          
png("With_IQR_No_gaps_with_repeats.png", width = 4 , height = 4, units = "in", res = 600)
print(With_IQR_No_gaps_with_repeats)
dev.off() 

png("With_IQR_No_gaps_No_repeats.png", width = 4 , height = 4, units = "in", res = 600)
print(With_IQR_No_gaps_No_repeats)
dev.off()

################## BOXPLOTS ####################
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr) 
library(patchwork)
library(gridExtra)
library(cowplot)
# Load data
df <- read.table("final_output.sorted.tsv", header = TRUE, sep = "\t")
df_filtered_no_gaps <- df[df$Num_Gaps == 0, ]
# Identify numeric columns (excluding "gene_status")
numeric_cols1 <- setdiff(names(df_filtered_no_gaps)[sapply(df_filtered_no_gaps, is.numeric)], c("Group", "Species_name", "IL34_status"))
df_filtered_no_gaps <- df_filtered_no_gaps[!(df_filtered_no_gaps$Group %in% c("Actinopteri", "Crocodilia")), ]
colnames(df_filtered_no_gaps)


df_filtered_no_gaps$SF3B3_MTSS2.COG4_SF3B3 <- df_filtered_no_gaps$SF3B3_MTSS2/df_filtered_no_gaps$COG4_SF3B3
df_filtered_no_gaps$SF3B3_MTSS2.MTSS2_VAC14 <- df_filtered_no_gaps$SF3B3_MTSS2/df_filtered_no_gaps$MTSS2_VAC14
df_filtered_no_gaps$COG4_SF3B3.MTSS2_VAC14 <- df_filtered_no_gaps$COG4_SF3B3/df_filtered_no_gaps$MTSS2_VAC14



df_melt_SF3B3_MTSS2.COG4_SF3B3 <- melt(df_filtered_no_gaps[,c("Group", "Species_name", "IL34_status", "SF3B3_MTSS2.COG4_SF3B3")], id.vars = c("Group", "Species_name", "IL34_status"))
df_melt_SF3B3_MTSS2.MTSS2_VAC14 <- melt(df_filtered_no_gaps[,c("Group", "Species_name", "IL34_status", "SF3B3_MTSS2.MTSS2_VAC14")], id.vars = c("Group", "Species_name", "IL34_status"))
df_melt_COG4_SF3B3.MTSS2_VAC14 <- melt(df_filtered_no_gaps[,c("Group", "Species_name", "IL34_status", "COG4_SF3B3.MTSS2_VAC14")], id.vars = c("Group", "Species_name", "IL34_status"))

# Create boxplot

# Define pairwise comparisons with Squamata
pairwise_comparisons <- list(
  c("Squamata", "Aves"), 
  c("Squamata", "Mammalia"), 
  c("Squamata", "Amphibia"),
  c("Squamata", "Chondrichthyes"),
  c("Squamata", "Testudines")
)

# Create boxplot with significance annotations

p_df_melt_SF3B3_MTSS2.COG4_SF3B3 <- ggplot(df_melt_SF3B3_MTSS2.COG4_SF3B3, aes(x = factor(Group), y = value)) +
  geom_boxplot(aes(fill = factor(Group))) +
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = pairwise_comparisons, 
                     method.args = list(alternative = "less"), 
                     label = "p.signif") +  # Show significance as '*'
  labs(x = "Group", y = "Gene distance ratio", fill = "Group") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))
p_df_melt_SF3B3_MTSS2.MTSS2_VAC14 <- ggplot(df_melt_SF3B3_MTSS2.MTSS2_VAC14, aes(x = factor(Group), y = value)) +
  geom_boxplot(aes(fill = factor(Group))) +
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = pairwise_comparisons, 
                     method.args = list(alternative = "less"), 
                     label = "p.signif") +  # Show significance as '*'
  labs(x = "Group", y = "Gene distance ratio", fill = "Group") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))
p_df_melt_COG4_SF3B3.MTSS2_VAC14 <- ggplot(df_melt_COG4_SF3B3.MTSS2_VAC14, aes(x = factor(Group), y = value)) +
  geom_boxplot(aes(fill = factor(Group)),outlier.shape = NA) +
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = pairwise_comparisons, 
                     method.args = list(alternative = "two.sided"), 
                     label = "p.signif") +  # Show significance as '*'
  labs(x = "Group", y = "Gene distance ratio", fill = "Group") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Remove individual legends
T1p1 <- p_df_melt_SF3B3_MTSS2.COG4_SF3B3 + theme(legend.position = "none")
T1p2 <- p_df_melt_SF3B3_MTSS2.MTSS2_VAC14 + theme(legend.position = "none")
T1p3 <- p_df_melt_COG4_SF3B3.MTSS2_VAC14 + theme(legend.position = "none")

# Extract common legend
legend <- get_legend(
  p_df_melt_SF3B3_MTSS2.COG4_SF3B3 + 
    theme(legend.position = "bottom") + 
    guides(fill = guide_legend(nrow = 1))  # Force legend to be in one row
)
# Arrange plots with the common legend
T1_plot<- grid.arrange(
  arrangeGrob(T1p1, T1p2, T1p3, ncol = 3),  # Arrange plots in a column
  legend,  # Add common legend below
  ncol = 1,  # Stack vertically
  heights = c(3, 0.2)  # Adjust legend size
)


##IL34_gene_status
T1gene_p_df_melt_SF3B3_MTSS2.COG4_SF3B3<-ggplot(df_melt_SF3B3_MTSS2.COG4_SF3B3, aes(x = factor(IL34_status), y = value)) +
  geom_boxplot(aes(fill = factor(IL34_status))) + 
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "less"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Gene Status", y = "Gene distance ratio", fill = "Gene Status") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom")  # Move legend to bottom

T1gene_p_df_melt_SF3B3_MTSS2.MTSS2_VAC14<-ggplot(df_melt_SF3B3_MTSS2.MTSS2_VAC14, aes(x = factor(IL34_status), y = value)) +
  geom_boxplot(aes(fill = factor(IL34_status))) +
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "less"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Gene Status", y = "Gene distance ratio", fill = "Gene Status") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom")

T1gene_p_df_df_melt_COG4_SF3B3.MTSS2_VAC14<-ggplot(df_melt_COG4_SF3B3.MTSS2_VAC14, aes(x = factor(IL34_status), y = value)) +
  geom_boxplot(aes(fill = factor(IL34_status))) + 
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Gene Status", y = "Gene distance ratio", fill = "Gene Status") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom")

T1Gp1 <- T1gene_p_df_melt_SF3B3_MTSS2.COG4_SF3B3 + theme(legend.position = "none")
T1Gp2 <- T1gene_p_df_melt_SF3B3_MTSS2.MTSS2_VAC14 + theme(legend.position = "none")
T1Gp3 <- T1gene_p_df_df_melt_COG4_SF3B3.MTSS2_VAC14 + theme(legend.position = "none")

legend <- get_legend(
  T1gene_p_df_melt_SF3B3_MTSS2.COG4_SF3B3 + 
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

Main_fig_Gene_plot<- grid.arrange(
  arrangeGrob(T1Gp2, T1Gp3, ncol = 1),  # Arrange plots in a column
  legend,  # Add common legend below
  ncol = 1,  # Stack vertically
  heights = c(4, 0.2)  # Adjust legend size
)
ggsave(filename = "Main_fig_Gene_plot.png", plot = Main_fig_Gene_plot, width = 2.60, height = 5.20,dpi = 600) 

ggsave(filename = "Group_plot.png", plot = T1_plot, width = 16, height = 9,dpi = 600) 


ggsave(filename = "Gene_plot.png", plot = Gene_T1_plot, width = 16, height = 9,dpi = 600) 

################## Phylogenetic signal ####################
library(phylolm)
library(ape)

df <- read.table("final_output.sorted.tsv", header = TRUE, sep = "\t")
df_filtered_no_gaps <- df[df$Num_Gaps == 0, ]
df_filtered_no_gaps <- df_filtered_no_gaps[!(df_filtered_no_gaps$Group %in% c("Actinopteri", "Crocodilia")), ]
head(df_filtered_no_gaps,3)
df_filtered_no_gaps <- df_filtered_no_gaps[,c(3,4,6,9,11)]
df_filtered_no_gaps$IL34_status[df_filtered_no_gaps$IL34_status == "Intact"] <- 1
df_filtered_no_gaps$IL34_status[df_filtered_no_gaps$IL34_status == "Loss"] <- 0
#View(df_filtered_no_gaps)
write.table(df_filtered_no_gaps, "phylolm.df_filtered_no_gaps.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
colnames(df_filtered_no_gaps)
row_sums <- rowSums(df_filtered_no_gaps[, c("COG4_SF3B3", "SF3B3_MTSS2", "MTSS2_VAC14")])

# Compute percentages
df_filtered_no_gaps$COG4_SF3B3 <- (df_filtered_no_gaps$COG4_SF3B3 / row_sums) * 100
df_filtered_no_gaps$SF3B3_MTSS2 <- (df_filtered_no_gaps$SF3B3_MTSS2 / row_sums) * 100
df_filtered_no_gaps$MTSS2_VAC14 <- (df_filtered_no_gaps$MTSS2_VAC14 / row_sums) * 100
head(df_filtered_no_gaps,5)

tree=read.tree("phylolm.df_filtered_no_gaps.species.nwk")
rownames(df_filtered_no_gaps)=df_filtered_no_gaps[,1]
head(df_filtered_no_gaps)
colnames(df_filtered_no_gaps)
str(df_filtered_no_gaps)
df_filtered_no_gaps$IL34_status <- as.numeric(df_filtered_no_gaps$IL34_status)
View(df_filtered_no_gaps)

# Identify tips with zero branch length
zero_length_tips <- tree$tip.label[tree$edge.length[match(1:length(tree$tip.label), tree$edge[,2])] == 0]
tree <- drop.tip(tree, c("Cygnus_olor", "Cygnus_atratus"))
df_filtered_no_gaps <- df_filtered_no_gaps[!df_filtered_no_gaps$Species_name %in% c("Cygnus_olor", "Cygnus_atratus"), ]
head(df_filtered_no_gaps,4)
num_tips <- length(tree$tip.label)
num_internal_nodes <- Nnode(tree)
total_nodes <- num_tips + num_internal_nodes
total_nodes
# Extract row names and tree tip labels
data_species <- rownames(df_filtered_no_gaps)
tree_species <- tree$tip.label

# Check for mismatches
setdiff(data_species, tree_species)  # Species in df_filtered_no_gaps but not in tree
setdiff(tree_species, data_species)  # Species in tree but not in df_filtered_no_gaps
###############MPLE###############
fitwmple=phyloglm(IL34_status~SF3B3_MTSS2 ,df_filtered_no_gaps,tree, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4,start.beta=NULL, start.alpha=NULL,boot = 2000, full.matrix = TRUE)
cc1=coef(fitwmple)
summary(fitwmple)
mpleinter=round(cc1[1],2)
mplew=round(cc1[2],2)
#####################IG10#########
fitwig10=phyloglm(IL34_status~SF3B3_MTSS2 ,df_filtered_no_gaps,tree, method = c("logistic_IG10"), btol = 20, log.alpha.bound = 4,start.beta=NULL, start.alpha=NULL,boot = 2000, full.matrix = TRUE)
cc2=coef(fitwig10)
summary(fitwig10)
ig10inter=round(cc2[1],2)
ig10w=round(cc2[2],2)
####################################
t(table(df_filtered_no_gaps$IL34_status,df_filtered_no_gaps$SF3B3_MTSS2))->SM
png("Phylogenetic_signal_of_deletion_plot.png",height=8,width=8,units = "in",res=600)
# in shell script cut -f3 -d' ' SM.tsv|tr '\n' ','
plot(rep(c(35.1052877866167,38.7640449438202,44.9738702735936,45.5138950154389,46.9730870037717,47.3333772826159,50.4181380484654,54.5740615868735,55.3638723724981,55.5536756619575,56.030701754386,58.5469107551487,59.8462115528882,66.2702746430449,67.8047748263586,67.8168291431503,68.1178701409915,69.0062198943086,70.8300869044051,73.0416992983522,73.2545737960456,73.480136020569,74.5584449490595,75.6795797767564,75.7317044106674,75.9137577002053,76.2267365191111,76.4604316546763,76.552183657399,76.6393442622951,76.7794043412418,76.7806848138565,76.9735798701689,77.2313535459367,77.658403986797,77.8548933772938,77.9474422984666,77.9610920034394,78.1549173194082,78.2821824381927,78.9144675276062,78.9664252797893,79.3247269116187,79.4141517983957,79.4664550999157,79.4823676597076,79.5345237437116,79.540439412837,79.6280991735537,79.6500820120284,79.7240779968497,79.7298769197255,79.7636551761695,79.7678727114211,79.8047544526418,79.8862096992685,79.8949818654252,79.9522931800932,80.1029444321452,80.2219162549607,80.2820766106904,80.5724674446136,80.7939865868479,80.9348047154533,81.0829916584193,81.5627423254117,81.6604148777512,81.7024296568423,81.8439905495008,82.326945316372,82.7801527450886,84.0045362466887,84.8997667154324,85.8762958324636,87.0002920845098,87.163866314135,87.6477771715125,87.9813602127707,88.4059065674686,88.4097001120018,88.4183078837683,88.4435537742151,88.4990445663564,88.9087478899811,89.0573020625416,89.1456956361098,89.1567344744153,90.2993032477146,90.4434909982443,90.4748565677944,91.0873938310237,91.6202380010623,91.6855372783593,92.0738583555363,92.0987691702118,92.2025370657817,92.2400710916672,92.2513363379101,92.2879039605552,92.2994324282865,92.3577014478492,92.3740510697032,92.4557003659088,92.4671401484719,92.516758213523,92.6376888068069,92.7452166837065,92.9027327815389,92.9262487853652,92.9817991161711,93.0484816782755,93.0684828896811,93.0894085281981,93.1097872810958,93.1156380242363,93.1370129692084,93.1730241230031,93.2283687943262,93.2337606335019,93.2567249934709,93.2577876362234,93.2678638966052,93.2945650685292,93.3141270212048,93.3395968819272,93.3666119467989,93.3757223716073,93.376365030357,93.3803384367446,93.3860901166618,93.3886288185826,93.3958328988258,93.4529693559425,93.5453009293773,93.5573649532045,93.5791332778678,93.5946194292715,93.608090059149,93.6214669051878,93.6353425683265,93.6506636080616,93.6656642763513,93.6855096798266,93.6869926373419,93.6973898281106,93.7195424528902,93.719616514413,93.75,93.7782404536192,93.7855277477919,93.8291102492045,94.0275407333595,94.0482598825099,94.155987197595,94.1937568798538,94.2118432026689,94.2139547377611,94.3634617436773,94.3889881265105,94.4012936404846,94.4861122558312,94.5330711232261,94.5944919759024,94.9941145697097,95.1654368338406,95.1769027698239,95.1805735379191,95.1929158760278,95.1948008678286,95.2943015711717,95.3146244629358,95.5275913692564,95.6043158003942,95.8363341186236,96.0207876127774,96.1735495631404,96.2594038796024,96.2856025015392,97.6535860635573),2),c(rep(1,179),rep(0,179)),cex=c(as.vector(SM[,2]),as.vector(SM[,1]))/3,pch=19, xlab="",ylab='',xlim=c(0,100),axes=F)
axis(1)
axis(2, at = seq(0,1,by=1), las = 1)
box()
mtext("Loss",side=2,line=1.5,at=0)
mtext(expression(paste("Relative distance (in %): ", frac(italic('SF3B3-MTSS2'), italic('COG4-SF3B3') + italic('SF3B3-MTSS2') + italic('MTSS2-VAC14')), " X 100")), side = 1, line = 4, at = 50)
mtext(substitute(paste(italic('IL34'), " gene status")),side=2,line=2,at=0.5)
mtext("Intact",side=2,line=1.5,at=1)
curve(plogis(cc1[1]+cc1[2]*x),col="red",add=TRUE,lwd=2)
textleg=paste("logistic_MPLE\nIntercept = ",mpleinter,"\nSlope = ",mplew,"\np = 0.002175\nn = 179")
legend(85,0.8,textleg,xjust = 0.5,yjust = 0.5,x.intersp = 0.2,y.intersp = 0.8,adj = c(0, 0.5),bty='n', text.col = "red")
curve(plogis(cc2[1]+cc2[2]*x),col="blue",add=TRUE,lwd=2)
text2leg=paste("logistic_IG10\nIntercept = ",ig10inter,"\nSlope = ",ig10w,"\np = 0.001806\nn = 179")
legend(40,0.4,text2leg,xjust = 0.5,yjust = 0.5,x.intersp = 0.2,y.intersp = 0.8,adj = c(0, 0.5),bty='n', text.col = "blue")
dev.off()


