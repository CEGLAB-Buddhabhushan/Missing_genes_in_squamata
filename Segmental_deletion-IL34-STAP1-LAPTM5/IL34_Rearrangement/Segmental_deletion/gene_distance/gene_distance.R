setwd("/home/ceglab358/BUDDHA/BMCBIO-EoI/IL34_Rearrangement/Segmental_deletion/gene_distance/")
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr) 
library(patchwork)
library(gridExtra)
library(cowplot)
# Load data
df1 <- read.table("1T_All_CDS.final_output.sorted.tsv", header = TRUE, sep = "\t")
colnames(df1)
# Identify numeric columns (excluding "gene_status")
numeric_cols1 <- setdiff(names(df1)[sapply(df1, is.numeric)], c("Group", "Species_name", "IL34_status"))
df1 <- df1[!(df1$Group %in% c("Crocodylia", "Actinopterygii")), ]

head(df1,10)
df1$SF3B3.MTSS2_FCSK.SF3B3 <- df1$SF3B3.MTSS2/df1$FCSK.SF3B3
df1$SF3B3.MTSS2_MTSS2.HYDIN <- df1$SF3B3.MTSS2/df1$MTSS2.HYDIN
df1$FCSK.SF3B3_MTSS2.HYDIN <- df1$FCSK.SF3B3/df1$MTSS2.HYDIN


df_melt_SF3B3.MTSS2_FCSK.SF3B3 <- melt(df1[,c("Group", "Species_name", "IL34_status", "SF3B3.MTSS2_FCSK.SF3B3")], id.vars = c("Group", "Species_name", "IL34_status"))
df_melt_SF3B3.MTSS2_MTSS2.HYDIN <- melt(df1[,c("Group", "Species_name", "IL34_status", "SF3B3.MTSS2_MTSS2.HYDIN")], id.vars = c("Group", "Species_name", "IL34_status"))
df_melt_FCSK.SF3B3_MTSS2.HYDIN <- melt(df1[,c("Group", "Species_name", "IL34_status", "FCSK.SF3B3_MTSS2.HYDIN")], id.vars = c("Group", "Species_name", "IL34_status"))

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

p_df_melt_SF3B3.MTSS2_FCSK.SF3B3 <- ggplot(df_melt_SF3B3.MTSS2_FCSK.SF3B3, aes(x = factor(Group), y = value)) +
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
p_df_melt_SF3B3.MTSS2_MTSS2.HYDIN <- ggplot(df_melt_SF3B3.MTSS2_MTSS2.HYDIN, aes(x = factor(Group), y = value)) +
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
p_df_melt_FCSK.SF3B3_MTSS2.HYDIN <- ggplot(df_melt_FCSK.SF3B3_MTSS2.HYDIN, aes(x = factor(Group), y = value)) +
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
T1p1 <- p_df_melt_SF3B3.MTSS2_FCSK.SF3B3 + theme(legend.position = "none")
T1p2 <- p_df_melt_SF3B3.MTSS2_MTSS2.HYDIN + theme(legend.position = "none")
T1p3 <- p_df_melt_FCSK.SF3B3_MTSS2.HYDIN + theme(legend.position = "none")

# Extract common legend
legend <- get_legend(
  p_df_melt_SF3B3.MTSS2_FCSK.SF3B3 + 
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
T1gene_p_df_melt_SF3B3.MTSS2_FCSK.SF3B3<-ggplot(df_melt_SF3B3.MTSS2_FCSK.SF3B3, aes(x = factor(IL34_status), y = value)) +
  geom_boxplot(aes(fill = factor(IL34_status))) + 
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "less"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Gene Status", y = "Gene distance ratio", fill = "Gene Status") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom")  # Move legend to bottom

T1gene_p_df_melt_SF3B3.MTSS2_MTSS2.HYDIN<-ggplot(df_melt_SF3B3.MTSS2_MTSS2.HYDIN, aes(x = factor(IL34_status), y = value)) +
  geom_boxplot(aes(fill = factor(IL34_status))) +
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "less"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Gene Status", y = "Gene distance ratio", fill = "Gene Status") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom")

T1gene_p_df_melt_FCSK.SF3B3_MTSS2.HYDIN<-ggplot(df_melt_FCSK.SF3B3_MTSS2.HYDIN, aes(x = factor(IL34_status), y = value)) +
  geom_boxplot(aes(fill = factor(IL34_status))) + 
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Gene Status", y = "Gene distance ratio", fill = "Gene Status") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom")

T1Gp1 <- T1gene_p_df_melt_SF3B3.MTSS2_FCSK.SF3B3 + theme(legend.position = "none")
T1Gp2 <- T1gene_p_df_melt_SF3B3.MTSS2_MTSS2.HYDIN + theme(legend.position = "none")
T1Gp3 <- T1gene_p_df_melt_FCSK.SF3B3_MTSS2.HYDIN + theme(legend.position = "none")

legend <- get_legend(
  T1gene_p_df_melt_SF3B3.MTSS2_FCSK.SF3B3 + 
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

################# 2nd case
df2 <- read.table("2T_All_CDS.final_output.sorted.tsv", header = TRUE, sep = "\t")
colnames(df2)
# Identify numeric columns (excluding "gene_status")
numeric_cols2 <- setdiff(names(df2)[sapply(df2, is.numeric)], c("Group", "Species_name", "IL34_status"))
df2 <- df2[!(df2$Group %in% c("Crocodylia", "Actinopterygii")), ]

head(df2,10)
df2$SF3B3.MTSS2_COG4.SF3B3 <- df2$SF3B3.MTSS2/df2$COG4.SF3B3
df2$SF3B3.MTSS2_MTSS2.VAC14 <- df2$SF3B3.MTSS2/df2$MTSS2.VAC14
df2$COG4.SF3B3_MTSS2.VAC14 <- df2$COG4.SF3B3/df2$MTSS2.VAC14

df_melt_SF3B3.MTSS2_COG4.SF3B3 <- melt(df2[,c("Group", "Species_name", "IL34_status", "SF3B3.MTSS2_COG4.SF3B3")], id.vars = c("Group", "Species_name", "IL34_status"))
df_melt_SF3B3.MTSS2_MTSS2.VAC14 <- melt(df2[,c("Group", "Species_name", "IL34_status", "SF3B3.MTSS2_MTSS2.VAC14")], id.vars = c("Group", "Species_name", "IL34_status"))
df_melt_COG4.SF3B3_MTSS2.VAC14 <- melt(df2[,c("Group", "Species_name", "IL34_status", "COG4.SF3B3_MTSS2.VAC14")], id.vars = c("Group", "Species_name", "IL34_status"))

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

p_df_melt_SF3B3.MTSS2_COG4.SF3B3 <- ggplot(df_melt_SF3B3.MTSS2_COG4.SF3B3, aes(x = factor(Group), y = value)) +
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
p_df_melt_SF3B3.MTSS2_MTSS2.VAC14 <- ggplot(df_melt_SF3B3.MTSS2_MTSS2.VAC14, aes(x = factor(Group), y = value)) +
  geom_boxplot(aes(fill = factor(Group)),outlier.shape = NA) +
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = pairwise_comparisons, 
                     method.args = list(alternative = "less"), 
                     label = "p.signif") +  # Show significance as '*'
  labs(x = "Group", y = "Gene distance ratio", fill = "Group") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))
p_df_melt_COG4.SF3B3_MTSS2.VAC14 <- ggplot(df_melt_COG4.SF3B3_MTSS2.VAC14, aes(x = factor(Group), y = value)) +
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
T2p1 <- p_df_melt_SF3B3.MTSS2_COG4.SF3B3 + theme(legend.position = "none")
T2p2 <- p_df_melt_SF3B3.MTSS2_MTSS2.VAC14 + theme(legend.position = "none")
T2p3 <- p_df_melt_COG4.SF3B3_MTSS2.VAC14 + theme(legend.position = "none")

# Extract common legend
legend <- get_legend(
  p_df_melt_SF3B3.MTSS2_COG4.SF3B3 + 
    theme(legend.position = "bottom") + 
    guides(fill = guide_legend(nrow = 1))  # Force legend to be in one row
)
# Arrange plots with the common legend
T2_plot<- grid.arrange(
  arrangeGrob(T2p1, T2p2, T2p3, ncol = 3),  # Arrange plots in a column
  legend,  # Add common legend below
  ncol = 1,  # Stack vertically
  heights = c(3, 0.2)  # Adjust legend size
)


##IL34_gene_status
T2_gene_p_df_melt_SF3B3.MTSS2_COG4.SF3B3<-ggplot(df_melt_SF3B3.MTSS2_COG4.SF3B3, aes(x = factor(IL34_status), y = value)) +
  geom_boxplot(aes(fill = factor(IL34_status))) + 
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "less"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Gene Status", y = "Gene distance ratio", fill = "Gene Status") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom")  # Move legend to bottom

T2_gene_p_df_melt_SF3B3.MTSS2_MTSS2.VAC14<-ggplot(df_melt_SF3B3.MTSS2_MTSS2.VAC14, aes(x = factor(IL34_status), y = value)) +
  geom_boxplot(aes(fill = factor(IL34_status))) +
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "less"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Gene Status", y = "Gene distance ratio", fill = "Gene Status") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom")

T2_gene_p_df_melt_COG4.SF3B3_MTSS2.VAC14<-ggplot(df_melt_COG4.SF3B3_MTSS2.VAC14, aes(x = factor(IL34_status), y = value)) +
  geom_boxplot(aes(fill = factor(IL34_status))) + 
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Gene Status", y = "Gene distance ratio", fill = "Gene Status") +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "bottom")

T2Gp1 <- T2_gene_p_df_melt_SF3B3.MTSS2_COG4.SF3B3 + theme(legend.position = "none")
T2Gp2 <- T2_gene_p_df_melt_SF3B3.MTSS2_MTSS2.VAC14 + theme(legend.position = "none")
T2Gp3 <- T2_gene_p_df_melt_COG4.SF3B3_MTSS2.VAC14 + theme(legend.position = "none")

legend <- get_legend(
  T2_gene_p_df_melt_SF3B3.MTSS2_COG4.SF3B3 + 
    theme(legend.position = "bottom") + 
    guides(fill = guide_legend(nrow = 1))  # Force legend to be in one row
)
# Arrange plots with the common legend
Gene_2T_plot<- grid.arrange(
  arrangeGrob(T2Gp1, T2Gp2, T2Gp3, ncol = 3),  # Arrange plots in a column
  legend,  # Add common legend below
  ncol = 1,  # Stack vertically
  heights = c(3, 0.2)  # Adjust legend size
)

Group_plot<- grid.arrange(
  arrangeGrob(T1p1, T1p2, T1p3,T2p1, T2p2, T2p3, ncol = 3),  # Arrange plots in a column
  legend,  # Add common legend below
  ncol = 1,  # Stack vertically
  heights = c(3, 0.2)  # Adjust legend size
)
ggsave(filename = "Group_plot.png", plot = Group_plot, width = 16, height = 9,dpi = 600) 

Gene_plot<- grid.arrange(
  arrangeGrob(T1Gp1, T1Gp2, T1Gp3, T2Gp1, T2Gp2, T2Gp3, ncol = 3),  # Arrange plots in a column
  legend,  # Add common legend below
  ncol = 1,  # Stack vertically
  heights = c(3, 0.2)  # Adjust legend size
)
ggsave(filename = "Gene_plot.png", plot = Gene_plot, width = 16, height = 9,dpi = 600) 

