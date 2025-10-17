##### Figure 4 scripts
#HEATMAP
library(tidyheatmaps)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
df_heatmap <- read.table("FS_Exonic.Final.Gallus_gallus.53_genes.focal.tsv", header = T, sep = '\t')
head(df_heatmap)
unique(df_heatmap$Class)
summary(df_heatmap$Normalized_value)
#C11orf74==IFTAP
#HRASLS==PLAAT1
#C10orf128==TMEM273
df_heatmap <- df_heatmap %>%
  mutate(Gene_name = recode(Gene_name,
                            "C11orf74" = "IFTAP",
                            "HRASLS"   = "PLAAT1",
                            "C10orf128" = "TMEM273"))
df_heatmap <- df_heatmap %>%
  arrange(Gene_name, Class, Order)
df_heatmap$Gene_name <- factor(df_heatmap$Gene_name, levels = unique(df_heatmap$Gene_name))
head(df_heatmap)
plot_A <-tidyheatmap(df = df_heatmap,
                     rows = Species_name,
                     columns = Gene_name,
                     values = Normalized_value,
                     annotation_row = c(Order, Class),
                     gaps_row = Order,
                     colors = c("yellow","purple"),
                     color_legend_min = 0,
                     color_legend_max = 1,
                     cellheight = 8,
                     cellwidth = 8 )
png(filename="Fig1.Focal_genes_exonic_chicken_heatmap.png",
    width = 12,
    height = 5,
    res = 900, units = "in")
print(plot_A)
dev.off()

# BOXPLOT ACROSS ORDERS (Vertical comparison)
df_bxp_V <- read.table("FS_Exonic.Final.Gallus_gallus.53_genes.focal.tsv", header = TRUE, sep = "\t")
shapiro.test(df_bxp_V$Normalized_value)
df_bxp_V$Class <- as.factor(as.character(df_bxp_V$Class))
df_bxp_V$Normalized_value <- as.numeric(as.character(df_bxp_V$Normalized_value))

df_bxp_V <- df_bxp_V %>%
  mutate(Group = case_when(
    Order %in% c("Anseriformes", "Accipitriformes", "Galliformes", 
                 "Passeriformes", "Psittaciformes", "Struthioniformes") ~ "Birds",
    TRUE ~ Order
  ))
df_bxp_V$Group <- factor(df_bxp_V$Group)
df_bxp_V_comparisons <- list( c("Birds", "Squamata"), c("Testudines", "Squamata"), c("Coelacanthiformes", "Squamata"), c("Crocodylia", "Squamata"), c("Sphenodontia", "Squamata"))
plot_B <- ggboxplot(df_bxp_V, x = "Group", y = "Normalized_value",
                    fill  = "Group")+ 
  stat_compare_means(comparisons = df_bxp_V_comparisons, method = "wilcox.test", 
                     method.args = list(alternative = "greater"), label = "p.signif") +
  stat_compare_means(label.y = 1.5) +
  labs(x = "Group", y = "Aligned region (Normalized with CDS length)", fill = "Group") +
  theme_classic() + 
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal") +
  guides(fill = guide_legend(nrow = 1)) +
 # scale_fill_manual(values=c("Birds"="#520120","Testudines"="#08403E","Squamata"="#706513","Coelacanthiformes"="#B57114","Crocodylia"="#962B09","Sphenodontia"="#8c564b")) 
  scale_fill_manual(values=c("Birds"="#1f77b4","Testudines"="#ff7f0e","Squamata"="#2ca02c","Coelacanthiformes"="#d62728","Crocodylia"="#9467bd","Sphenodontia"="#8c564b")) 

png("Fig1.Comparison_across_orders.bxp.png", width = 8, height = 4, units = "in", res = 900)
print(plot_B)
dev.off()
png("Supp.Comparison_across_orders.bxp.png", width = 8, height = 64, units = "in", res = 900)
print(plot_B + facet_wrap(~ Gene_name, ncol = 2))
dev.off()

# BOXPLOT ACROSS FLANKING GENES OF SAQUMATES

df_bxp_flank<- read.table("Pairwise_bxp.tsv", header = T, sep = '\t')
head(df_bxp_flank)
# Convert to long format
df_bxp_flank_long <- df_bxp_flank %>%
  pivot_longer(cols = c(Left, Focal, Right),
               names_to = "Gene_position",
               values_to = "Normalized_value")
head(df_bxp_flank_long)

df_bxp_flank_long %>%
  group_by(Gene_position) %>%
  summarise(
    shapiro_p = shapiro.test(Normalized_value)$p.value
  )
unique(df_long$Gene_position)
df_bxp_flank_long <- df_bxp_flank_long %>%
  mutate(Group = case_when(
    order %in% c("Anseriformes", "Accipitriformes", "Galliformes", 
                 "Passeriformes", "Psittaciformes", "Struthioniformes") ~ "Birds",
    TRUE ~ order
  ))
head(df_bxp_flank_long)

df_bxp_flank_long_comparisons <- list( c("Left", "Focal"), c("Right", "Focal"))
plot_C <- ggboxplot(df_bxp_flank_long, x = "Gene_position", y = "Normalized_value",
                    fill  = "Gene_position")+ 
  stat_compare_means(comparisons = df_bxp_flank_long_comparisons, method = "wilcox.test", 
                     method.args = list(alternative = "greater"), label = "p.signif") +
  stat_compare_means(label.y = 1.2) +
  labs(x = "Gene_position", y = "Aligned region (Normalized with CDS length)", fill = "Gene_position") +
  theme_classic() + 
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal") +
  scale_fill_manual(values=c("Left"="#005E54","Right"="#C2BB00","Focal"="#E1523D"))+
  #  scale_fill_manual(values=c("Left"="#1f77b4","Right"="#ff7f0e","Focal"="#2ca02c"))+
  
  facet_wrap(~Group, ncol=1)



png("Fig1.Comparison_of_flanking_region_cross_orders.bxp.png", width = 7.2, height = 16.2, units = "in", res = 900)
print(plot_C)
dev.off()

png("Supp.Comparison_of_flanking_region_cross_orders.bxp.png", width = 8, height = 48, units = "in", res = 900)
print(ggboxplot(df_bxp_flank_long, x = "Gene_position", y = "Normalized_value",
                fill  = "Gene_position")+ 
        stat_compare_means(comparisons = df_bxp_flank_long_comparisons, method = "wilcox.test", 
                           method.args = list(alternative = "greater"), label = "p.signif") +
        stat_compare_means(label.y = 1.2) +
        labs(x = "Gene_position", y = "Aligned region (Normalized with CDS length)", fill = "Gene_position") +
        theme_classic() + 
        theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal") +
        scale_fill_manual(values=c("Left"="#005E54","Right"="#C2BB00","Focal"="#E1523D"))+ facet_wrap(~ focal_gene, ncol = 1))
dev.off()
# BOXPLOT WITHIN SQUAMATES
df_heatmap_Squamata <- df_heatmap[df_heatmap$Order == "Squamata", ]
head(df_heatmap_Squamata)
#shapiro.test(df_heatmap_Squamata$Normalized_value)
#hist(df_focal_comp$Normalized_value, breaks=30,
#     main="Histogram of Normalized values", col="lightblue")
#qqnorm(df_focal_comp$Normalized_value); qqline(df_focal_comp$Normalized_value, col="red")
Mannual_chains <- c("Anolis_sagrei", 
                    "Eublepharis_macularius", 
                    "Euleptes_europaea", 
                    "Podarcis_raffonei", 
                    "pogona_vitticeps_zw", 
                    "pogona_vitticeps_zz", 
                    "Rhineura_floridana")

# Add grouping column
df_heatmap_Squamata$Group <- ifelse(df_heatmap_Squamata$Species_name %in% Mannual_chains, 
                                    "Mannual_chains", "ENSEMBL_chains")
wilcox.test(Normalized_value ~ Group, data = df_heatmap_Squamata)
head(df_heatmap_Squamata)
df_heatmap_Squamata_comparisons <- list( c("Mannual_chains", "ENSEMBL_chains"))
plot_D <- ggboxplot(df_heatmap_Squamata, x = "Group", y = "Normalized_value",
                    fill  = "Group")+ 
  stat_compare_means(comparisons = df_heatmap_Squamata_comparisons, method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"), 
                     label = "p.signif") +
  stat_compare_means(label.y = 1.4) +
  labs(x = "Group", y = "Aligned region (Normalized with CDS length)", fill = "Group") +
  theme_classic() + 
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal") +
  scale_fill_manual(values=c("Mannual_chains"="#1f77b4","ENSEMBL_chains"="#ff7f0e"))


png("Fig1.Comparison_of_squamates_chains.bxp.png", width = 4, height = 4, units = "in", res = 900)
print(plot_D)
dev.off()














#heatmap
library(tidyheatmaps)
library(dplyr)

df_focal<- read.table("./plotting/FS_Exonic.Final.Gallus_gallus.53_genes.focal.tsv", header = T, sep = '\t')
head(df_focal)
unique(df$Class)
summary(df_focal$Normalized_value)
#C11orf74==IFTAP
#HRASLS==PLAAT1
#C10orf128==TMEM273
df_focal <- df_focal %>%
  mutate(Gene_name = recode(Gene_name,
                            "C11orf74" = "IFTAP",
                            "HRASLS"   = "PLAAT1",
                            "C10orf128" = "TMEM273"))
df_focal <- df_focal %>%
  arrange(Gene_name, Class, Order)
df_focal$Gene_name <- factor(df_focal$Gene_name, levels = unique(df_focal$Gene_name))

View(df_focal)
head(df_focal)
Focal_genes_exonic_chicken_heatmap<-tidyheatmap(df = df_focal,
                                                rows = Species_name,
                                                columns = Gene_name,
                                                values = Normalized_value,
                                                annotation_row = c(Order, Class),
                                                gaps_row = Order,
                                                colors = c("yellow","purple"),
                                                color_legend_min = 0,
                                                color_legend_max = 1,
                                                cellheight = 8,
                                                cellwidth = 8,
)
png(filename="Focal_genes_exonic_chicken_heatmap.png",
    width = 9,
    height = 5,
    res = 900, units = "in")
print(Focal_genes_exonic_chicken_heatmap)
dev.off()



library(dplyr)
library(tidyr)
library(pheatmap)

# Wide format: Species as rows, Genes as columns
heatmap_matrix <- df_focal %>%
  select(Species_name, Gene_name, Normalized_value) %>%
  pivot_wider(names_from = Gene_name, values_from = Normalized_value) %>%
  as.data.frame()

# Make Species_name rownames
rownames(heatmap_matrix) <- heatmap_matrix$Species_name
heatmap_matrix <- heatmap_matrix[ , -1]   # remove Species_name col

# Ensure it's numeric
heatmap_matrix <- as.matrix(heatmap_matrix)

# Heatmap
png(filename="Pheatmap.Focal_genes_exonic_chicken_heatmap.png",
    width = 12,
    height = 5,
    res = 900, units = "in")
pheatmap(
  heatmap_matrix,
  color = colorRampPalette(c("yellow", "purple"))(10),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  cellwidth = 8,
  cellheight = 8,
  show_rownames = TRUE,
  show_colnames = TRUE
)
dev.off()

png(filename="Color.Focal_genes_exonic_chicken_heatmap.png",
    width = 12,
    height = 5,
    res = 900, units = "in")
tidyheatmap(df = df_focal,
            rows = Species_name,
            columns = Gene_name,
            values = Normalized_value,
            annotation_row = c(Order, Class),
            gaps_row = Order,
            colors = c("#145afc","#ffffff","#ee4445"),
            color_legend_min = 0,
            color_legend_max = 1,
            cellheight = 8,
            cellwidth = 8,
            
            
)
dev.off()

### compare ENSEMBL chain file with new chain files with  BUSCO > 97
library(ggpubr)
df_focal_comp <- df_focal[df_focal$Order == "Squamata", ]
head(df_focal_comp)
#shapiro.test(df_focal_comp$Normalized_value)
#hist(df_focal_comp$Normalized_value, breaks=30,
#     main="Histogram of Normalized values", col="lightblue")
#qqnorm(df_focal_comp$Normalized_value); qqline(df_focal_comp$Normalized_value, col="red")
Mannual_chains <- c("Anolis_sagrei", 
                    "Eublepharis_macularius", 
                    "Euleptes_europaea", 
                    "Podarcis_raffonei", 
                    "pogona_vitticeps_zw", 
                    "pogona_vitticeps_zz", 
                    "Rhineura_floridana")

# Add grouping column
df_focal_comp$Group <- ifelse(df_focal_comp$Species_name %in% Mannual_chains, 
                              "Mannual_chains", "ENSEMBL_chains")
wilcox.test(Normalized_value ~ Group, data = df_focal_comp)

library(ggplot2)
library(ggpubr)

# Boxplot with Wilcoxon p-value
png(filename="plot_ENS_MAN_boxplot.png",
    width = 7,
    height = 7,
    res = 900, units = "in")
ggplot(df_focal_comp, aes(x = Group, y = Normalized_value, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("Mannual_chain" = "skyblue", "ENSEMBL_chains" = "lightgreen")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom") +
  labs(y = "Normalized Value",
       x = "")
dev.off()

### Boxplot of left, focal, right gene combination
library(ggplot2)
library(tidyr)
library(dplyr)
df_bxp_flank<- read.table("./plotting/Pairwise_bxp.tsv", header = T, sep = '\t')
head(df_bxp_flank)
# Convert to long format
df_long <- df_bxp_flank %>%
  pivot_longer(cols = c(Left, Focal, Right),
               names_to = "Gene_position",
               values_to = "value")
head(df_long)
pairwise.wilcox.test(df_long$value, df_long$Gene_position, p.adjust.method = "BH")

df_long %>%
  group_by(Gene_position) %>%
  summarise(
    shapiro_p = shapiro.test(value)$p.value
  )
unique(df_long$Gene_position)
my_comparisons <- list( c("Left", "Focal"), c("Right", "Focal"))
#across orders
png("Comparison_of_flanking_region_cross_orders.bxp.png", width = 16, height = 9, units = "in", res = 900)
ggboxplot(df_long, x = "Gene_position", y = "value",
          color = "Gene_position")+ 
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 1.3) +
  facet_wrap(~order)
dev.off()

#across genes of squamates
df_long_squamates <- df_long[df_long$order=="Squamata", ]
head(df_long_squamates)
unique(df_long_squamates$order)
png("Comparison_of_flanking_region_in_squmates_acorss_genes.bxp.png", width = 32, height = 18, units = "in", res = 900)
ggboxplot(df_long_squamates, x = "Gene_position", y = "value",
          color = "Gene_position")+ 
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 1.3) +
  facet_wrap(~focal_gene)
dev.off()





########## Pairwise comparison across orders#################
# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# Read your TSV table
df <- read.delim("./plotting/FS_Exonic.Final.Gallus_gallus.53_genes.focal.tsv", header = TRUE, sep = "\t")
shapiro.test(df$Normalized_value)

df$Class <- as.factor(as.character(df$Class))
df$Normalized_value <- as.numeric(as.character(df$Normalized_value))
pairwise_comparisons <- combn(levels(df$Class), 2, simplify = FALSE)



# Plot
png("Comparison_of_alined_region.bxp.png", width = 9, height = 8, units = "in", res = 900)  

ggplot(df, aes(x = Class, y = Normalized_value, fill = Class)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", method = "wilcox.test", comparisons = pairwise_comparisons) +
  theme_bw()+
  theme(legend.position = "bottom")
dev.off()


### FOR MAIN FIGURE
head(df)
unique(df$Order)
df_main <- df %>%
  mutate(Group = case_when(
    Order %in% c("Anseriformes", "Accipitriformes", "Galliformes", 
                 "Passeriformes", "Psittaciformes", "Struthioniformes") ~ "Birds",
    TRUE ~ Order
  ))
head(df_main)
df_main$Group <- factor(df_main$Group)
unique(df_main$Group)

pairwise_comparisons1 <- combn(levels(df_main$Group), 2, simplify = FALSE)

df_main_comparisons <- list( c("Birds", "Squamata"), c("Testudines", "Squamata"), c("Coelacanthiformes", "Squamata"), c("Crocodylia", "Squamata"), c("Sphenodontia", "Squamata"))
png("Comparison_of_flanking_region_cross_orders.bxp.png", width = 16, height = 9, units = "in", res = 900)
ggboxplot(df_main, x = "Group", y = "Normalized_value",
          fill  = "Group")+ 
  stat_compare_means(comparisons = df_main_comparisons)+ 
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "greater"), 
                     label = "p.signif") +
  stat_compare_means(label.y = 1.4) +
  labs(x = "Group", y = "Aligned region (Normalized with gene length)", fill = "Group") +
  theme_classic() + 
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal") +
  scale_fill_manual(values=c("Birds"="#1f77b4","Testudines"="#ff7f0e","Squamata"="#2ca02c","Coelacanthiformes"="#d62728","Crocodylia"="#9467bd","Sphenodontia"="#8c564b"))
dev.off()



png("Comparison_of_alined_region_acorss_clades.bxp.png", width = 9, height = 8, units = "in", res = 900)  
ggplot(df_main, aes(x = Group, y = Normalized_value, fill = Group)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", method = "wilcox.test", comparisons = pairwise_comparisons1) +
  theme_bw()+
  theme(legend.position = "bottom")
dev.off()


###bxplot with alternative
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr) 
library(patchwork)
library(gridExtra)
library(cowplot)
# Create boxplot
# Create boxplot with significance annotations
P1_ST<-ggplot(df_main[df_main$Group== c("Squamata", "Testudines"),], aes(x = Group,  y = Normalized_value)) +
  geom_boxplot(aes(fill = factor(Group))) + 
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "greater"), 
                     label = "p.signif") +   # Show significance as '*'
  labs(x = "Group", y = "Aligned region (Normalized with gene length)", fill = "Gene Status") +
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


###############
orders <- unique(df$Order)
results <- data.frame(Group1 = character(), Group2 = character(),
                      W_statistic = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

for (i in 1:(length(orders) - 1)) {
  for (j in (i + 1):length(orders)) {
    group1 <- orders[i]
    group2 <- orders[j]
    
    x <- df$Normalized_value[df$Order == group1]
    y <- df$Normalized_value[df$Order == group2]
    
    test <- wilcox.test(x, y, alternative = "greater")  # group1 > group2
    
    results <- rbind(results, data.frame(
      Group1 = group1,
      Group2 = group2,
      W_statistic = test$statistic,
      p_value = test$p.value
    ))
  }
}

# multiple testing
results$p_adj <- p.adjust(results$p_value, method = "BH")

print(results)
write.csv(results, "squamata_vs_others_wilcox_results.csv", row.names = FALSE)




##### for focal and syntenic genes heatmap
file_list <- list.files(path = "./plotting/", 
                        pattern = "\\.FS_Exonic\\.Final\\.Gallus_gallus\\.21_genes\\.tsv$", 
                        full.names = TRUE)

# Loop over files
for(file in file_list) {
  
  df_heatmap <- read.table(file, header = TRUE, sep = "\t")
  
  # Preserve gene order
  df_heatmap$Gene_name <- factor(df_heatmap$Gene_name, levels = unique(df_heatmap$Gene_name))
  #C11orf74==IFTAP
  #HRASLS==PLAAT1
  #C10orf128==TMEM273
  
  df_heatmap <- df_heatmap %>%
    mutate(Gene_name = recode(Gene_name,
                              "C11orf74" = "IFTAP",
                              "HRASLS"   = "PLAAT1",
                              "C10orf128" = "TMEM273"))
  df_heatmap <- df_heatmap[df_heatmap$Order == "Squamata", ]
  
  # Generate alluvial plot
  p <-tidyheatmap(df = df_heatmap,
                  rows = Species_name,
                  columns = Gene_name,
                  values = Normalized_value,
                  colors = c("yellow","purple"),
                  color_legend_min = 0,
                  color_legend_max = 1,
                  cellheight = 8,
                  cellwidth = 8,
  )
  
  print(p)
  
  # Optional: save to PDF
  ggsave(filename = paste0(tools::file_path_sans_ext(basename(file)), "_heatmap.pdf"),
         plot = p,
         width = 12,
         height = 8)
}



#### aluvial plot
library(dplyr)
library(ggplot2)
library(ggalluvial)

# List all files matching pattern
file_list <- list.files(path = "./plotting/", 
                        pattern = "\\.FS_Exonic\\.Final\\.Gallus_gallus\\.21_genes\\.tsv$", 
                        full.names = TRUE)

# Loop over files
for(file in file_list) {
  
  df_alluvial <- read.table(file, header = TRUE, sep = "\t")
  
  # Preserve gene order
  df_alluvial$Gene_name <- factor(df_alluvial$Gene_name, levels = unique(df_alluvial$Gene_name))
  
  
  # Generate alluvial plot
  p <- ggplot(df_alluvial,
              aes(x = Gene_name, 
                  stratum = Class, 
                  alluvium = Species_name, 
                  y = Normalized_value,
                  fill = Class, 
                  label = Class)) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "darkgray") +
    geom_stratum() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(
        angle = 45, 
        hjust = 1, 
        face = "italic"  # italicize gene names
      )
    ) +
    facet_wrap(~ Class, nrow = 4)
  
  # Print the plot
  print(p)
  
  # Optional: save to PDF
  ggsave(filename = paste0(tools::file_path_sans_ext(basename(file)), "_alluvial.pdf"),
         plot = p,
         width = 12,
         height = 8)
}



