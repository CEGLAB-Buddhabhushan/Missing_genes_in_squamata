setwd("../sorted/")
getwd()
library(ape)
library(ggtree)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tidyverse)

df<- read.table("FEL_MEME_RELAX_aBSREL_BUSTED.out.tsv", sep = '\t', header = T)
tree <- read.tree("Species.nwk")

head(df,3)
df$relax_pvalue_fdr   <- p.adjust(df$relax_pvalue, method = "fdr")
df$absrel_pvalue_fdr  <- p.adjust(df$absrel_pvalue, method = "fdr")
df$busted_pvalue_fdr  <- p.adjust(df$busted_pvalue, method = "fdr")
df$RELAX_status <- "NS"  
df$RELAX_status[df$relax_pvalue_fdr < 0.05 & df$relax_kvalue < 1] <- "RELAXATION"
df$RELAX_status[df$relax_pvalue_fdr < 0.05 & df$relax_kvalue > 1] <- "INTENSIFICATION"
df$aBSREL_status <- "NS"  
df$aBSREL_status[df$absrel_pvalue_fdr < 0.05 ] <- "Positive Selection"
df$BUSTED_status <- "NS"  
df$BUSTED_status[df$busted_pvalue_fdr < 0.05 ] <- "Positive Selection"
write.csv(df, "selection_results_processed.csv", row.names = FALSE)
head(df,5)
colnames(df)

species_order <- c("Danio_rerio", "Oryzias_melastigma", "Poecilia_formosa", "Poecilia_reticulata",
                   "Latimeria_chalumnae", "Xenopus_tropicalis", "Rana_temporaria", "Hyla_sarda", "Bufo_bufo",
                   "Ornithorhynchus_anatinus", "Phascolarctos_cinereus", "Loxodonta_africana", "Homo_sapiens",
                   "Mus_musculus", "Rhinolophus_ferrumequinum", "Orcinus_orca", "Equus_caballus", "Canis_lupus",
                   "Sphenodon_punctatus", "Rhineura_floridana", "Thamnophis_elegans", "Varanus_komodoensis",
                   "Pogona_vitticeps", "Anolis_carolinensis", "Dermochelys_coriacea", "Caretta_caretta", "Alligator_mississippiensis",
                   "Gallus_gallus", "Cygnus_olor", "Calypte_anna", "Taeniopygia_guttata","Rissa_tridactyla" ) 


df_long <- df %>%
  select(Species_name, FEL_positive_sites, MEME_Episodic_div_site) %>%
  pivot_longer(cols = c(FEL_positive_sites, MEME_Episodic_div_site),
               names_to = "FEL_type", values_to = "Sites") %>%
  mutate(
    Sites = ifelse(FEL_type == "MEME_Episodic_div_site", -Sites, Sites),  # Make MEME negative
    Species_name = factor(Species_name, levels = rev(species_order))          # Align with tree
  )

# Plot
plot_dot_chart <- ggdotchart(df_long, x = "Species_name", y = "Sites",
           color = "FEL_type",
           palette = c("#00AFBB", "#FC4E07"),
           sorting = "none",
           add = "segments",
           add.params = list(color = "lightgray", size = 1),
           group = "FEL_type",
           dot.size = 5,
           label = abs(df_long$Sites),  # Show positive values as labels
           font.label = list(color = "white", size = 8, vjust = 0.5),
           ggtheme = theme_pubr(base_size = 14)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray") +
  labs(y = "Number of Selection Sites (+FEL / -MEME)",
       x = "Species") +
  coord_flip()+
  theme_classic() +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )


df1<-read.table("HYPHY_FEL_merged_data.tsv", header = T)
df2<- read.table("HYPHY_MEME_merged_data.tsv", header = T)
df2_1 <- df2[, c(1:6)]
df2_2 <- df2[, c(7:8)]
df2_1$Selection <- "Episodic diversifying selection" 
df2_final <- cbind(df2_1, df2_2)[,c(1,6,7,9)]
df1_final <-df1[,c(1,5,6,8)]
head(df2_final)
head(df1_final)
data<- rbind(df2_final, df1_final)
data$Species <- factor(data$Species, levels = species_order)
plot_sites1 <- ggplot(data, aes(x = Codon, colour = factor(Selection), shape = factor(Selection))) +
  geom_point(aes(y = as.numeric(factor(Selection))), size = 2) +
  geom_segment(aes(xend = Codon, y = 0, yend = as.numeric(factor(Selection)))) +
  facet_grid(vars(Species)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  ) +
  scale_color_manual(values = c(
    "Positive" = "#00AFBB",
    "Negative" = "lightgrey",
    "Episodic diversifying selection" = "#FC4E07"
  )) +
  scale_shape_manual(values = c(
    "Positive" = 17,
    "Negative" = 21,
    "Episodic diversifying selection" = 15
  ))

plot_sites2 <- ggplot(data, aes(x = Codon, colour = factor(Selection), shape = factor(Selection))) +
  geom_point(aes(y = as.numeric(factor(Selection))), size = 2) +
  geom_segment(aes(xend = Codon, y = 0, yend = as.numeric(factor(Selection)))) +
  facet_grid(vars(Species)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    strip.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  ) +
  scale_color_manual(values = c(
    "Positive" = "#00AFBB",
    "Negative" = "lightgrey",
    "Episodic diversifying selection" = "#FC4E07"
  )) +
  scale_shape_manual(values = c(
    "Positive" = 17,
    "Negative" = 21,
    "Episodic diversifying selection" = 15
  ))
#####
library(grid)
png("combined_plot.png", width = 12, height = 9, units = "in", res = 900)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(plot_sites2, vp = define_region(row = 1, col = 1:4))   # Span over two columns
print(plot_dot_chart, vp = define_region(row = 1, col = 5))
dev.off()

library(grid)
png("combined_plot_with_species_name.png", width = 12, height = 9, units = "in", res = 900)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(plot_sites1, vp = define_region(row = 1, col = 1:4))   # Span over two columns
print(plot_dot_chart, vp = define_region(row = 1, col = 5))
dev.off()
