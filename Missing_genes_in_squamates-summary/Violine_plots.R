setwd("../../")
getwd()
df <- read.table("./Squamata_genone_blast/GENE_wise/Gene_species_Max_coverage.tsv", header = T, sep = '\t')
df <- na.omit(df)
library(ggplot2)
library(dplyr)
library(vioplot)
library(RColorBrewer)   
# Summarize the data
df_summary <- df %>%
  group_by(Gene_name) %>%
  summarise(
    median_val = median(Max_cov, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
df$Gene_name <- factor(df$Gene_name, levels = df_summary$Gene_name)
n_colors <- nlevels(df$Gene_name)  # Replace 'group_column' with your actual column representing groups
palette <- colorRampPalette(brewer.pal(9, "Set1"))(n_colors)

png("Max_coverage.png", width = 16, height = 9, units = 'in', res = 900)
vioplot(Max_cov ~ Gene_name, data = df,
        names = rep("", nrow(df_summary)),
        xlab = "", ylab = "Max. Coverage", col = palette,colMed = "black", ylim = c(-5, 100))

x_pos <- 1:nrow(df_summary)
axis(1, at = x_pos, labels = FALSE)
text(x = x_pos, y = par("usr")[3] - 2, labels = df_summary$Gene_name,
     srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
text(x = x_pos,
     y = 96,  
     labels = round(df_summary$median_val, 0),
     cex = 1, col = "red")+
text(x = x_pos,
     y = 99, 
     labels = df_summary$n,
     cex = 0.7)
abline(h = 70, lty = 2, col = "red")
title(xlab = "Gene Names")
dev.off()
##main figure plot
n_genes <- nrow(df_summary)

# Set margins and plot area tightly
par(mar = c(0.5, 0.5, 0, 0)) 
png("Main_fig.Max_coverage.png", width = 16, height = 3, units = 'in', res = 900)
vioplot(Max_cov ~ Gene_name, data = df,
        names = rep("", nrow(df_summary)),
        xlab = "", ylab = "Max. Coverage", col = palette,colMed = "black", ylim = c(-5, 100)
        ,xlim = c(2, 52))

text(x = x_pos, y = par("usr")[3] - 1, labels = df_summary$Gene_name,
     srt = 30, adj = 1.1, xpd = TRUE, cex = 0.6)
abline(h = 70, lty = 2, col = "red")
dev.off()
###### Seq divegence ####
df <- read.table("./GC_content_seq_divergence/Identity_wrt_human.tsv", header = T, sep = '\t')
df <- na.omit(df)
library(ggplot2)
library(dplyr)
library(vioplot)
library(RColorBrewer)   
# Summarize the data
df_summary <- df %>%
  group_by(Gene_name) %>%
  summarise(
    median_val = median(Seq.Identity, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
df$Gene_name <- factor(df$Gene_name, levels = df_summary$Gene_name)
n_colors <- nlevels(df$Gene_name)  # Replace 'group_column' with your actual column representing groups
palette <- colorRampPalette(brewer.pal(9, "Set1"))(n_colors)

png("Seq_identity.png", width = 16, height = 9, units = 'in', res = 900)
vioplot(Seq.Identity ~ Gene_name, data = df,
        names = rep("", nrow(df_summary)),
        xlab = "", ylab = "Seq. Identity", col = palette,colMed = "black", ylim = c(-5, 100))

x_pos <- 1:nrow(df_summary)
axis(1, at = x_pos, labels = FALSE)
text(x = x_pos, y = par("usr")[3] - 2, labels = df_summary$Gene_name,
     srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
text(x = x_pos,
     y = 10,  
     labels = round(df_summary$median_val, 0),
     cex = 1, col = "red")+
  text(x = x_pos,
       y = 5, 
       labels = df_summary$n,
       cex = 0.7)
abline(h = 70, lty = 2, col = "red")
title(xlab = "Gene Names")
dev.off()

##### GC content #####
df <- read.table("./GC_content_seq_divergence/GC_content.tsv", header = T, sep = '\t')
df <- na.omit(df)
library(ggplot2)
library(dplyr)
library(vioplot)
library(RColorBrewer)   
# Summarize the data
df_summary <- df %>%
  group_by(Gene_name) %>%
  summarise(
    median_val = median(GC_content, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
df$Gene_name <- factor(df$Gene_name, levels = df_summary$Gene_name)
n_colors <- nlevels(df$Gene_name)  # Replace 'group_column' with your actual column representing groups
palette <- colorRampPalette(brewer.pal(9, "Set1"))(n_colors)

png("GC_content.png", width = 16, height = 9, units = 'in', res = 900)
vioplot(GC_content ~ Gene_name, data = df,
        names = rep("", nrow(df_summary)),
        xlab = "", ylab = "GC-content (%)", col = palette,colMed = "black", ylim = c(20, 80))

x_pos <- 1:nrow(df_summary)
axis(1, at = x_pos, labels = FALSE)
text(x = x_pos, y = par("usr")[3] - 2, labels = df_summary$Gene_name,
     srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
text(x = x_pos,
     y = 30,  
     labels = round(df_summary$median_val, 0),
     cex = 1, col = "red")+
  text(x = x_pos,
       y = 25, 
       labels = df_summary$n,
       cex = 0.7)
abline(h = 55, lty = 2, col = "red")
title(xlab = "Gene Names")
dev.off()

##### GC strech ####
df <- read.table("./GC_content_seq_divergence/GC_Stretch.tsv", header = T, sep = '\t')
df <- na.omit(df)
library(ggplot2)
library(dplyr)
library(vioplot)
library(RColorBrewer)   
# Summarize the data
df_summary <- df %>%
  group_by(Gene_name) %>%
  summarise(
    median_val = median(GC_Stretch, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
summary(df$GC_Stretch)
df$Gene_name <- factor(df$Gene_name, levels = df_summary$Gene_name)
n_colors <- nlevels(df$Gene_name)  # Replace 'group_column' with your actual column representing groups
palette <- colorRampPalette(brewer.pal(9, "Set1"))(n_colors)

png("GC_Stretch.png", width = 16, height = 9, units = 'in', res = 900)
vioplot(GC_Stretch ~ Gene_name, data = df,
        names = rep("", nrow(df_summary)),
        xlab = "", ylab = "GC-Stretch", col = palette,colMed = "black", ylim = c(2.6, 7))

x_pos <- 1:nrow(df_summary)
axis(1, at = x_pos, labels = FALSE)
text(x = x_pos, y = par("usr")[3] - 2, labels = df_summary$Gene_name,
     srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
text(x = x_pos,
     y = 2.8,  
     labels = round(df_summary$median_val, 0),
     cex = 1, col = "red")+
  text(x = x_pos,
       y = 2.7, 
       labels = df_summary$n,
       cex = 0.7)
abline(h = 3.936, lty = 2, col = "red")
title(xlab = "Gene Names")
dev.off()



############ ALL IN ONE ################
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)  # For combining ggplots

### Helper function to create violin plots ###
plot_violin <- function(df, value_col, ylab_text, hline = NULL, ylims = NULL, y_med = NULL, y_n = NULL) {
  df_summary <- df %>%
    group_by(Gene_name) %>%
    summarise(
      median_val = median(.data[[value_col]], na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  df$Gene_name <- factor(df$Gene_name, levels = df_summary$Gene_name)
  
  gg <- ggplot(df, aes(x = Gene_name, y = .data[[value_col]], fill = Gene_name)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_text(data = df_summary, aes(x = Gene_name, y = y_med, label = round(median_val, 0)), 
              inherit.aes = FALSE, color = "red", size = 3) +
    geom_text(data = df_summary, aes(x = Gene_name, y = y_n, label = n), 
              inherit.aes = FALSE, size = 2.5) +
    labs(x = "Gene Name", y = ylab_text) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          legend.position = "none") +
    scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(df$Gene_name)))) +
    coord_cartesian(ylim = ylims)
  
  if (!is.null(hline)) {
    gg <- gg + geom_hline(yintercept = hline, linetype = "dashed", color = "red")
  }
  return(gg)
}

### Plot 1: Max Coverage ###
df1 <- read.table("./Squamata_genone_blast/GENE_wise/Gene_species_Max_coverage.tsv", header = T, sep = '\t') %>% na.omit()
p1 <- plot_violin(df1, "Max_cov", "Max. Coverage", hline = 70, ylims = c(-5, 100), y_med = 90, y_n = 99)+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank())


### Plot 2: Seq Identity ###
df2 <- read.table("./GC_content_seq_divergence/Identity_wrt_human.tsv", header = T, sep = '\t') %>% na.omit()
p2 <- plot_violin(df2, "Seq.Identity", "Sequence Identity (%)", hline = 70, ylims = c(-5, 100), y_med = 15, y_n = 5)+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank())


### Plot 3: GC Content ###
df3 <- read.table("./GC_content_seq_divergence/GC_content.tsv", header = T, sep = '\t') %>% na.omit()
p3 <- plot_violin(df3, "GC_content", "GC-content (%)", hline = 55, ylims = c(20, 80), y_med = 30, y_n = 25)+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank())


### Plot 4: GC Stretch ###
df4 <- read.table("./GC_content_seq_divergence/GC_Stretch.tsv", header = T, sep = '\t') %>% na.omit()
p4 <- plot_violin(df4, "GC_Stretch", "GC-Stretch", hline = 3.936, ylims = c(2.6, 7), y_med = 3.1, y_n = 2.7)
### Combine all plots ###
final_plot <- (p1 / p2 / p3 / p4) + plot_layout(ncol = 1)

### Save to file ###
ggsave("All_violin_plots.png", final_plot, width = 16, height = 9, dpi = 900)

