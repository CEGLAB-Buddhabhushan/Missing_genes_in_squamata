setwd("../GENE_wise/")
getwd()
df <- read.table("Gene_species_Max_coverage.tsv", header = T, sep = '\t')
df <- na.omit(df)
library(ggplot2)
library(dplyr)
library(vioplot)
# Summarize the data
df_summary <- df %>%
  group_by(Gene_name) %>%
  summarise(
    median_val = median(Max_cov, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

df$Gene_name <- factor(df$Gene_name, levels = df_summary$Gene_name)

png("Max_coverage.png", width = 16, height = 9, units = 'in', res = 900)
vioplot(Max_cov ~ Gene_name, data = df,
        names = rep("", nrow(df_summary)),
        xlab = "", ylab = "Max. Coverage", col = "lightblue",colMed = "black", ylim = c(-5, 100))

x_pos <- 1:nrow(df_summary)
axis(1, at = x_pos, labels = FALSE)
text(x = x_pos, y = par("usr")[3] - 2, labels = df_summary$Gene_name,
     srt = 45, adj = 1, xpd = TRUE, cex = 0.8)
text(x = x_pos,
     y = df_summary$median_val,  
     labels = df_summary$median_val,
     cex = 1, col = "red", adj = c(-0.5,0.5))+
text(x = x_pos,
     y = 95, 
     labels = df_summary$n,
     cex = 0.5)
abline(h = 70, lty = 2, col = "red")
title(xlab = "Gene Names")
dev.off()

##### GC content ####
setwd("../../GC_content_seq_divergence/")
df <- read.table("Identity_wrt_human.tsv", header = T, sep = '\t')
df <- na.omit(df)
library(ggplot2)
library(dplyr)
library(vioplot)
head(df)
# Summarize the data
df_summary <- df %>%
  group_by(Gene_name) %>%
  summarise(
    median_val = median(Seq.Identity, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

df$Gene_name <- factor(df$Gene_name, levels = df_summary$Gene_name)

png("Seq.Identity.png", width = 16, height = 9, units = 'in', res = 900)
vioplot(Seq.Identity ~ Gene_name, data = df,
        names = rep("", nrow(df_summary)),
        xlab = "", ylab = "Seq. Identity", col = "lightblue",colMed = "black", ylim = c(-5, 100))

x_pos <- 1:nrow(df_summary)
axis(1, at = x_pos, labels = FALSE)
text(x = x_pos, y = par("usr")[3] - 2, labels = df_summary$Gene_name,
     srt = 45, adj = 1, xpd = TRUE, cex = 0.8)
text(x = x_pos,
     y = df_summary$median_val,  
     labels = df_summary$median_val,
     cex = 1, col = "red", adj = c(0,0.5))+
  text(x = x_pos,
       y = 95, 
       labels = df_summary$n,
       cex = 0.5)
abline(h = 70, lty = 2, col = "red")
title(xlab = "Gene Names")
dev.off()
