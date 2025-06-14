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

png(filename="MAIN_FIGURE_COV_VIOLINEap.png", width = 16, height = 3,units = "in",res = 900)
vioplot(Max_cov ~ Gene_name, data = df,
        names = rep("", nrow(df_summary)),
        xlab = "", ylab = "Max. Coverage", col = palette,colMed = "black", ylim = c(0, 100)
        ,xlim = c(2.15, 51.85))
abline(h = 70, lty = 2, col = "red")
dev.off()
