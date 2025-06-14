library(ggpubr)
library(dplyr)
library(ggplot2)
getwd()
setwd("/home/sonal/Documents/New_gc/output")
data_gc=read.table("GC_content-GC_Stretch.tsv", header = T, sep = '\t')
View(data_gc)
data_gc$Gene_name <- toupper(data_gc$Gene_name)
desired_groups <- c("Aves", "Crocodylia", "Mammalia", "Testudines", "Amphibia")
filtered_data <- data_gc[data_gc$Group %in% desired_groups, ]
dev.new()

library(vioplot)
library(dplyr)

png("blue_gc_content_plot.png", width = 6.5, height = 2.4, units = "in", res = 250)

genes_uni <- unique(filtered_data$Gene_name)
groups_uni <- unique(filtered_data$Group)
groups_uni
group_colors <- c("Amphibia" = "#F8766D", "Aves" = "#B79F00", 
                  "Crocodylia" = "#00BA38", "Mammalia" = "#00BFC4", 
                  "Testudines" = "#F564E3")

ncols <- 4
nrows <- 1
layout_matrix <- matrix(1:(nrows*ncols), nrow = nrows, ncol = ncols, byrow = TRUE)
layout(rbind(layout_matrix, rep(nrows*ncols+1, ncols)), 
       heights = c(rep(4, nrows), 1))

par(mar = c(2, 4, 2, 1))

for(gene in genes_uni) {
  gene_data <- filtered_data %>% 
    filter(Gene_name == gene) %>%
    split(.$Group)
  plot_data <- lapply(groups_uni, function(g) {
    gene_data[[g]]$GC_content[!is.na(gene_data[[g]]$GC_content)]
  })
  names(plot_data) <- groups_uni
  
  
  plot_data <- plot_data[sapply(plot_data, length) > 0]
  
  
  y_lim <- range(filtered_data$GC_content, na.rm = TRUE) * c(0.95, 1.05)
  
  
  vioplot::vioplot(plot_data, col = group_colors[names(plot_data)],
                   ylim = y_lim, main = gene,
                   border = "black", rectCol = "brown",
                   drawRect = TRUE, xaxt="n")
  mtext("GC content (%)", side = 2, line = 2.5, cex= 0.5)
  
  abline(h = 50, lty = 2, col = "red")
  sapply(seq_along(plot_data), function(i) {
    points(i, mean(plot_data[[i]]), pch = 16, col = "red")
    text(i, mean(plot_data[[i]]), 
         labels = round(mean(plot_data[[i]]), 2), 
         pos = 3, cex = 0.6, col="blue"
    )
  })
}


par(mar = c(0,0,0,0))
plot.new()
legend("center", legend = names(group_colors), fill = group_colors,
       ncol = 5, bty = "n", title = "Taxonomic Groups")

dev.off()

dev.new()
