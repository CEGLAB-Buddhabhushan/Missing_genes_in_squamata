library(ggpubr)
library(dplyr)
library(ggplot2)
getwd()
setwd("/home/sonal/Documents/New_gc/output")
data_gc=read.table("GC_content-GC_Stretch.tsv", header = T, sep = '\t')
View(data_gc)
data_gc$Gene_name <- toupper(data_gc$Gene_name)
Group = unique(data_gc$Group[data_gc$Gene_name=="STAP1"])
Group


colors = c("#ff7f00",  "#e31a1c","#33a02c", "#0000FF", "#FFD700",
           "#F0E68C", "#800080", "brown", "#fb9a99", "#fdbf6f")
shapes <- c(19, 17, 18, 15, 0, 1, 2, 21, 4, 5, 6, 7, 8, 9, 10, 11)
color_mapping <- setNames(colors, Group)
shape_mapping <- setNames(shapes, Group)
fill_mapping <- setNames(colors, Group)

#Scatter Plot
TNIP2_scatter_plot <- ggplot() +
  geom_point(data = data_gc[data_gc$Gene_name=="TNIP2", ], aes(x = GC_Stretch, y = GC_content), alpha = 0.02, color = "gray") +
  geom_point(data = data_gc[data_gc$Gene_name=="TNIP2", ], aes(x = GC_Stretch, y = GC_content, color = Group, shape = Group), size = 3) + 
  scale_color_manual(values = color_mapping, name = "Group") +
  scale_shape_manual(values = shape_mapping, name = "Group") +
  labs(title = "GC-content vs Avg. G/C-stretch length",
       x = "Average length of G/C-stretches",
       y = "GC content (in %)") +
  theme_classic() +
  #scale_y_continuous(limits = c(30, 80))+
  #scale_x_continuous(limits = c(2.7, 7)) +
  
  theme(legend.position = "bottom", legend.box = "horizontal",
        axis.text = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, title = NULL),
         shape = guide_legend(nrow = 2, byrow = TRUE, title = NULL))

##
LAPTM5_scatter_plot <- ggplot() +
  geom_point(data = data_gc[data_gc$Gene_name=="LAPTM5", ], aes(x = GC_Stretch, y = GC_content), alpha = 0.02, color = "gray") +
  geom_point(data = data_gc[data_gc$Gene_name=="LAPTM5", ], aes(x = GC_Stretch, y = GC_content, color = Group, shape = Group), size = 3) + 
  scale_color_manual(values = color_mapping, name = "Group") +
  scale_shape_manual(values = shape_mapping, name = "Group") +
  labs(title = "GC-content vs Avg. G/C-stretch length",
       x = "Average length of G/C-stretches",
       y = "GC content (in %)") +
  theme_classic() +
  scale_y_continuous(limits = c(30, 80))+
  scale_x_continuous(limits = c(2.7, 7)) +
  
  theme(legend.position = "bottom", legend.box = "horizontal",
        axis.text = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, title = NULL),
         shape = guide_legend(nrow = 2, byrow = TRUE, title = NULL))

##
STAP1_scatter_plot <- ggplot() +
  geom_point(data = data_gc[data_gc$Gene_name=="STAP1", ], aes(x = GC_Stretch, y = GC_content), alpha = 0.02, color = "gray") +
  geom_point(data = data_gc[data_gc$Gene_name=="STAP1", ], aes(x = GC_Stretch, y = GC_content, color = Group, shape = Group), size = 3) + 
  scale_color_manual(values = color_mapping, name = "Group") +
  scale_shape_manual(values = shape_mapping, name = "Group") +
  labs(title = "GC-content vs Avg. G/C-stretch length",
       x = "Average length of G/C-stretches",
       y = "GC content (in %)") +
  theme_classic() +
  scale_y_continuous(limits = c(30, 80))+
  scale_x_continuous(limits = c(2.7, 7)) +
  
  theme(legend.position = "bottom", legend.box = "horizontal",
        axis.text = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, title = NULL),
         shape = guide_legend(nrow = 2, byrow = TRUE, title = NULL))
TNIP2_scatter_plot
STAP1_scatter_plot
LAPTM5_scatter_plot

#### Boxplots : 
TNIP2_boxplot<- ggplot(data_gc[data_gc$Gene_name=="TNIP2", ], aes(x = Group, y = GC_content)) +
  geom_boxplot() +
  stat_summary(fun = "median", geom = "text", aes(label = paste(round(after_stat(y), 2))), vjust = -0.5, position = position_dodge(width = 0.75), size = 4)+
  labs(title = "GC-Content by Group", x = "Groups", y = "GC %") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red", lwd = 1)
##STAP1
STAP1_boxplot<- ggplot(data_gc[data_gc$Gene_name=="STAP1", ], aes(x = Group, y = GC_content)) +
  geom_boxplot() +
  stat_summary(fun = "median", geom = "text", aes(label = paste(round(after_stat(y), 2))), vjust = -0.5, position = position_dodge(width = 0.75), size = 4)+
  labs(title = "GC-Content by Group", x = "Groups", y = "GC %") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red", lwd = 1)

LAPTM5_boxplot<- ggplot(data_gc[data_gc$Gene_name=="LAPTM5", ], aes(x = Group, y = GC_content)) +
  geom_boxplot() +
  stat_summary(fun = "median", geom = "text", aes(label = paste(round(after_stat(y), 2))), vjust = -0.5, position = position_dodge(width = 0.75), size = 4)+
  labs(title = "GC-Content by Group", x = "Groups", y = "GC %") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red", lwd = 1)

TNIP2_boxplot
STAP1_boxplot
LAPTM5_boxplot

TNIP2=ggarrange(TNIP2_scatter_plot,TNIP2_boxplot, ncol = 1, nrow = 2, labels = c("A", "C"), common.legend = TRUE, legend = "bottom", align = "hv")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE, title = NULL))
STAP1=ggarrange(STAP1_scatter_plot,STAP1_boxplot, ncol = 1, nrow = 2, labels = c("A", "C"), common.legend = TRUE, legend = "bottom", align = "hv")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE, title = NULL))
LAPTM5=ggarrange(LAPTM5_scatter_plot,LAPTM5_boxplot, ncol = 1, nrow = 2, labels = c("A", "C"), common.legend = TRUE, legend = "bottom", align = "hv")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE, title = NULL))

TNIP2
STAP1
LAPTM5
