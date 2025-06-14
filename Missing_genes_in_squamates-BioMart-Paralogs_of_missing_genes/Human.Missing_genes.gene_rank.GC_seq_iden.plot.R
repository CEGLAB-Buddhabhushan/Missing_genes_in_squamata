getwd()
############# Gene rank plot #################
library(ggplot2)
library(dplyr)
library(tidyverse)
library(forcats)
library(ggpubr)
# Read the main file with missing genes and paralog info
missing_paralog <- read.delim("Human.Missing_genes-Paralog_info.mart_export.txt", sep = "\t", stringsAsFactors = FALSE)
head(missing_paralog)
# Read the file with paralog GC content info
paralog_gc <- read.delim("Human.paralog_GC_content.mart_export.txt", sep = "\t", stringsAsFactors = FALSE)
head(paralog_gc)
# Rename columns in GC content file for clarity and to match merge key
colnames(paralog_gc)[colnames(paralog_gc) == "Gene.stable.ID"] <- "Human.paralogue.gene.stable.ID"
colnames(paralog_gc)[colnames(paralog_gc) == "Gene...GC.content"] <- "Paralogue.GC.content"

# Merge the two data frames based on Human paralogue gene stable ID
merged_data <- merge(missing_paralog, paralog_gc[, c("Human.paralogue.gene.stable.ID", "Paralogue.GC.content")],
                     by.x = "Human.paralogue.gene.stable.ID", by.y = "Human.paralogue.gene.stable.ID",
                     all.x = TRUE)
head(merged_data)
# Write the merged data to a new file
write.table(merged_data, "Human.Missing_genes-Paralog_info_with_GC.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



data <- read.table("./Human.Missing_genes-Paralog_info_with_GC.tsv", header = T, sep = "\t")
data<-na.omit(data)
df<- data[,c(1,2,4,6,15,19)]
head(df)
df_unique <- df[!duplicated(df), ] 
head(df_unique)
##############
gene_gc_ranked <- df_unique %>%
  select(Gene.stable.ID, gene_gc = Gene...GC.content) %>%
  distinct() %>%
  arrange(gene_gc) %>%
  mutate(Gene_stable_ID_rank = row_number())
length(unique(gene_gc_ranked$Gene_stable_ID_rank))
head(gene_gc_ranked)
# Step 2: Rank unique paralogue GC content
paralog_gc_ranked <- df_unique %>%
  select(Human.paralogue.gene.stable.ID, gene_gc = Paralogue.GC.content) %>%
  distinct() %>%
  arrange(gene_gc) %>%
  mutate(Human.paralogue_gene_stable_ID_rank = row_number())
length(unique(paralog_gc_ranked$Human.paralogue_gene_stable_ID_rank))
head(paralog_gc_ranked)



df_plot <- df_unique %>%
  left_join(gene_gc_ranked, by = "Gene.stable.ID") %>%
  left_join(paralog_gc_ranked, by = "Human.paralogue.gene.stable.ID")

head(df_plot)
View(df_plot)
##############
df_plot <- df_plot %>%
  mutate(
    identity_bin_raw = cut(
      Paralogue..id..target.Human.gene.identical.to.query.gene,
      breaks = c(0, 25, 50, 75, 100),
      include.lowest = TRUE,
      right = FALSE  # [0–25), [25–50), etc.
    )
  )

# Step 2: Create readable labels
df_plot <- df_plot %>%
  mutate(
    identity_bin_labeled = recode_factor(identity_bin_raw,
                                         "[0,25)" = "0–25%",
                                         "[25,50)" = "25–50%",
                                         "[50,75)" = "50–75%",
                                         "[75,100]" = "75–100%"
    )
  )

# Step 3: Reorder levels (for legend ordering)
df_plot <- df_plot %>%
  mutate(
    identity_bin = fct_relevel(identity_bin_labeled, "0–25%", "25–50%", "50–75%", "75–100%")
  )


##### for all####
all_label_df <- df_plot %>%
  group_by(Gene.name) %>%
  slice_max(order_by = Paralogue..id..target.Human.gene.identical.to.query.gene, n = 1) %>%
  ungroup()
library(ggplot2)
library(ggpubr)

All_with_label<-ggplot(df_plot, aes(x = Gene_stable_ID_rank, y = Human.paralogue_gene_stable_ID_rank)) +
  geom_point(aes(
    color = Gene.name,
    size = identity_bin
  ), alpha = 0.5) +
  geom_text(
    data = all_label_df,
    aes(label = Gene.name),
    size = 2.5,hjust=1.5,
    check_overlap = TRUE  # avoids overlapping labels
  ) +
  geom_hline(yintercept = 165, linetype = "dashed", color = "blue") +
  annotate("text", x = 5, y = 169 ,
           label = "GC-content = 55.30 %", color = "blue", hjust = 1, size = 4) +
  geom_vline(xintercept = 36, linetype = "dashed", color = "red") +
  annotate("text", x = 36.7, y = 20,
           label = "GC-content = 55.56 %", color = "red", angle = 90, vjust = -0.5, size = 4) +
  labs(
    x = "Gene stable ID (ranked by GC content)",
    y = "Human paralogue transcript ID (ranked by GC content)",
    size = "% Sequence identity (Human gene vs Paralog)"
  ) +
  theme_classic2() +
  theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 5, byrow = TRUE))


png("Human.Missing_genes.gene_rank.GC_seq_iden.plot.png", units="in", width=17, height=9.5, res=1500)
All_with_label
dev.off()



write.table(df_plot, file = "human.df_plot.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
