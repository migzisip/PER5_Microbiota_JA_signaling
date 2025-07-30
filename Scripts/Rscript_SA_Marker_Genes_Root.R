
library('RColorBrewer')
library('ComplexHeatmap')
library('scales')
library('circlize')
library('tidyr')
library('tidyverse')
library('magrittr')


setwd('D:/Projects/PER5_RNA_seq/Root_DEGs/Marker_Genes')

marker_genes_df <- read.csv("Marker_genes_Roots_df.csv")

deganno_logFC <- marker_genes_df %>%
  select(matches('Gene|log2')) %>%
  column_to_rownames(var = 'Gene')

scale_deganno_logFC <- deganno_logFC %>%
  t %>%
  scale %>%
  t
scale_deganno_logFC %<>% .[complete.cases(.), ]



arrange_df <- scale_deganno_logFC[, c(1, 3, 2, 4)]


p2<-Heatmap(matrix = arrange_df,
        name = 'Z-score',
        cluster_columns = FALSE,
        column_order = 1:4,
        column_title = NULL,
        show_column_names = TRUE,
        column_split = factor(colnames(arrange_df), levels = colnames(arrange_df)),
        row_split = marker_genes_df$Marker)


p2

pdf('marker_genes_Root.pdf', width = 8, height = 10)
draw(p2)
#heatmap_legend_side = "bottom",  # Puts DEG and Z-score legend horizontally at the bottom
#annotation_legend_side = "bottom")  # Moves class annotation legend to the bottom
dev.off()






