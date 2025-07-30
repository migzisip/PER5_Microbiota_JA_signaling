######################hierarchical clustering####################

setwd('D:/Projects/PER5_RNA_seq/Root_DEGs/')

library('readr')
library('magrittr')
library('tibble')
library('gplots')
library('dendextend')
library('dynamicTreeCut')
library('ggplot2')
library('tidyr')
library('DESeq2')
library('dplyr')
library('RColorBrewer')
library('gridExtra')
library('cluster')
library('scales')
library('ComplexHeatmap')
library('magrittr')
library('tidyverse')
library('circlize')


deganno <- read_csv('root_All_DEGs.csv',
                    col_types = cols(Chromosome = col_character()))


deganno_logFC <- deganno %>%
  select(matches('ID|log2')) %>%
  column_to_rownames(var = 'ID')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#scale
scale_deganno_logFC <- deganno_logFC %>%
  t %>%
  scale %>%
  t
scale_deganno_logFC %<>% .[complete.cases(.), ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~cluster~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Cluster rows by Pearson correlation
hr <- scale_deganno_logFC %>%
  t %>%
  cor(method = 'pearson') %>%
  {1 - .} %>%
  as.dist %>%
  hclust(method = 'complete')

## Clusters columns by Spearman correlation
hc <- scale_deganno_logFC %>%
  cor(method = 'spearman') %>%
  {1 - .} %>%
  as.dist %>%
  hclust(method = 'complete')

cairo_pdf('heatmap.pdf')
heatmap.2(scale_deganno_logFC,
          Rowv = as.dendrogram(hr),
          Colv = as.dendrogram(hc),
          col = redgreen(100),
          scale = 'row',
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = 'Heatmap.2',
          trace = 'none')
dev.off()

hc %>%
  as.dendrogram(method = 'average') %>%
  plot(main = 'Sample Clustering',
       ylab = 'Height',
       mar = c(20, 2, 2, 2))


hr %>%
  as.dendrogram(method = 'average') %>%
  plot(leaflab = 'none',
       main = 'Gene Clustering',
       ylab = 'Height')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

##~~~~~~~~~~~~~~~~~~~~~~~~~K-means cluster~~~~~~~~~~~~~~~~~~~~~~
## choose groups
## 1. sum of squared error
wss <- (nrow(scale_deganno_logFC) - 1) * sum(apply(scale_deganno_logFC, 2, var))

for (i in 2:20) {
  wss[i] <- sum(kmeans(scale_deganno_logFC,
                       centers=i,
                       algorithm = 'MacQueen')$withinss)
}

ggplot(tibble(k = 1:20, wss = wss), aes(k, wss)) +
  geom_point(colour = '#D55E00', size = 3) +
  geom_line(linetype = 'dashed') +
  xlab('Number of clusters') +
  ylab('Sum of squared error')
ggsave('kmeans_sse.pdf')
ggsave('kmeans_sse.jpg')


## 2. Akaike information criterion
kmeansAIC = function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + 2*m*k)
}

aic <- numeric(20)
for (i in 1:20) {
  fit <- kmeans(x = scale_deganno_logFC, centers = i, algorithm = 'MacQueen')
  aic[i] <- kmeansAIC(fit)
}

ggplot(tibble(k = 1:20, aic = aic), aes(k, wss)) +
  geom_point(colour = '#009E73', size = 3) +
  geom_line(linetype = 'dashed') +
  xlab('Number of clusters') +
  ylab('Akaike information criterion')
ggsave('kmeans_AIC.pdf')
ggsave('kmeans_AIC.jpg')


set.seed(123)
clNum <- 8
kClust8 <- kmeans(scale_deganno_logFC, centers = clNum, algorithm= 'MacQueen', nstart = 1000, iter.max = 20)


cl <- kClust8$cluster
prefix <- paste0('kmeans', clNum)


clusterGene <- scale_deganno_logFC %>%
  as.data.frame %>%
  rownames_to_column(var = 'ID') %>%
  as_tibble %>%
  {
    cl <- as.data.frame(cl) %>%
      rownames_to_column(var = 'ID')
    inner_join(., cl)
  }


cl <- as.data.frame(cl) %>%
  rownames_to_column(var = 'ID')

kmeans8 <- inner_join(deganno, cl, by = "ID")


write.csv(kmeans8, file='Kmeans8.csv', row.names = FALSE)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~select DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fcsig <- kmeans8 %>%
  select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ -1,
                                 TRUE ~ 0)))
padjsig <- kmeans8 %>%
  select(ends_with('padj')) %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

heatsig <- (padjsig * fcsig) %>%
  as_tibble %>%
  abs %>%
  rowSums %>%
  {. >= 1} %>%
  which %>%
  dplyr::slice(kmeans8, .) #%>%
#inner_join(kmeansRes)

sigMat <- (padjsig * fcsig) %>%
  as_tibble %>%
  setNames(names(.) %>% substr(., start = 1, stop = nchar(.) - 5)) %>%
  mutate(ID = deganno$ID) %>%
  inner_join(heatsig %>% select(ID), .) %T>%
  {(sum(.$ID == heatsig$ID) == nrow(.)) %>% print} %>%
  transmute_at(.var = vars(contains('vs')),
               list(~ case_when(. == -1 ~ 'down',
                                . == 0 ~'no',
                                . == 1 ~ 'up'))) %>%
  as.matrix


write.csv(heatsig, file='kmeans8_sig.csv', row.names = FALSE)

# heatsig %>%
#   mutate_at(c('Gene', 'Description'), .funs = list(~if_else(is.na(.), '', .))) %>%
#   write_csv('kmeans12_sig.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


deganno_logFC_DEGs <- heatsig %>%
  select(matches('ID|log2')) %>%
  column_to_rownames(var = 'ID')

#scale
scale_deganno_logFC_DEGs <- deganno_logFC_DEGs %>%
  t %>%
  scale %>%
  t
scale_deganno_logFC_DEGs %<>% .[complete.cases(.), ]

scale_deganno_logFC_matrix <- as.data.frame(scale_deganno_logFC_DEGs) %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(heatsig %>% select(ID, Gene, Description, cl)) %>%
  select(ID, Gene, Description, everything())


colnames<- c("ID", "Gene", "Description", "per5_M_vs_Ws4_M",
             "Ws4_SC_vs_Ws4_M" ,   "per5_SC_vs_Ws4_SC",
              "per5_SC_vs_per5_M", "cl")

colnames(scale_deganno_logFC_matrix) <- colnames

write.csv(scale_deganno_logFC_matrix, file='kmeans8_sig_z_score.csv', row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#heat plot - with scaling



col_fun <- colorRamp2(seq(-3, 3, length.out = 100), colorRampPalette(c("blue", "white", "red"))(100))

# col_tissue_root <- HeatmapAnnotation(Tissue = c(rep(c('Root'), each = 1), rep(c('Shoot'), each = 4)),
#                                   col = list(Tissue = c('Root' = 'brown', 'Shoot' = 'green')),
#                                   gp = gpar(col = 'black'))


col_tissue_root <- HeatmapAnnotation(
  Tissue = rep('Root', 4),
  col = list("Tissue" = c('Root' = 'brown')),
  gp = gpar(col = 'black'),
  annotation_label = NULL, # hide the Tissue in heatmap
annotation_legend_param = list(
  Tissue = list(
    title_gp = gpar(fontsize = 16,  fontface = "bold"),  # Title font size
    title_position = "lefttop",
    #ncol = 1, 
    nrow = 1,
    labels_gp = gpar(fontsize = 14) # Labels font size
    )
  )
)


ht_z_score_root <- Heatmap(matrix = scale_deganno_logFC_matrix %>% select(contains('_')),
                                name = 'Z-score',
                                width = unit(10, "cm"),
                                #row_order = order(scaleC$cl) %>% rev,
                                row_split = scale_deganno_logFC_matrix$cl,
                                row_gap = unit(4, "mm"),
                                column_order = 1:4,
                                column_split = rep(c('Root'), c(4)),
                                show_column_names = T,
                                col = col_fun,
                                column_names_rot = -45,
                                top_annotation = c(col_tissue_root),
                                use_raster = FALSE,
                                heatmap_legend_param = list(
                                  at = -3:3,  # Set custom legend breaks
                                  labels = as.character(-3:3)
                                )) 
ht_z_score_root

png('kmeans8_heatmap_Z-score_root_PER5.jpg', width = 6, height = 10, units = 'in', res = 300)
draw(ht_z_score_root)
dev.off()

pdf('kmeans8_heatmap_Z-score_root_PER5.pdf', width = 6, height = 10)
draw(ht_z_score_root)
dev.off()


Heatmap(matrix = scale_deganno_logFC_matrix)

##Combine heatmap

ht_z_score_root_combinePlot <- Heatmap(matrix = scale_deganno_logFC_matrix %>% select(contains('_')),
                           name = 'Z-score',
                           width = unit(10, "cm"),
                           cluster_columns = F,
                           #row_order = order(scaleC$cl) %>% rev,
                           row_split = scale_deganno_logFC_matrix$cl,
                           row_gap = unit(4, "mm"),
                           #column_title_gp = gpar(fontsize = 18,  fontface = "bold"), #font size
                           column_title = NULL,  # Remove column titles
                           column_order = 1:4,
                           column_split = rep(c('Root'), c(4)),
                           column_gap = unit(2, 'mm'),
                           show_column_names = T,
                           col = col_fun,
                           column_names_rot = -45,
                           top_annotation = c(col_tissue_root),
                           row_title_gp = gpar(fontsize = 16, fontface = "bold"),
                           use_raster = FALSE,
                           heatmap_legend_param = list(
                             at = -3:3,  # Set custom legend breaks
                             direction = "horizontal",
                             title_position = "lefttop",
                             labels = as.character(-3:3),
                             title_gp = gpar(fontsize = 16,  fontface = "bold"), # Increase legend title font size
                             labels_gp = gpar(fontsize = 14) 
                           )) +
  Heatmap(sigMat,
          heatmap_width = unit(3, "cm"),
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 14),
          column_split = paste0('Root_', 1:4),
          column_names_rot = -45,
          column_title = NULL,
          border_gp = gpar(col = "black"),
          heatmap_legend_param = list(
            title = 'DEGs',
            title_gp = gpar(fontsize = 16,  fontface = "bold"), # Legend title font size
            labels_gp = gpar(fontsize = 14), # Legend label font size
            title_position = "lefttop",
            nrow = 1
          ),
          cluster_columns = FALSE,
          use_raster = FALSE)
  

ht_z_score_root_combinePlot


pdf('kmeans8_heatmap_Z-score_root_PER5_compbinePlot_v2.pdf', width = 8, height = 10)
draw(ht_z_score_root_combinePlot)
dev.off()


pdf('kmeans8_heatmap_Z-score_root_PER5_compbinePlot_v3.pdf', width = 8, height = 10)
draw(ht_z_score_root_combinePlot, 
     heatmap_legend_side = "bottom",  # Puts DEG and Z-score legend horizontally at the bottom
     annotation_legend_side = "bottom")  # Moves class annotation legend to the bottom
dev.off()





# Prepare DEG count data for barplot annotation
deg_counts <- data.frame(
  up = c(52, 252, 543, 771),    
  down = c(83, 297, 505, 651)
)

# Create barplot annotation for the second heatmap
deg_bar_annotation <- HeatmapAnnotation(
  DEG_Barplot = anno_barplot(
    deg_counts,
    bar_width = 0.8,
    gp = gpar(fill = c("red", "green")),
    axis_param = list(side = "right")
  ),
  #annotation_name_side = "top",
  height = unit(1, "cm") #3
)

# Main heatmap
ht_z_score_root_combinePlot_1 <- Heatmap(
  matrix = scale_deganno_logFC_matrix %>% select(contains('_')),
  name = 'Z-score',
  width = unit(10, "cm"),
  cluster_columns = F,
  row_split = scale_deganno_logFC_matrix$cl,
  row_gap = unit(4, "mm"),
  column_title = NULL,
  column_order = 1:4,
  column_split = rep(c('Root'), c(4)),
  column_gap = unit(2, 'mm'),
  show_column_names = T,
  col = col_fun,
  column_names_rot = 45,
  top_annotation = col_tissue_root,
  row_title_gp = gpar(fontsize = 16, fontface = "bold"),
  use_raster = FALSE,
  heatmap_legend_param = list(
    at = -3:3,
    direction = "horizontal",
    title_position = "lefttop",
    labels = as.character(-3:3),
    title_gp = gpar(fontsize = 16,  fontface = "bold"),
    labels_gp = gpar(fontsize = 14)
  )
) +
  # Second heatmap with DEG barplot annotation
  Heatmap(
    sigMat,
    heatmap_width = unit(3, "cm"),
    col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
    column_names_gp = gpar(fontsize = 14),
    column_split = paste0('Root_', 1:4),
    column_names_rot = -45,
    column_title = NULL,
    border_gp = gpar(col = "black"),
    
    # Add barplot on top!
    top_annotation = deg_bar_annotation,
    
    heatmap_legend_param = list(
      title = 'DEGs',
      title_gp = gpar(fontsize = 16,  fontface = "bold"),
      labels_gp = gpar(fontsize = 14),
      title_position = "lefttop",
      nrow = 1
    ),
    cluster_columns = FALSE,
    use_raster = FALSE
  )

ht_z_score_root_combinePlot_1

pdf('kmeans8_heatmap_Z-score_root_PER5_compbinePlot_v4.pdf', width = 8, height = 10)
draw(ht_z_score_root_combinePlot_1, 
     heatmap_legend_side = "bottom",  # Puts DEG and Z-score legend horizontally at the bottom
     annotation_legend_side = "bottom")  # Moves class annotation legend to the bottom
dev.off()


