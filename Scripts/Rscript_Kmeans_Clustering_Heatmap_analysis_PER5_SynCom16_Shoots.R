######################hierarchical clustering####################

setwd('D:/Projects/PER5_RNA_seq/Shoot_DEGs/')

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


deganno <- read_csv('shoot_All_DEGs.csv',
                    col_types = cols(Chromosome = col_character()))


deganno_logFC <- deganno %>%
  select(matches('ID|log2')) %>%
  column_to_rownames(var = 'ID')

# col_names <- c("SynCom_Live", "Actinomycetales_Live", "Flavobacteriales_Live", 
#                "Bacillales_Live", "Rhizobiales_Live", "Caulobacterales_Live",
#                "Sphingomonadales_Live", "Burkholderiales_Live",
#                "Xanthomonadales_Live", "Pseudomonadales_Live",
#                "SynCom_HK", "Actinomycetales_HK", "Flavobacteriales_HK",
#                "Bacillales_HK", "Rhizobiales_HK", "Caulobacterales_HK",
#                "Sphingomonadales_HK", "Burkholderiales_HK",
#                "Xanthomonadales_HK", "Pseudomonadales_HK")

#colnames(deganno_logFC) <- col_names

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

# col_tissue_shoot <- HeatmapAnnotation(Tissue = c(rep(c('Root'), each = 1), rep(c('Shoot'), each = 4)),
#                                      col = list(Tissue = c('Root' = 'brown', 'Shoot' = 'green')),
#                                      gp = gpar(col = 'black'))


col_tissue_shoot <- HeatmapAnnotation(
  Tissue = rep('Shoot', 4),
  col = list("Tissue" = c('Shoot' = 'green')),
  #gp = gpar(col = 'black'),
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






ht_z_score_shoot <- Heatmap(matrix = scale_deganno_logFC_matrix %>% select(contains('_')),
                           name = 'Z-score',
                           width = unit(10, "cm"),
                           #row_order = order(scaleC$cl) %>% rev,
                           row_split = scale_deganno_logFC_matrix$cl,
                           row_gap = unit(4, "mm"),
                           column_order = 1:4,
                           column_split = rep(c('Shoot'), c(4)),
                           show_column_names = T,
                           col = col_fun,
                           column_names_rot = -45,
                           top_annotation = c(col_tissue_shoot),
                           use_raster = FALSE,
                           heatmap_legend_param = list(
                             at = -3:3,  # Set custom legend breaks
                             labels = as.character(-3:3)
                           )) 
ht_z_score_shoot

png('kmeans8_heatmap_Z-score_shoot_PER5.jpg', width = 8, height = 10, units = 'in', res = 1200)
draw(ht_z_score_shoot)
dev.off()


pdf('kmeans8_heatmap_Z-score_shoot_PER5.pdf', width = 8, height = 10)
draw(ht_z_score_shoot)
dev.off()




##Combine heatmap

#ht_z_score_root <- 

ht_z_score_shoot_combinePlot <- Heatmap(matrix = scale_deganno_logFC_matrix %>% select(contains('_')),
                                       name = 'Z-score',
                                       width = unit(10, "cm"),
                                       cluster_columns = F,
                                       #row_order = order(scaleC$cl) %>% rev,
                                       row_split = scale_deganno_logFC_matrix$cl,
                                       row_gap = unit(4, "mm"),
                                       #column_title_gp = gpar(fontsize = 18,  fontface = "bold"), #font size
                                       column_title = NULL,  # Remove column titles
                                       column_order = 1:4,
                                       column_split = rep(c('Shoot'), c(4)),
                                       column_gap = unit(2, 'mm'),
                                       show_column_names = T,
                                       col = col_fun,
                                       column_names_rot = -45,
                                       top_annotation = c(col_tissue_shoot),
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


ht_z_score_shoot_combinePlot



png('kmeans8_heatmap_Z-score_root_PER5_compbinePlot.jpg', width = 6, height = 10, units = 'in', res = 1200)
draw(ht_z_score_root_combinePlot)
dev.off()


pdf('kmeans8_heatmap_Z-score_root_PER5_compbinePlot.pdf', width = 12, height = 13)
draw(ht_z_score_root_combinePlot)
dev.off()

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
  up = c(9, 348, 315, 1090),    
  down = c(5, 82, 89, 399)
)

# Create barplot annotation for the second heatmap
deg_bar_annotation <- HeatmapAnnotation(
  DEG_Barplot = anno_barplot(
    deg_counts,
    bar_width = 0.8,
    gp = gpar(fill = c("red", "blue")),
    axis_param = list(side = "right")
  ),
  #annotation_name_side = "top",
  height = unit(1, "cm") #3
)

# Main heatmap
ht_z_score_shoot_combinePlot_1 <- Heatmap(
  matrix = scale_deganno_logFC_matrix %>% select(contains('_')),
  name = 'Z-score',
  width = unit(10, "cm"),
  cluster_columns = F,
  row_split = scale_deganno_logFC_matrix$cl,
  row_gap = unit(4, "mm"),
  column_title = NULL,
  column_order = 1:4,
  column_split = rep(c('Shoot'), c(4)),
  column_gap = unit(2, 'mm'),
  show_column_names = T,
  col = col_fun,
  column_names_rot = 45,
  top_annotation = col_tissue_shoot,
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
    column_split = paste0('Shoot_', 1:4),
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

ht_z_score_shoot_combinePlot_1

pdf('kmeans8_heatmap_Z-score_shoot_PER5_compbinePlot.pdf', width = 8, height = 10)
draw(ht_z_score_shoot_combinePlot_1, 
     heatmap_legend_side = "bottom",  # Puts DEG and Z-score legend horizontally at the bottom
     annotation_legend_side = "bottom")  # Moves class annotation legend to the bottom
dev.off()

















_________




col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun_1 <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

col_order <-  HeatmapAnnotation(Order = c(rep(c("SynCom", "Actinomycetales", "Flavobacteriales", 
                                                "Bacillales", "Rhizobiales", "Caulobacterales", "Sphingomonadales", 
                                                "Burkholderiales", "Xanthomonadales", "Pseudomonadales"), each = 1),
                                          rep(c("SynCom", "Actinomycetales", "Flavobacteriales", 
                                                "Bacillales", "Rhizobiales", "Caulobacterales", "Sphingomonadales", 
                                                "Burkholderiales", "Xanthomonadales", "Pseudomonadales"), each = 1)),
                                col = list(Order = c("SynCom" = "#17becf",
                                                     "Actinomycetales" = "#1f77b4",
                                                     "Flavobacteriales" = "#9467bd",
                                                     "Bacillales" = "#ff7f0e",
                                                     "Rhizobiales" = "#7f7f7f",
                                                     "Caulobacterales" = "#d62728",
                                                     "Sphingomonadales" = "#bcbd22",
                                                     "Burkholderiales" = "#2ca02c",
                                                     "Xanthomonadales" = "#ff9896",
                                                     "Pseudomonadales" = "#e377c2" )),
                                gp = gpar(col = 'black'))

col_class <-  HeatmapAnnotation(class = c(rep(c("SynCom", "Actinobacteria", "Flavobacteriia", 
                                                "Bacilli", "Alphaproteobacteria", "Alphaproteobacteria", "Alphaproteobacteria", 
                                                "Betaproteobacteria", "Gammaproteobacteria", "Gammaproteobacteria"), each = 1),
                                          rep(c("SynCom", "Actinobacteria", "Flavobacteriia", 
                                                "Bacilli", "Alphaproteobacteria", "Alphaproteobacteria", "Alphaproteobacteria", 
                                                "Betaproteobacteria", "Gammaproteobacteria", "Gammaproteobacteria"), each = 1)),
                                col = list(class = c("SynCom" = "#999999",
                                                     "Actinobacteria" = "#FF9999",
                                                     "Flavobacteriia" = "#99FFFF",
                                                     "Bacilli" = "#FFCC66",
                                                     "Alphaproteobacteria" = "#66FF00",
                                                     "Alphaproteobacteria" = "#66FF00",
                                                     "Alphaproteobacteria" = "#66FF00",
                                                     "Betaproteobacteria" = "#00C00C",
                                                     "Gammaproteobacteria" = "#006006",
                                                     "Gammaproteobacteria" = "#006006" )),
                                gp = gpar(col = 'black'))

#Top annotation class and order
ht_scale <- Heatmap(matrix = scale_deganno_logFC_matrix %>% select(contains('_')),
                    name = 'z-score',
                    cluster_columns = T,
                    show_column_dend = F,
                    cluster_column_slices = F,
                    row_gap = unit(4, "mm"),
                    row_split = scale_deganno_logFC_matrix$cl,
                    show_column_names = F,
                    column_split = factor(rep(c("Root - Live bacteria", "Root - Heat killed bacteria"), c(10, 10)),
                                          levels = c("Root - Live bacteria", "Root - Heat killed bacteria")),
                    column_order = 1 : 20,
                    column_gap = unit(2, 'mm'),
                    col = col_fun,
                    top_annotation = c(col_class), 
                    #column_title_gp = gpar(fontsize = 12),
                    use_raster = FALSE,
                    heatmap_legend_param = list(
                      at = -2:2,  # Set custom legend breaks
                      labels = as.character(-2:2)
                    ))

ht_scale

png('kmeans12_heatmap_sig_DEG_Scaled_class_clusterbyCols_scale3.jpg', width = 12, height = 13, units = 'in', res = 1200)
draw(ht_scale)
dev.off()






#group based on phylogeny
indices <- c(1, 2, 4, 5, 3, 6, 7, 8, 9, 10, 11, 12, 14, 15, 13, 16, 17, 18, 19, 20, 21, 22)
arranged_df <- scale_deganno_logFC_matrix[indices]

col_class_1 <-  HeatmapAnnotation(class = c(rep(c("SynCom", "Flavobacteriia", "Bacilli", "Actinobacteria",  
                                                  "Alphaproteobacteria", "Alphaproteobacteria", "Alphaproteobacteria", 
                                                  "Betaproteobacteria", "Gammaproteobacteria", "Gammaproteobacteria"), each = 1),
                                            rep(c("SynCom", "Flavobacteriia", "Bacilli", "Actinobacteria",  
                                                  "Alphaproteobacteria", "Alphaproteobacteria", "Alphaproteobacteria", 
                                                  "Betaproteobacteria", "Gammaproteobacteria", "Gammaproteobacteria"), each = 1)),
                                  col = list(class = c("SynCom" = "#999999",
                                                       "Flavobacteriia" = "#99FFFF",
                                                       "Bacilli" = "#FFCC66",
                                                       "Actinobacteria" = "#FF9999",
                                                       "Alphaproteobacteria" = "#66FF00",
                                                       "Alphaproteobacteria" = "#66FF00",
                                                       "Alphaproteobacteria" = "#66FF00",
                                                       "Betaproteobacteria" = "#00C00C",
                                                       "Gammaproteobacteria" = "#006006",
                                                       "Gammaproteobacteria" = "#006006" )),
                                  gp = gpar(col = 'black'))

ht_scale_phylo <- Heatmap(matrix = arranged_df %>% select(contains('_')),
                          name = 'z-score',
                          cluster_columns = F,
                          #show_column_dend = F,
                          cluster_column_slices = F,
                          row_gap = unit(4, "mm"),
                          row_split = arranged_df$cl,
                          show_column_names = F,
                          column_split = factor(rep(c("Root - Live bacteria", "Root - Heat killed bacteria"), c(10, 10)),
                                                levels = c("Root - Live bacteria", "Root - Heat killed bacteria")),
                          column_order = 1 : 20,
                          column_gap = unit(2, 'mm'),
                          col = col_fun,
                          top_annotation = c(col_class_1), 
                          #column_title_gp = gpar(fontsize = 12),
                          use_raster = FALSE,
                          heatmap_legend_param = list(
                            at = -2:2,  # Set custom legend breaks
                            labels = as.character(-2:2)
                          ))

ht_scale_phylo

png('kmeans12_heatmap_sig_DEG_Scaled_class_PhyloBased_scale2.jpg', width = 12, height = 13, units = 'in', res = 1200)
draw(ht_scale_phylo)
dev.off()





#heatplot - no scaling
ht_NO_scale<- Heatmap(matrix = heatsig %>% select(contains('_log2')),
                      name = 'log2 FC',
                      cluster_columns = F,
                      #show_column_dend = F,
                      row_gap = unit(4, "mm"),
                      row_split = heatsig$cl,
                      show_column_names = F,
                      #column_split = rep(c("Root - Live bacteria", "Root - Heat killed bacteria"), c(10, 10)),
                      column_order = 1 : 20,
                      col = col_fun,
                      top_annotation = c(col_order), 
                      #column_title_gp = gpar(fontsize = 12),
                      use_raster = FALSE)

ht_NO_scale

png('kmeans12_heatmap_sig_DEG_NoScaled.jpg', width = 12, height = 13, units = 'in', res = 1200)
draw(ht_NO_scale)
dev.off()


#heatplot - no scaling
col_order <-  HeatmapAnnotation(Order = c(rep(c("SynCom", "Actinomycetales", "Flavobacteriales", 
                                                "Bacillales", "Rhizobiales", "Caulobacterales", "Sphingomonadales", 
                                                "Burkholderiales", "Xanthomonadales", "Pseudomonadales"), each = 1),
                                          rep(c("SynCom", "Actinomycetales", "Flavobacteriales", 
                                                "Bacillales", "Rhizobiales", "Caulobacterales", "Sphingomonadales", 
                                                "Burkholderiales", "Xanthomonadales", "Pseudomonadales"), each = 1)),
                                col = list(Order = c("SynCom" = "#17becf",
                                                     "Actinomycetales" = "#1f77b4",
                                                     "Flavobacteriales" = "#9467bd",
                                                     "Bacillales" = "#ff7f0e",
                                                     "Rhizobiales" = "#7f7f7f",
                                                     "Caulobacterales" = "#d62728",
                                                     "Sphingomonadales" = "#bcbd22",
                                                     "Burkholderiales" = "#2ca02c",
                                                     "Xanthomonadales" = "#ff9896",
                                                     "Pseudomonadales" = "#e377c2" )),
                                gp = gpar(col = 'black'))

ht_NO_scale<- Heatmap(matrix = heatsig %>% select(contains('_log2')),
                      name = 'log2 FC',
                      cluster_columns = F,
                      #show_column_dend = F,
                      row_gap = unit(4, "mm"),
                      row_split = heatsig$cl,
                      show_column_names = F,
                      #column_split = rep(c("Shoot - Live bacteria", "Shoot - Heat killed bacteria"), c(10, 10)),
                      column_order = 1 : 20,
                      col = col_fun,
                      top_annotation = c(col_order), 
                      #column_title_gp = gpar(fontsize = 12),
                      use_raster = FALSE)

ht_NO_scale

png('kmeans12_heatmap_sig_DEG_NoScaled.jpg', width = 12, height = 13, units = 'in', res = 1200)
draw(ht_NO_scale)
dev.off()

#Combine plot

col_class_com <-  HeatmapAnnotation(class = c(rep(c("SynCom", "Actinobacteria", "Flavobacteriia", 
                                                    "Bacilli", "Alphaproteobacteria", "Alphaproteobacteria", "Alphaproteobacteria", 
                                                    "Betaproteobacteria", "Gammaproteobacteria", "Gammaproteobacteria"), each = 1),
                                              rep(c("SynCom", "Actinobacteria", "Flavobacteriia", 
                                                    "Bacilli", "Alphaproteobacteria", "Alphaproteobacteria", "Alphaproteobacteria", 
                                                    "Betaproteobacteria", "Gammaproteobacteria", "Gammaproteobacteria"), each = 1)),
                                    col = list(class = c("SynCom" = "#999999",
                                                         "Actinobacteria" = "#FF9999",
                                                         "Flavobacteriia" = "#99FFFF",
                                                         "Bacilli" = "#FFCC66",
                                                         "Alphaproteobacteria" = "#66FF00",
                                                         "Alphaproteobacteria" = "#66FF00",
                                                         "Alphaproteobacteria" = "#66FF00",
                                                         "Betaproteobacteria" = "#00C00C",
                                                         "Gammaproteobacteria" = "#006006",
                                                         "Gammaproteobacteria" = "#006006" )),
                                    gp = gpar(col = "black"), # Border color
                                    annotation_legend_param = list(
                                      class = list(
                                        title_gp = gpar(fontsize = 16,  fontface = "bold"),  # Title font size
                                        nrow = 3,
                                        labels_gp = gpar(fontsize = 14) # Labels font size
                                      )
                                    )
)


ht_scale_combine <- Heatmap(matrix = scale_deganno_logFC_matrix %>% select(contains('_')),
                            heatmap_width = unit(20, "cm"),
                            name = 'Z-score',
                            cluster_columns = F,
                            #show_column_dend = F,
                            cluster_column_slices = F,
                            row_gap = unit(4, "mm"),
                            row_split = scale_deganno_logFC_matrix$cl,
                            show_column_names = F,
                            column_split = factor(rep(c("Shoot - Live bacteria", "Shoot - Heat killed bacteria"), c(10, 10)),
                                                  levels = c("Shoot - Live bacteria", "Shoot - Heat killed bacteria")),
                            column_order = 1 : 20,
                            column_gap = unit(2, 'mm'),
                            col = col_fun_1,
                            top_annotation = c(col_class_com), 
                            column_title_gp = gpar(fontsize = 18,  fontface = "bold"), #font size
                            row_title_gp = gpar(fontsize = 16, fontface = "bold"),
                            use_raster = FALSE,
                            heatmap_legend_param = list(
                              at = -3:3,  # Set custom legend breaks
                              labels = as.character(-3:3),
                              title_gp = gpar(fontsize = 16,  fontface = "bold"), # Increase legend title font size
                              labels_gp = gpar(fontsize = 14) 
                            )) +
  Heatmap(sigMat,
          heatmap_width = unit(14, "cm"),
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 14),
          column_names_rot = 45,
          heatmap_legend_param = list(
            title = 'DEGs',
            title_gp = gpar(fontsize = 16,  fontface = "bold"), # Legend title font size
            labels_gp = gpar(fontsize = 14) # Legend label font size
          ),
          cluster_columns = FALSE,
          use_raster = FALSE)


ht_scale_combine


#Export jpeg
png('kmeans10_heatmap_sig_DEG_Scaled_combineHeatplot_shoot_v2.jpg', width = 17, height = 14, units = 'in', res = 600)
draw(ht_scale_combine)
dev.off()

png('kmeans10_heatmap_sig_DEG_Scaled_combineHeatplot_shoot_test.jpg', width = 17, height = 12, units = 'in', res = 600)
draw(ht_scale_combine)
dev.off()

png('kmeans10_heatmap_sig_DEG_Scaled_combineHeatplot_shoot_test2.jpg', width = 17, height = 13, units = 'in', res = 600)
draw(ht_scale_combine)
dev.off()

png('kmeans10_heatmap_sig_DEG_Scaled_combineHeatplot_shoot_test3.jpg', width = 17, height = 14, units = 'in', res = 600)
draw(ht_scale_combine)
dev.off()


# Exporting to a PDF file
pdf('kmeans10_heatmap_sig_DEG_Scaled_combineHeatplot_shoot.pdf', width = 16, height = 13)
draw(ht_scale_combine)
dev.off()






#new hetmap layout - Feb 16, 2025
Heatmap(matrix = scale_deganno_logFC_matrix %>% select(contains('_')),
        heatmap_width = unit(20, "cm"),
        name = 'Z-score',
        cluster_columns = F,
        #show_column_dend = F,
        cluster_column_slices = F,
        row_gap = unit(4, "mm"),
        row_split = scale_deganno_logFC_matrix$cl,
        show_column_names = F,
        column_split = factor(rep(c("Shoot - Live bacteria", "Shoot - Heat killed bacteria"), c(10, 10)),
                              levels = c("Shoot - Live bacteria", "Shoot - Heat killed bacteria")),
        column_order = 1 : 20,
        column_gap = unit(2, 'mm'),
        col = col_fun_1,
        top_annotation = c(col_class_com), 
        column_title_gp = gpar(fontsize = 18,  fontface = "bold"), #font size
        row_title_gp = gpar(fontsize = 16, fontface = "bold"),
        use_raster = FALSE,
        heatmap_legend_param = list(
          at = -3:3,  # Set custom legend breaks
          labels = as.character(-3:3),
          direction = "horizontal",
          title_position = "lefttop",
          title_gp = gpar(fontsize = 16,  fontface = "bold"), # Increase legend title font size
          labels_gp = gpar(fontsize = 14) 
        )) +
  Heatmap(sigMat,
          heatmap_width = unit(14, "cm"),
          #border = TRUE,
          #border_gp = gpar(col = "black"),
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          show_column_names = F,
          column_title = NULL,  # Remove column titles
          #column_names_gp = gpar(fontsize = 14),
          #column_names_rot = 45,
          #column_split = factor(rep(1:20, each = ncol(sigMat) / 20)),  # Create 20 column groups
          column_split = factor(rep(c("Shoot - Live bacteria", "Shoot - Heat killed bacteria"), c(10, 10)),
                                levels = c("Shoot - Live bacteria", "Shoot - Heat killed bacteria")),
          column_gap = unit(1.25, 'mm'),
          top_annotation = c(col_class_com), 
          heatmap_legend_param = list(
            title = 'DEGs',
            title_gp = gpar(fontsize = 16,  fontface = "bold"),
            title_position = "lefttop",
            labels_gp = gpar(fontsize = 14), # Legend label font size
            nrow = 1
          ),
          cluster_columns = FALSE,
          use_raster = FALSE)

png('kmeans10_heatmap_sig_DEG_Scaled_combineHeatplot_shoot_New_Layout_021625_v2.jpg', width = 17, height = 14, units = 'in', res = 600)
draw(p1)
dev.off()


png('kmeans10_heatmap_sig_DEG_Scaled_combineHeatplot_shoot_New_Layout_021625_v3.jpg', width = 17, height = 14, units = 'in', res = 600)
draw(p1, 
     heatmap_legend_side = "bottom",  # Puts DEG and Z-score legend horizontally at the bottom
     annotation_legend_side = "bottom")  # Moves class annotation legend to the bottom

dev.off()

