
#GO analysis - GOseq

meanWs4per5 <- function(v) {
  
  require('magrittr')
  
  res <- v %>%
    split(rep(1 : 8, each = 3)) %>%
    sapply(mean, na.rm = TRUE)
  
  return(res)
}

#######################GO analysis############################
setwd('D:/Projects/PER5_RNA_seq/Root_DEGs/')


library('goseq')
library('GO.db')
library('foreach')
library('doMC')
library('KEGGAPI')
library('BioCycAPI')
library('magrittr')
library('dplyr')
library('tibble')
library('readr')
library('stringr')
library('enrichplot')

#registerDoMC(12)

load('athGO.RData')
load('athKEGG.RData')
load('athBioCyc.RData')

###############################cluster profiler#####################
library('org.At.tair.db')
library('clusterProfiler')
library('magrittr')
library('tidyverse')
library('RColorBrewer')
library('GOSemSim')
library('DOSE')
library('RColorBrewer')
library('enrichplot')
library('stringr')


kmeansRes_sig<- read_csv('kmeans8_sig.csv',
                         col_types = cols(Chromosome = col_character()))
#kmeansBkg <- read_csv('Root_All_DEGs.csv',
                      #col_types = cols(Chromosome = col_character()))

savepath <- 'D:/Projects/PER5_RNA_seq/Root_DEGs/GO_1//'

setwd(savepath)

prefix <- 'kmeans8'

for (i in kmeansRes_sig $cl %>% unique) {
  
  ## BP
  goBP <- enrichGO(gene = kmeansRes_sig  %>% filter(cl == i) %>% .$ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1) %>% unlist %>% unique,
                   OrgDb = 'org.At.tair.db',
                   keyType= 'TAIR',
                   ont = 'BP',
                   universe = keys(org.At.tair.db),
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1)
  
  
  goBPSim <- clusterProfiler::simplify(goBP,
                                       cutoff = 0.5,
                                       by = 'p.adjust',
                                       select_fun = min)
  ## check and plot
  write.csv(as.data.frame(goBPSim),
            paste0(prefix, '_cluster', i, '_cp_BP.csv') %>% file.path(savepath, .))
  
  ## KEGG
  kk2 <- enrichKEGG(gene = kmeansRes_sig  %>% filter(cl == i) %>% .$ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1) %>% unlist %>% unique,
                    organism = 'ath',
                    pvalueCutoff = 0.05)
  
  write.csv(as.data.frame(kk2),
            paste0(prefix, '_cluster', i, '_cp_KEGG.csv') %>% file.path(savepath, .))
}

#organ specific cluster may not have highly significant DEGs e.g. cluster 9 and thus filtered out
kall <- lapply(kmeansRes_sig $cl %>% unique %>% .[!(. %in% c(9, 10))], function(x) {
  
  eachG <- kmeansRes_sig  %>% filter(cl == x) %>% .$ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1) %>% unlist %>% unique
  
  return(eachG)
  
}) %>%
  set_names(kmeansRes_sig $cl %>% unique %>% .[!(. %in% c(9, 10))] %>% paste0('cluster', .))

kallGOBP <- compareCluster(geneCluster = kall,
                           fun = 'enrichGO',
                           OrgDb = 'org.At.tair.db',
                           keyType= 'TAIR',
                           ont = 'BP',
                           universe = keys(org.At.tair.db),
                           pAdjustMethod = 'BH',
                           pvalueCutoff=0.01,
                           qvalueCutoff=0.1)

kallGOBPSim <- clusterProfiler::simplify(kallGOBP,
                                         cutoff = 0.5,
                                         by = 'p.adjust',
                                         select_fun = min) 
#dotplot(kallGOBP, showCategory = 10)

# dotplot(kallGOBPSim, showCategory = 5) +
#   theme(axis.text.y = element_text(size=12)) +
#   scale_y_discrete(labels = function(y) str_wrap(y, width = 60))
# 
dotplot(kallGOBPSim, showCategory = 5) +
  scale_y_discrete(labels=function(x) str_wrap(x, width=100)) +
  theme(axis.text.y = element_text(size=16, face = 'bold'),
        axis.text.x = element_text(size=15),
        legend.text=element_text(size= 15),
        legend.title = element_text(size = 16, face="bold"))



#Final plot

options(enrichplot.colours = c("blue", "red"))

# Modify kallGOBPSim to add -log10(p.adjust)
kallGOBPSim@compareClusterResult <- kallGOBPSim@compareClusterResult %>%
  mutate(log10_FDR = -log10(p.adjust))  # Compute log-transformed p.adjust

# Plot with corrected color mapping
dotplot(kallGOBPSim, showCategory = 5, color = 'log10_FDR') +  # Use log_padj as color
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  #scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adjust)") +  # Custom color gradient
  theme(
    axis.text.y = element_text(size = 18, face = 'bold'),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.title = element_text(size = 18, face = 'bold'),
    legend.text = element_text(size = 17),
    legend.title = element_text(size = 18, face = "bold"),
    panel.border = element_rect(size = 2, color = "black"),  # Thicker border
    panel.grid.major = element_line(color = "black", size = 0.1)  # Black grid lines
    #panel.grid.minor = element_line(color = "black", size = 0.1)   # Minor grid lines visible
  )


#dotplot(kallGOBPSim, showCategory = 10) + coord_flip() +  theme(axis.text.x = element_text(angle = 90, size=10))

ggsave('kmeans8_Ws4per5_cp_BP_dotplot_top5_root.tiff', width = 12, height = 14, dpi = 300)
ggsave('kmeans8_Ws4per5_cp_BP_dotplot_top5_root.pdf', width = 12, height = 14, dpi = 300)


# ggsave('kmeans8_Ws4per5_cp_BP_dotplot_roots.jpg', width = 15, height = 15)
# ggsave('kmeans8_Ws4per5_cp_BP_dotplot.pdf', width = 15, height = 15)

kallGOBP %>%
  as.data.frame %>%
  write_csv('kmeans8_Ws4per5_cp_BP.csv')

save(kallGOBP, file = 'kmeans11_Ws4per5_cp_BP.RData')

# https://rdrr.io/github/GuangchuangYu/enrichplot/man/emapplot.html
AtGO <- godata('org.At.tair.db', ont="BP")

kallGOBP_pairwise<-pairwise_termsim(kallGOBP, method = "JC", semData = AtGO, showCategory = 200)

#display.brewer.all()

#colour_pie<- brewer.pal(9,"Set1")

emapplot(kallGOBP_pairwise,
         showCategory = 20,
         pie='count',
         pie_scale=1.5,
         col= colour_pie,
         cex_label_category=1.0,
         repel=TRUE,
         layout='nicely')
#ggsave('kmeans11_Ws4per5_cp_BP_network_20cat.jpg', width = 15, height = 12)
ggsave('kmeans8_Ws4per5_cp_BP_network_20cat_Root.pdf', width = 15, height = 12)


#Selected cluster of interest

selected_GO <- c(
  'jasmonic acid mediated signaling pathway',
  'regulation of jasmonic acid mediated signaling pathway',
  'cellular response to jasmonic acid stimulus',
  'cellular response to fatty acid'
)


emapplot(
  kallGOBP_pairwise,
  showCategory = selected_GO, 
  cex_label_category = 2,
  pie.params = list(pie = 'count', pie_scale = 1.5),
  col = colour_pie,
  #cex.params = list(category_label = 1.0),
  repel = TRUE,
  layout.params = list(layout = 'nicely')) +
  theme(
    text = element_text(size = 14), # Adjusts the font size of all plot text
    legend.text = element_text(size = 14), # Adjusts legend text size
    legend.title = element_text(size = 16) # Adjusts legend title size
  )




#take gene symbol
kallGOBPSim_r <- setReadable(kallGOBPSim, 'org.At.tair.db')
kallGOBP_r <- setReadable(kallGOBP, 'org.At.tair.db')

#Layout KK
cnetplot(kallGOBP_r, showCategory = selected_GO,  categorySize=1,  
         vertex.label.font=10, vertex.label.cex = 2, cex_label_category = 3, cex_label_gene = 2,
         layout="kk", repel = TRUE) +
  theme(
    text = element_text(size = 14), # Adjusts the font size of all plot text
    legend.text = element_text(size = 12), # Adjusts legend text size
    legend.title = element_text(size = 14) # Adjusts legend title size
  )



000000# Layout dh
cnetplot(kallGOBP_r, 
         showCategory = selected_GO,  
         categorySize = 1,  
         vertex.label.font = 10, 
         vertex.label.cex = 1.8, 
         cex_label_category = 2.5, 
         cex_label_gene = 2,  # Hide gene labels for clarity
         layout = "dh",     # Try "tree" or "circle" for non-overlapping layout
         repel = TRUE) +
  theme(
    text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

ggsave('kmeans8_Ws4per5_cp_BP_CNETPlot_Root.pdf', width = 10, height = )





#Plot KEGG Pathway

kallKEGG <- compareCluster(geneCluster = kall,
                           fun = 'enrichKEGG',
                           organism = 'ath',
                           pvalueCutoff = 0.05)
dotplot(kallKEGG)











# ## top 5
# interesGO <- list(JA = c('GO:0009753','GO:0009694'),
#                   ET = c('GO:0009755'),
#                   SA = c('GO:0009696', 'GO:0009863'),
#                   defense = c('GO:0006955', 'GO:0098542', 'GO:0009814'),
#                   hypoxia = c('GO:0071456', 'GO:0070483'))


# topGONum <- 20
# 
# anno <- read_csv('U:/10 RNAseq/Ws4 per5 5686/heatmap/Ensembl_ath_Anno.csv',
#                  col_types = cols(Chromosome = col_character())) %>%
#   mutate_all(list(~replace(., is.na(.), ''))) %>%
#   mutate(GeneID = strsplit(ID, split = '.', fixed = TRUE) %>%
#            sapply('[[', 1) %>%
#            unlist) %>%
#   mutate(Gene = if_else(nchar(Gene) == 0, GeneID, Gene)) %>%
#   dplyr::select(ID, GeneID, Gene, Description) %>%
#   dplyr::slice(which(!duplicated(.)))
# 
# scaleC2 <- rawC %>%
#   select(contains('_')) %>%
#   t %>%
#   scale %>%
#   t %>%
#   as_tibble %>%
#   bind_cols(rawC %>% select(ID, cl)) %>%
#   mutate(GeneID = ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[', 1))
# 
# cpBP <- fortify(kallGOBP,showCategory = topGONum) %>%
#   as_tibble %>%
#   mutate(geneName = sapply(geneID, function(x) {
#     strsplit(x, split = '/', fixed = TRUE) %>%
#       unlist %>%
#       tibble(GeneID = .) %>%
#       inner_join(anno) %>%
#       .$Gene %>%
#       paste(collapse = '/') %>%
#       str_replace('C/VIF2', 'C-VIF2') ## replace genes with '/'
#   }))
# 
# for (i in seq_along(interesGO)) {
#   
#   interesGene <- cpBP %>%
#     filter(ID %in% interesGO[[i]]) %>%
#     .$geneID %>%
#     strsplit(split = '/', fixed = TRUE) %>%
#     unlist %>%
#     unique
#   
#   interesMat <- scaleC %>%
#     dplyr::filter(GeneID %in% interesGene) %>%
#     dplyr::filter(!(cl %in% c(9:10))) %>%
#     inner_join(anno)
#   
#   matcol <- colorRamp2(seq(min(scaleC %>% select(contains('_'))), max(scaleC %>% select(contains('_'))), length = 100), colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100))
#   
#   dim(interesMat) %>% print
#   
#   defenseGene <- c('AT2G19190',
#                    'AT1G14550',
#                    'AT1G18570',
#                    'AT2G30750',
#                    'AT1G73805',
#                    'AT4G28460',
#                    'AT4G37290',
#                    'AT3G48090',
#                    'AT3G07040',
#                    'AT3G50950',
#                    'AT4G11170',
#                    'AT5G41540')
#   defenseGeneIdx <- match(defenseGene, interesMat$GeneID) %>% .[!is.na(.)]
#   defenseAnno <- rowAnnotation(foo = anno_mark(at = defenseGeneIdx, labels = interesMat$Gene[defenseGeneIdx]))
#   
#   ht_list <- Heatmap(matrix = interesMat %>%
#                        select(contains('_')) %>%
#                        apply(1, meanWs4per5) %>%
#                        t,
#                      name = 'Scaled Counts',
#                      row_order = order(interesMat$cl) %>% rev,
#                      row_split = interesMat$cl,
#                      row_gap = unit(2, "mm"),
#                      column_order = 1 : 10,
#                      column_split = rep(c('Mock/HKSynCom', 'Non-sup', 'Sup'), c(6, 2, 2)),
#                      show_column_names = FALSE,
#                      col = matcol,
#                      use_raster = FALSE,
#                      right_annotation = defenseAnno)
# }
