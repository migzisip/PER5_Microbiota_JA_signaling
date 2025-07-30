# Shoot Transcriptome analysis

#Load the required libraries
library('tximport')
library('rhdf5')
library('magrittr')
library('DESeq2')
library('tidyverse')

#=================User define function==================#
checkZeros <- function(v, threshold) {
  res <- ifelse(sum(v == 0) > threshold, FALSE, TRUE)
  return(res)
}

checkMinExpression  <- function(v, threshold) {
  
  require('magrittr')
  
  res <- v %>%
    split(rep(1 : 8, each = 3)) %>%
    sapply(checkZeros, threshold) %>%
    all
  
  return(res)
}
#======================================================#

#Load Arabidopsis annotation file (Csv format)
anno <- read_csv('D:/Projects/PER5_RNA_seq/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), '')))


#==========Load the kallisto alignment data============#
annoSample <- read_csv("D:/Projects/PER5_RNA_seq/Shoot_DEGs/list_samples.csv")


wd <- "D:/Projects/PER5_RNA_seq/02.Alignment"
slabel <- annoSample$Sample %>%
  paste0("_kallisto")

files <- file.path(wd, slabel, "abundance.h5")
names(files) <- annoSample$Sample
kres <- tximport(files, type = "kallisto", txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#===============================Normalization=========================#
setwd('D:/Projects/PER5_RNA_seq/Shoot_DEGs/')

## sampleTable
sampleInfo <- annoSample[match(colnames(kres$counts), annoSample$Sample), ] %>%
  as.data.frame
rownames(sampleInfo) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres, sampleInfo, ~ Condition)


## remove 0|0|x, 0|0|0
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkMinExpression, 1) %>%
  degres[., ]

degres <- DESeq(degres)

## count transformation
rld <- rlog(degres)
ntd <- normTransform(degres)

##===============Batch effects===========###


##===============Batch effects===========###
library('sva')
library('ggplot2')

dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}
mod <- model.matrix(~ Condition, colData(degres))
mod0 <- model.matrix(~ 1, colData(degres))

## ## manual detect surrogate variance
## svnum <- 4
## svseq <- svaseq(dat, mod, mod0, n.sv = svnum)

## auto detect sv
svobj <- sva(dat, mod, mod0)
svnum <- svobj$sv %>% ncol

SVA_plot <- svobj$sv %>%
  set_colnames(paste0('sv', seq_len(svnum))) %>%
  as_tibble %>%
  gather(key = 'sv', value = 'value') %>%
  mutate(Condition = colData(degres) %>%
           .$Condition %>%
           rep(svnum) %>%
           as.character,
         sample = rep(colnames(degres), svnum)) %>%
  mutate(group = paste(sv, Condition, sep = '_')) %>%
  ggplot(aes(sample, value, colour = sv, group = group)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))

SVA_plot

ggsave(SVA_plot, file = 'auto_sv.jpg', width = 12, height = 8, dpi = 300)
ggsave(SVA_plot, file = 'auto_sv.pdf', width = 12, height = 8, dpi = 300)

#======================================================================#

#================================DEGs Mock vs Conditions==================================#
degres$sv1 <- svobj$sv[, 1]
degres$sv2 <- svobj$sv[, 2]


design(degres) <- ~sv1 + ~sv2 + Condition  # Only include condition factor

## Perform Differential Expression Analysis with DESeq
degres <- DESeq(degres)


## Define Contrast for SynCom vs Mock
## root
## shoot
cond <- list(
  c("per5_L_M", "Ws4_L_M"),
  c("Ws4_L_SC", "Ws4_L_M"),
  c("per5_L_SC", "Ws4_L_SC"),
  c("per5_L_SC", "per5_L_M"))


resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(contrast = c('Condition', x)) %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs = list(~paste0(paste(x, collapse = '_vs_'), '_', .)))
                 }) %>%
  bind_cols

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), counts(degres, normalize = TRUE), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene:Description, Ws4_L_M1 : per5_L_SC_vs_per5_L_M_log2FoldChange) %>%  # Select relevant columns
  arrange(per5_L_SC_vs_per5_L_M_log2FoldChange)  

#Export DEG tables
res %>%
  filter(per5_L_SC_vs_per5_L_M_padj <= 0.05,
         per5_L_SC_vs_per5_L_M_log2FoldChange %>% abs > log2(1.5)) %>%
  write_csv("per5_L_M_SC_DEGs.csv")

res %>%
  filter(
    Ws4_L_SC_vs_Ws4_L_M_padj <= 0.05,
    Ws4_L_SC_vs_Ws4_L_M_log2FoldChange %>% abs() > log2(1.5)
  ) %>%
  write_csv("Ws4_L_M_SC_DEGs.csv")

write_csv(res, "Shoot_All_DEGs.csv")
#==================================================================#


#===============================PCA===============================#
library('ggrepel')
library('ggplot2')
library('RColorBrewer')
library('limma')
library('sva')
library(dplyr)

#~~~~~~~~~~~~~~~~~~#
## remove low count
dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}

group <- sampleInfo$Condition
design <- model.matrix(~ group)
rldData <- dat %>%
  removeBatchEffect(covariates = svobj$sv,
                    design = design)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###

#PCA ALL
allIdx <- 1:12
selectIdx <- allIdx

pca <- prcomp(t(rldData[, selectIdx]))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[, 1]
pca2 <- pca$x[, 2]
pcaData <- tibble(
  PC1 = pca1,
  PC2 = pca2,
  Condition = colData(rld)$Condition[selectIdx],
  SynCom = colData(rld)$SynCom[selectIdx],
  Tissue = colData(rld)$Tissue[selectIdx],
  Line = colData(rld)$Line[selectIdx],
  Sample = colData(rld)$Sample[selectIdx],
)


#Plot PCA

ggplot(
  pcaData,
  aes(x = PC1, y = PC2, colour = Condition, label = Sample)
) +
  geom_point(aes(shape = SynCom), size = 4) +
  coord_fixed(ratio = 1) +
  scale_colour_manual(values=c("#00cd00","#00cd00","#6BB100","#6BB100"))+
  #scale_colour_manual(values=c("#00cd00","#00cd00","#6BB100","#6BB100"))+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggforce::geom_mark_ellipse(aes(colour=Condition),label.fontsize=5,label.width=NULL) +
  #geom_text_repel(force=20) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.text.align = 0,
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 16, face = "bold"),
        legend.text=element_text(size= 15),
        legend.title = element_text(size = 17, face = "bold"))


ggsave("k_pca_sva_shoot.pdf")
ggsave("k_pca_sva_shoot.jpg")







