
# Misregulation of the jasmonate signaling pathway leads to altered plant microbiota interaction and plant stress responses

This repository contains scripts and datasets for the 2025 bioRxiv study by Lu et al. 2025, bioRxiv  
https://www.biorxiv.org/content/10.1101/2025.03.29.646076v1.full




## 📂 Data

The root-associated bacterial genomes (FASTA) and annotation files (GFF) used for RNA-seq read mapping, quantification, and genome analysis were obtained from:  
🔗 [NCBI BioProject PRJNA297942](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA297942)



## Scripts  

---

## 🧪 Scripts for RNA-seq Analysis

This section contains scripts for:

- RNA-seq preprocessing and quality control
- Read alignment and quantification
- Differential gene expression analysis
- Functional enrichment and visualization

### 1.1 Experimental Design and Raw Data

Arabidopsis thaliana ecotypes Col-0 and Ws4, as well as mutant lines including **per5**, were used in transcriptomic experiments.  
Plants were germinated with synthetic bacterial community **At-16SC1** on half-strength MS medium for 14 days. Roots and shoots were harvested separately and RNA was extracted using the RNeasy Plant Mini Kit (Qiagen).

RNA-seq libraries were prepared by Novogene-Europe and sequenced on an Illumina platform (paired-end 150 bp), with an average depth of ~3 Gb per sample. Each biological condition had 3 replicates.  

Raw reads will be deposited in GEO under accession **[GSEXXXXX]** (to be updated upon release).

---

### 1.2 Read Pre-processing and Alignment

Quality Checking and Trimming: [pre_processing_Script.sh](Scripts/pre_processing_Script.sh)
Alignment: [align_kallisto_script.sh](Scripts/align_kallisto_script.sh)

Raw RNA-seq reads were processed using **fastp** (v0.23.2) with default parameters for paired-end reads. Adapter sequences and low-quality bases were trimmed.  

High-quality reads were pseudo-aligned to the *Arabidopsis thaliana* TAIR10.58 transcriptome using **kallisto** (Bray et al., 2016). On average, ~11.6 million paired-end reads were obtained per sample.

Transcript-level quantifications were imported using **tximport** (Soneson et al., 2016) into R for downstream analysis.

### 1.3 Normalization and Clustering

Roots: [Rscript_Kmeans_Clustering_Heatmap_analysis_PER5_SynCom16_Roots.R](Scripts/Rscript_Kmeans_Clustering_Heatmap_analysis_PER5_SynCom16_Roots.R)  
Shoots: [Rscript_Kmeans_Clustering_Heatmap_analysis_PER5_SynCom16_Shoots.R](Scripts/Rscript_Kmeans_Clustering_Heatmap_analysis_PER5_SynCom16_Shoots.R)

Low-abundance transcripts absent in at least two replicates were removed.  
Counts were log2-transformed and batch effects were corrected using surrogate variable analysis (SVA) and the `removeBatchEffect()` function in the **limma** package (Ritchie et al., 2015).

Normalized counts were scaled and z-transformed by transcript.  
Z-scores were used for **k-means clustering** (k = 8), with cluster number determined by SSE and AIC.

### Output:

---

### 1.4 Differential Expression and GO Enrichment

Differential expression was analyzed using **DESeq2** (Love et al., 2014). Key comparisons included:

1. per5 mock vs WT mock  
2. SynCom-treated WT vs WT mock  
3. SynCom-treated per5 vs SynCom-treated WT  
4. SynCom-treated per5 vs per5 mock  

Transcripts with adjusted p-value ≤ 0.05 and |log2FC| ≥ log2(1.5) were considered significant.

Roots: [Rscript_DEG_analysis_PER5_SynCom16_Roots.R](Scripts/Rscript_DEG_analysis_PER5_SynCom16_Roots.R)  
Shoots: [script_DEG_analysis_PER5_SynCom16_Shoots.R](Scripts/Rscript_DEG_analysis_PER5_SynCom16_Shoots.R)

**GO enrichment** was performed using **GOseq** (Young et al., 2010), accounting for transcript length bias.  
Significant Biological Process terms (adj. p < 0.05) were visualized with **clusterProfiler** (Yu et al., 2012).

Roots: [Rscript_GO_analysis_PER5_SynCom16_Roots.R](Scripts/Rscript_GO_analysis_PER5_SynCom16_Roots.R)  
Shoots: [Rscript_GO_analysis_PER5_SynCom16_Shoots.R](Scripts/Rscript_GO_analysis_PER5_SynCom16_Shoots.R)

---

## 🌿 Marker Gene Analysis for JA and SA Pathways

This script generates a heatmap of selected marker genes involved in **jasmonic acid (JA)** and **salicylic acid (SA)** signaling pathways based on RNA-seq data from Arabidopsis roots and shoots.  

Z-score normalized log₂ fold-change values for different genotypes and SynCom treatments are visualized to highlight hormone-responsive expression patterns. Marker genes were grouped by functional annotation, and plotted using the `ComplexHeatmap` package.

JA: 

 Root SA: [Rscript_SA_Marker_Genes_Root.R](Scripts/Rscript_SA_Marker_Genes_Root.R)  
Shoot SA: [Rscript_SA_Marker_Genes_Shoot.R](Scripts/Rscript_SA_Marker_Genes_Shoot.R)

📊 Output: [SA marker PER5 Root.pdf](<Data/SA marker PER5 Root.pdf>);   [SA marker PER5 Shoot.pdf](<Data/SA marker PER5 Shoot.pdf>)

*More details and usage instructions can be found in each script folder.*

---

## 📦 Dependencies

Key tools and packages used:

- `kallisto`
- `DESeq2`
- `tximport`
- `edgeR`
- `ggplot2`, `Complexheatmap`, `clusterProfiler` (for visualization and enrichment)
- `R` (≥4.2.0)
---

## 📘 Citation

If you use this code or data, please cite:

Lu et al., 2025. *Misregulation of the jasmonate signaling pathway leads to altered plant microbiota interaction and plant stress responses*. bioRxiv. https://doi.org/10.1101/2025.03.29.646076

---

## 🛠 Contact

For questions or feedback, please open an issue or contact the corresponding author listed in the preprint.  

## 📚 References
Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology, 34(5), 525–527. https://doi.org/10.1038/nbt.3519

Soneson, C., Love, M. I., & Robinson, M. D. (2016). Differential analyses for RNA-seq: Transcript-level estimates improve gene-level inferences. F1000Research, 4, 1521. https://doi.org/10.12688/f1000research.7563.2

Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47. https://doi.org/10.1093/nar/gkv007

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8

Young, M. D., Wakefield, M. J., Smyth, G. K., & Oshlack, A. (2010). Gene ontology analysis for RNA-seq: Accounting for selection bias. Genome Biology, 11(2), R14. https://doi.org/10.1186/gb-2010-11-2-r14

Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). clusterProfiler: An R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology, 16(5), 284–287. https://doi.org/10.1089/omi.2011.0118

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

Yates, A. D., et al. (2022). Ensembl 2022. Nucleic Acids Research, 50(D1), D988–D995. https://doi.org/10.1093/nar/gkab1049

Leek, J. T., Johnson, W. E., Parker, H. S., Jaffe, A. E., & Storey, J. D. (2012). The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics, 28(6), 882–883. https://doi.org/10.1093/bioinformatics/bts034