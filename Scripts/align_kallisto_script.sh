#!/bin/bash

kallisto_index=/RAID1/working/R324/Per5_RNA_seq/Data_Analysis/CDNA_reference_ATH/ATH.kallisto.idx
OutDir=/RAID1/working/R324/Per5_RNA_seq/Data_Analysis/02.Alignment
input_dir=/RAID1/working/R324/Per5_RNA_seq/Data_Analysis/01.Raw_data/trimmed_fastq

for R1 in ${input_dir}/*_R1_fp.fq.gz

do
# Extract the base name of the sample (removing _R1_fp.fq.gz)
  sample=$(basename "$R1" _R1_fp.fq.gz)
  
  # Define the corresponding R2 file
  R2="${input_dir}/${sample}_R2_fp.fq.gz"
   
	echo "Mapping ${sample} to ATH cDNA"
	kallisto quant \
		--threads 32 \
		--index ${kallisto_index} \
		--output-dir ${OutDir}/${sample}_kallisto \
		$R1 $R2 &> ${OutDir}/${sample}_kallisto.log
done