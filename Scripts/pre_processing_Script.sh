#!/bin/bash

# Directories
RAW_DIR="/RAID1/working/R324/Per5_RNA_seq/Data_Analysis/01.Raw_data" # Directory with raw fastq files
QC_BEFORE_DIR="/RAID1/working/R324/Per5_RNA_seq/Data_Analysis/01.Raw_data/fastQC_Before_Trim"
TRIMMED_DIR="/RAID1/working/R324/Per5_RNA_seq/Data_Analysis/01.Raw_data/trimmed_fastq"
QC_AFTER_DIR="/RAID1/working/R324/Miguel/RNA_seq_data/2.Quality_Check/fastQC_After_Trim"
MULTIQC_DIR="/RAID1/working/R324/Per5_RNA_seq/Data_Analysis/01.Raw_data/multiQC"

# Ensure directories exist
mkdir -p "$QC_BEFORE_DIR" "$TRIMMED_DIR" "$QC_AFTER_DIR" "$MULTIQC_DIR"

# Step 1: Pre-trimming quality check with FastQC
echo "Running FastQC on raw FASTQ files..."
fastqc --threads 36 --nogroup "${RAW_DIR}"/*.fq.gz -o "$QC_BEFORE_DIR"

# Step 2: Trimming with Fastp
echo "Trimming FASTQ files with Fastp..."
for R1 in ${RAW_DIR}/*_1.fq.gz; do
  base_name=$(basename "$R1" _1.fq.gz)
  R2="${RAW_DIR}/${base_name}_2.fq.gz"

  echo "Processing sample: $base_name"

  fastp \
    --thread 8 \
    --overrepresentation_analysis \
    --correction \
    --html "${TRIMMED_DIR}/${base_name}.html" \
    --json "${TRIMMED_DIR}/${base_name}.json" \
    -i "$R1" -I "$R2" \
    -o "${TRIMMED_DIR}/${base_name}_R1_fp.fq.gz" -O "${TRIMMED_DIR}/${base_name}_R2_fp.fq.gz" \
    &> "${TRIMMED_DIR}/${base_name}_fastp.log"
done

# Step 3: Post-trimming quality check with FastQC
echo "Running FastQC on trimmed FASTQ files..."
fastqc --threads 36 --nogroup "${TRIMMED_DIR}"/*_fp.fq.gz -o "$QC_AFTER_DIR"

# Step 4: MultiQC for summary
echo "Running MultiQC for combined report..."
multiqc "$QC_BEFORE_DIR" "$QC_AFTER_DIR" "$TRIMMED_DIR" -o "$MULTIQC_DIR"

echo "Pipeline completed. Reports are in: $MULTIQC_DIR"

