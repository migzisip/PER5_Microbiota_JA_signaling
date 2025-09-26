#!/bin/bash

# Set the working directory where the FASTQ files are
INPUT_DIR="flash2_results"
OUTPUT_DIR="flash2_results/usearch/"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through all the extendedFrags FASTQ files
for FILE in ${INPUT_DIR}/*/*.extendedFrags.fastq
do
    # Extract the base sample name (without path and extension)
    BASENAME=$(basename "$FILE" .extendedFrags.fastq)
    
    # Run usearch filter
    usearch -fastq_filter "$FILE" \
            -fastqout "${OUTPUT_DIR}/${BASENAME}.extendedFrags.filtered.fastq" \
            -fastq_maxns 0 \
            -fastq_minlen 200

    echo "Processed: $BASENAME"
done

echo "All 96 samples processed!"