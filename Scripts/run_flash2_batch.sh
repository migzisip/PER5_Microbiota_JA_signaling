#!/bin/bash

# Input and output folders
INPUT_DIR="3demultiplexed-seqs-trimmed"
OUTPUT_PARENT_DIR="flash2_results"

# FLASH2 parameters
MAX_OVERLAP=300
MISMATCH_RATIO=0.25

# Create the parent output directory if it doesn't exist
mkdir -p $OUTPUT_PARENT_DIR

# Loop through all R1 files in the input directory
for R1 in ${INPUT_DIR}/*_R1.fastq.gz
do
    # Get the sample base name (remove path and _R1.fastq.gz)
    SAMPLE_NAME=$(basename ${R1} _R1.fastq.gz)

    # Define R2 filename based on the SAMPLE_NAME
    R2=${INPUT_DIR}/${SAMPLE_NAME}_R2.fastq.gz

    # Check if the matching R2 file exists
    if [[ ! -f "$R2" ]]; then
        echo "Warning: Missing R2 file for ${SAMPLE_NAME}. Skipping..."
        continue
    fi

    echo "Processing sample: $SAMPLE_NAME"

    # Create a separate output folder for each sample (optional)
    SAMPLE_OUTPUT_DIR=${OUTPUT_PARENT_DIR}/${SAMPLE_NAME}
    mkdir -p $SAMPLE_OUTPUT_DIR

    # Run FLASH2 for the sample
    flash2 $R1 $R2 \
        -M $MAX_OVERLAP \
        -x $MISMATCH_RATIO \
        -o ${SAMPLE_NAME} \
        -d $SAMPLE_OUTPUT_DIR

    echo "Finished processing $SAMPLE_NAME"
done

echo "All 96 samples processed!"
