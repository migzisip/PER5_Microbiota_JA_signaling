#!/bin/bash

# Set your working directory
WORK_DIR="/mnt/c/Users/user/Desktop/symcom/miseq/FLASH2/3demultiplexed-seqs-trimmed/"

# Change to that directory
cd "$WORK_DIR"

# Loop through wsper groups (1 to 48 or however many you have)
for i in {1..96}
do
    # Find the R1 and R2 files for each wsperX
    R1_FILE=$(ls wsper${i}_*_R1_001.fastq.gz 2>/dev/null)
    R2_FILE=$(ls wsper${i}_*_R2_001.fastq.gz 2>/dev/null)

    # If R1 file exists, rename it
    if [[ -f "$R1_FILE" ]]; then
        mv "$R1_FILE" wsper${i}_R1.fastq.gz
        echo "Renamed $R1_FILE to wsper${i}_R1.fastq.gz"
    else
        echo "No R1 file found for wsper${i}"
    fi

    # If R2 file exists, rename it
    if [[ -f "$R2_FILE" ]]; then
        mv "$R2_FILE" wsper${i}_R2.fastq.gz
        echo "Renamed $R2_FILE to wsper${i}_R2.fastq.gz"
    else
        echo "No R2 file found for wsper${i}"
    fi
done

echo "Renaming complete!"
