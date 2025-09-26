##set the path to the working directory #miseq
cd /your working directory/
#activate the qiime2
conda activate qiime2-amplicon-2024.10

#import the data and mapping
qiime tools import \
    --type MultiplexedPairedEndBarcodeInSequence \
    --input-path seq \ ##the raw files of Read 1 and Read 2
    --output-path multiplexed-seqs.qza

#demultiplexing
 qiime cutadapt demux-paired   \  
    --i-seqs multiplexed-seqs.qza \    
    --m-forward-barcodes-file metadata.txt \   
    --m-forward-barcodes-column BarcodeSequence_f \    
    --m-reverse-barcodes-file metadata.txt   \  ##using the reverse complement of the reverse barcode
    --m-reverse-barcodes-column BarcodeSequence_r \    ##using the reverse complement of the reverse barcode
    --o-per-sample-sequences 2demultiplexed.qza   \  
    --o-untrimmed-sequences 2unmapped.qza     \
    --p-anchor-forward-barcode \
    --p-anchor-reverse-barcode \
    --p-mixed-orientation \
    --p-error-rate 0 \
    --p-cores 5 \
    --verbose

qiime demux summarize \
	--i-data 2demultiplexed.qza \
	--o-visualization 2demultiplexed.qzv

qiime tools export \
	--input-path 2demultiplexed.qza \
	--output-path 2demultiplexed-seqs-trimmed ## output of the demultiplexed files with assigned barcodes

### Removing the PCR primers out of the reads
echo "Step3: Removing primer sequences..."
qiime cutadapt trim-paired \
	--i-demultiplexed-sequences 2demultiplexed.qza \
	--p-front-f AACMGGATTAGATACCC \ 
	--p-front-r ACGTCATCCCCACCTTCC \ 
	--p-cores 5 \
	--p-error-rate 0.25 \
	--o-trimmed-sequences 3demultiplexed-seqs-trimmed.qza \
	--p-discard-untrimmed \
	--p-overlap 10 \
    --verbose

    #--p-front-f AACMGGATTAGATACCC \ #799F V5
    #--p-front-r ACGTCATCCCCACCTTCC \ # 1192R V7

qiime tools export \
	--input-path 3demultiplexed-seqs-trimmed.qza \
	--output-path 3demultiplexed-seqs-trimmed  ## output of the demultiplexed files with assigned barcodes without primer sequences 

qiime demux summarize   \
--i-data  3demultiplexed-seqs-trimmed.qza   \
--o-visualization 3trimmed-seqs-summary.qzv

##Rename all the sequence file to make them shorter and clearer for later use
dos2unix rename_fastqs.sh
chmod +x rename_fastqs.sh ## make the script executable
./rename_fastqs.sh

##merge the demultiplexed paired end sequences using flash2
#using run_flash2_batch.sh
chmod +x run_flash2_batch.sh ## make the script executable
./run_flash2_batch.sh
#if it cannot work, likely due to the failed format conversion from windows to Unix, try using:
#convert the file from DOS format to Unix format and run again
dos2unix run_flash2_batch.sh
chmod +x run_flash2_batch.sh ## make the script executable
./run_flash2_batch.sh
#INPUT_DIR="3demultiplexed-seqs-trimmed"
#OUTPUT_PARENT_DIR="flash2_results

##run usearch to filter out low quality sequences ##to organize all the merged sequences in the same folder
dos2unix run_usearch_filter.sh
chmod +x run_usearch_filter.sh ## make the script executable
./run_usearch_filter.sh

## merge all samples together
cd /your working directory/FLASH2/flash2_results/
echo "sample-id,absolute-filepath,direction" > /merge-manifest
find /usearch/*.extendedFrags_filter.fastq|while read line;
do
	name=`basename $line .extendedFrags_filter.fastq`
	echo "$name,$line,forward" >> merge-manifest
done

##proceed with rbec rather than DADA2 (rbec and DADA2 serve the same function)
##rbec shall be run in R
##run the script "rbec.R" in the folder /Script/

