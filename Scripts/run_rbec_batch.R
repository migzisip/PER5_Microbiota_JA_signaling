# run_rbec_batch.R

library("Rbec")
library("ShortRead")

input_dir <- "C:/Users/user/Desktop/symcom/miseq/FLASH2/flash2_results/usearch"
ref_fasta <- "C:/Users/user/Desktop/symcom/miseq/rbec/ref_corrected.fasta"
output_dir <- "C:/Users/user/Desktop/symcom/miseq/rbec/results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fastq_files <- list.files(input_dir, pattern = "\\.extendedFrags\\.filtered\\.fastq$", full.names = TRUE)
fastq_files <- fastq_files[!grepl("wsper12\\.extendedFrags\\.filtered\\.fastq$", fastq_files)]

for (fq in fastq_files) {
  sample_name <- tools::file_path_sans_ext(basename(fq))
  sample_name <- sub("\\.extendedFrags_filtered$", "", sample_name)
  out_path <- file.path(output_dir, sample_name)
  dir.create(out_path, showWarnings = FALSE)
  
  cat("Running Rbec on:", sample_name, "\n")
  
  Rbec(
    fastq = fq,
    reference = ref_fasta,
    outdir = out_path,
    threads = 4,
    sampling_size = 5000 )
}
