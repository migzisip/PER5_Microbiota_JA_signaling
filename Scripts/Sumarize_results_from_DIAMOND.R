# build_ja_match_matrices.R
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tools)
})

# --------------------------
# Config
# --------------------------
base_dir      <- "D:/Loi/Per5 MS/diamond_results"
mapping_path  <- file.path(base_dir, "ja_uniprot_map.tsv")
pattern       <- "^[0-9]{2}_Root.*\\.tsv$"  # 01_Root..., 02_Root..., ..., 16_Root...

# Output files (match Python names)
out_all_hits          <- file.path(base_dir, "ja_hits_long_all.csv")
out_max_pident        <- file.path(base_dir, "ja_match_matrix_max_pident.csv")
out_counts            <- file.path(base_dir, "ja_match_counts.csv")
out_presence          <- file.path(base_dir, "ja_match_matrix_presence.csv")

# --------------------------
# Load mapping
# --------------------------
map_df <- read_tsv(mapping_path, show_col_types = FALSE) %>%
  select(agi, gene, pathway, uniprot_acc) %>%
  distinct()

# --------------------------
# Helper: read one TSV (BLAST-like, no header)
# --------------------------
read_one_hits <- function(fp) {
  # If your TSVs actually have a header, set col_names = TRUE and remove 'set_names' below.
  raw <- read_tsv(fp, col_names = FALSE, show_col_types = FALSE)
  stopifnot(ncol(raw) >= 12)
  raw <- raw %>%
    set_names(c("target","query","pident","length","mismatch","gapopen",
                "qstart","qend","sstart","send","evalue","bitscore")) %>%
    mutate(sample = file_path_sans_ext(basename(fp))) %>%
    # Extract AGI / gene / UniProt accession from the "target" field like:
    # "AT2G06050|OPR3|UniProt:Q9FUP0|pathway=JA"
    mutate(
      agi               = str_extract(target, "^[^|]+"),
      gene_from_target  = str_match(target, "^[^|]+\\|([^|]+)")[,2],
      uniprot_acc_target= str_match(target, "UniProt:([^|]+)")[,2]
    )
  # Ensure numeric pident (coerce if necessary)
  raw <- raw %>%
    mutate(
      pident   = suppressWarnings(as.numeric(pident)),
      length   = suppressWarnings(as.numeric(length)),
      mismatch = suppressWarnings(as.numeric(mismatch)),
      gapopen  = suppressWarnings(as.numeric(gapopen)),
      qstart   = suppressWarnings(as.numeric(qstart)),
      qend     = suppressWarnings(as.numeric(qend)),
      sstart   = suppressWarnings(as.numeric(sstart)),
      send     = suppressWarnings(as.numeric(send)),
      evalue   = suppressWarnings(as.numeric(evalue)),
      bitscore = suppressWarnings(as.numeric(bitscore))
    )
  raw
}

# --------------------------
# Collect & read all result files
# --------------------------
files <- list.files(base_dir, pattern = pattern, full.names = TRUE)
if (length(files) == 0L) stop("No files matched pattern in base_dir: ", base_dir)

hits <- files %>% map_dfr(read_one_hits)

# Write ALL hits (long format)
write_csv(hits, out_all_hits)

# --------------------------
# Summaries by AGI × sample
# --------------------------
# Max % identity
max_pident <- hits %>%
  group_by(agi, sample) %>%
  summarise(pident = max(pident, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = pident) %>%
  right_join(map_df, by = "agi") %>%   # keep all JA genes
  relocate(agi, gene, pathway, uniprot_acc) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)))

write_csv(max_pident, out_max_pident)

# Count of hits
counts <- hits %>%
  count(agi, sample, name = "n_hits") %>%
  pivot_wider(names_from = sample, values_from = n_hits) %>%
  right_join(map_df, by = "agi") %>%
  relocate(agi, gene, pathway, uniprot_acc) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)))

write_csv(counts, out_counts)

# Presence/absence (1 if any hit, else 0)
presence <- hits %>%
  count(agi, sample, name = "n_hits") %>%
  mutate(value = as.integer(n_hits > 0)) %>%
  select(-n_hits) %>%
  pivot_wider(names_from = sample, values_from = value) %>%
  right_join(map_df, by = "agi") %>%
  relocate(agi, gene, pathway, uniprot_acc) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)))

write_csv(presence, out_presence)

message("Done:")
message("  - All hits:            ", out_all_hits)
message("  - Max %id matrix:      ", out_max_pident)
message("  - Hit counts:          ", out_counts)
message("  - Presence/absence:    ", out_presence)
