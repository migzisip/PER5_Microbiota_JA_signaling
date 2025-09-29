# -------------------------------------------------------------
# ComplexHeatmap visualization for JA matrices (e.g., ja_match_counts.csv)
# - Columns (genes): fixed order following JA biosynthesis/signaling/catabolism
# - Column labels: "GENE - pathway"  (ASCII to avoid font/device issues)
# - Column splits: pathway groups (in biological order)
# - Rows (samples): fixed order 01_ ... 16_ based on numeric prefix
# - Uses ORIGINAL values (no log transform)
# - Saves as SVG
# -------------------------------------------------------------

suppressPackageStartupMessages({
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is not installed. Install via: BiocManager::install('ComplexHeatmap')")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is not installed. Install via: install.packages('circlize')")
  }
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ComplexHeatmap); library(circlize); library(grid)
})

# -----------------------------
# Config
# -----------------------------
infile <- "ja_match_counts.csv"    # Change to your matrix file if needed
out_svg <- "ja_counts_heatmap_ordered.svg"

# Desired biological order of pathway blocks (must match your 'pathway' strings)
# NOTE: use ASCII arrow "->" to avoid Windows SVG encoding warnings
pathway_levels <- c(
  "JA biosynthesis - lipoxygenase",
  "JA biosynthesis - allene oxide synthase",
  "JA biosynthesis - allene oxide cyclase",
  "JA biosynthesis - OPDA reductase",
  "JA biosynthesis - OPC-8:0 CoA ligase",
  "JA biosynthesis - peroxisomal beta-oxidation",
  "JA conjugation - JA->JA-Ile",
  "JA signaling - receptor (F-box)",
  "JA signaling - JAZ repressor",
  "JA signaling - bHLH TF",
  "JA-Ile catabolism - CYP94"
)

# Desired gene ordering within the full pathway (left-to-right)
desired_gene_order <- c(
  # LOX family
  "LOX2","LOX3","LOX4","LOX6",
  # AOS / AOC / OPR3 / OPCL1
  "AOS","AOC1","AOC2","AOC3","AOC4","OPR3","OPCL1",
  # Peroxisomal beta-oxidation
  "ACX1","ACX2","ACX3","ACX4","ACX5","ACX6","MFP2","KAT2","KAT1","KAT5",
  # Conjugation
  "JAR1",
  # Signaling (COI1, JAZs, MYCs)
  "COI1","JAZ1","JAZ2","JAZ3","JAZ4","JAZ10","JAZ12","MYC2","MYC3","MYC4",
  # Catabolism
  "CYP94B3","CYP94C1","CYP94B1"
)

# -----------------------------
# Load data
# -----------------------------
df <- readr::read_csv(infile, show_col_types = FALSE)

meta_cols <- c("agi","gene","pathway","uniprot_acc")
present_meta <- intersect(meta_cols, names(df))
sample_cols <- setdiff(names(df), present_meta)
stopifnot(length(sample_cols) > 0)

# Normalize pathway label variants (convert any arrows to ASCII ->)
if ("pathway" %in% names(df)) {
  df$pathway <- df$pathway |>
    # normalize various arrow encodings to ASCII
    str_replace_all("\u2192", "->") |>         # →  -> 
    str_replace_all("JA\\?JA-Ile", "JA->JA-Ile")
  df$pathway <- factor(df$pathway, levels = pathway_levels, ordered = TRUE)
}

# Enforce gene order
if ("gene" %in% names(df)) {
  df$gene <- as.character(df$gene)
  df$.__gene_order__ <- factor(df$gene, levels = desired_gene_order, ordered = TRUE)
  df <- df %>% arrange(.__gene_order__)
  df$.__gene_order__ <- NULL
}

# Build matrix: genes x samples (in our enforced order)
mat <- as.matrix(df[, sample_cols])
mat[is.na(mat)] <- 0

# Transpose -> rows = samples, cols = genes
mat_t <- t(mat)

# Order rows (samples) by numeric prefix 01_ ... 16_
sample_names <- rownames(mat_t)
prefix_num <- suppressWarnings(as.integer(str_extract(sample_names, "^[0-9]+")))
ord_rows <- order(ifelse(is.na(prefix_num), Inf, prefix_num), sample_names)
mat_t <- mat_t[ord_rows, , drop = FALSE]

# Column labels: "GENE - pathway" (ASCII only to avoid device warnings)
col_labels <- if (all(c("gene","pathway") %in% names(df))) {
  paste0(df$gene, " - ", as.character(df$pathway))
} else if ("gene" %in% names(df)) {
  df$gene
} else if ("agi" %in% names(df)) {
  df$agi
} else {
  colnames(mat_t)
}
# sanitize any lingering Unicode arrows/dashes
col_labels <- col_labels |>
  str_replace_all("\u2192", "->") |>
  str_replace_all("\u2014", "-")   # em dash -> hyphen
colnames(mat_t) <- col_labels

# Column splits by pathway (in biological order)
col_split <- if ("pathway" %in% names(df)) df$pathway else NULL

# Color mapping using ORIGINAL values
rng <- range(mat_t, na.rm = TRUE)
# Protect against constant matrix
if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
  rng <- c(0, 1)
}
col_fun <- colorRamp2(c(rng[1], mean(rng), rng[2]), c("snow", "lavenderblush3", "plum3"))

# Build heatmap (no clustering, fixed biological order and sample order)
ht <- Heatmap(
  mat_t,
  name = "count",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 7.5),
  column_names_rot = 60,
  column_split = col_split,     # keep splits
  column_title = NULL,          # hide large titles at the top
  border = TRUE,
  heatmap_legend_param = list(title = "count")
)

# Save as SVG
svg(out_svg, width = 14, height = 8)
draw(ht)
dev.off()

message("Saved: ", out_svg)