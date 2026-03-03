library(rstudioapi)
# Get the path of the current script
path <- rstudioapi::getActiveDocumentContext()$path
# Set the working directory to the folder containing that script
setwd(dirname(path))
# Verify the change
getwd()

library(phyloseq)
library(dplyr)
library(ggplot2)
library(biomformat)
library(ggtext)
library(forcats)
library(limma)
library(ggrepel)

## microbial profiling for per5 project - Color-matched Labels
# --- 1. Global Settings & Data Loading ---
setwd("C://Users//user//Desktop//per5 microbial profiling//LFC_LFC")
alpha <- 0.05
pseudo <- 1e-8

phylum_colors <- c(
  "Actinobacteria"      = "#F6AEB0", "Firmicutes"          = "#F69230",
  "Bacteroidetes"       = "#A4E0E2", "Alphaproteobacteria" = "#AFDFC9",
  "Betaproteobacteria"  = "#16946D", "Gammaproteobacteria" = "#708A63",
  "Other"               = "grey"
)

shape_values <- c(
  "Both Significant"           = 8, 
  "Treatment Significant Only"  = 15, 
  "Genotype Significant Only"  = 17, 
  "Non-significant"            = 16 
)

# Load Data
otu_df    <- read.table("rbec_absolute_no_input.txt", header = TRUE, sep = "\t", row.names = 1)
metadata  <- read.table("metadata_no_input.txt", header = TRUE, sep = "\t", row.names = 1)
taxonomy  <- read.table("taxonomy.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

physeq <- phyloseq(otu_table(otu_df, taxa_are_rows = TRUE), 
                   sample_data(metadata), 
                   tax_table(as.matrix(taxonomy)))

# --- 2. Data Processing ---

get_combined_stats <- function(ps, comp, target_treatment) {
  ps_sub <- prune_samples(sample_data(ps)$compartment == comp, ps)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  Y <- as(otu_table(ps_sub), "matrix")
  if (!taxa_are_rows(ps_sub)) Y <- t(Y)
  logY <- log2(Y + pseudo)
  
  meta_sub <- as(sample_data(ps_sub), "data.frame")
  meta_sub$treatment <- relevel(factor(meta_sub$treatment), ref = "mock")
  meta_sub$genotype  <- relevel(factor(meta_sub$genotype),  ref = "Ws4")
  
  design <- model.matrix(~ treatment + genotype, data = meta_sub)
  fit <- eBayes(lmFit(logY, design))
  
  treat_coef <- paste0("treatment", target_treatment)
  
  res <- data.frame(
    taxon = rownames(logY),
    LFC_Treatment = fit$coefficients[, treat_coef],
    LFC_Genotype  = fit$coefficients[, "genotypeper5"],
    q_treatment   = topTable(fit, coef = treat_coef, number = Inf, sort.by = "none")$adj.P.Val,
    q_genotype    = topTable(fit, coef = "genotypeper5", number = Inf, sort.by = "none")$adj.P.Val,
    Compartment   = comp,
    Treatment     = target_treatment
  )
  return(res)
}

tasks <- list(
  list(comp = "root",  treat = "DCB"),
  list(comp = "shoot", treat = "DCB"),
  list(comp = "root",  treat = "MeJA"),
  list(comp = "shoot", treat = "MeJA")
)

all_stats <- lapply(tasks, function(t) get_combined_stats(physeq, t$comp, t$treat)) %>%
  bind_rows() %>%
  mutate(
    Group = factor(paste0(Treatment, "_", Compartment), 
                   levels = c("DCB_root", "DCB_shoot", "MeJA_root", "MeJA_shoot")),
    Phylum = case_when(
      taxon %in% c("Microbacteriaceae", "Intrasporangiaceae", "Streptomycetaceae", "Mycobacteriaceae") ~ "Actinobacteria",
      taxon %in% c("Bacillaceae") ~ "Firmicutes",
      taxon %in% c("Flavobacteriaceae") ~ "Bacteroidetes",
      taxon %in% c("Bradyrhizobiaceae", "Rhizobiaceae", "Caulobacteraceae", "Phyllobacteriaceae", "Hyphomicrobiaceae") ~ "Alphaproteobacteria",
      taxon %in% c("Alcaligenaceae", "Comamonadaceae", "Oxalobacteraceae") ~ "Betaproteobacteria",
      taxon %in% c("Pseudomonadaceae", "Xanthomonadaceae") ~ "Gammaproteobacteria",
      TRUE ~ "Other"
    ),
    Sig_Status = case_when(
      q_treatment < alpha & q_genotype < alpha ~ "Both Significant",
      q_treatment < alpha ~ "Treatment Significant Only",
      q_genotype < alpha  ~ "Genotype Significant Only",
      TRUE ~ "Non-significant"
    )
  )

# --- 3. Plotting with Color-Matched Text ---

faceted_plot <- ggplot(all_stats, aes(x = LFC_Genotype, y = LFC_Treatment)) +
  # Quadrant backgrounds
  geom_rect(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "white", alpha = 0.02) +
  geom_rect(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "white", alpha = 0.02) +
  geom_rect(xmin = -Inf, xmax =0 , ymin = 0, ymax = Inf, fill = "gray60", alpha = 0.01) +
  geom_rect(xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "gray60", alpha = 0.01)+
  # Points
  geom_point(aes(color = Phylum, shape = Sig_Status), alpha = 0.7, size = 2.5) +
  # Labels - NOW WITH aes(color = Phylum)
  geom_text_repel(aes(label = taxon, color = Phylum), # <--- Color mapping added here
                  size = 2, 
                  fontface = "bold.italic", 
                  max.overlaps = 25, 
                  box.padding = 0.5,
                  segment.color = NA, # Makes the connecting line slightly transparent
                  show.legend = FALSE) + # Prevents "a" appearing in the legend
  # Scales
  scale_shape_manual(values = shape_values) +
  scale_color_manual(values = phylum_colors) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  theme_test() +
  facet_wrap(~ Group, ncol = 2, scales = "free") +
  labs(
    title = "Microbial Response Correlation: Genotype vs Treatment",
    x = expression(log[2]*FC~(italic(per5)-mock~vs~italic(Ws4)-mock)),
    y = expression(log[2]*FC~(italic(Ws4)-Treatment~vs~italic(Ws4)-mock)),
    shape = "Significance (p.adj < 0.05)",
    color = "Taxonomy"
  ) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

print(faceted_plot)

# Saving as a high-res PDF to ensure space for labels
ggsave("per5_LFC_LFC_colored_labels_try1.pdf", width = 13, height = 8)
