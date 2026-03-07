library(rstudioapi)
# Get the path of the current script
path <- rstudioapi::getActiveDocumentContext()$path
# Set the working directory to the folder containing that script
setwd(dirname(path))
# Verify the change

library(vegan)
library(phyloseq)
#Load data + build phyloseq
otu_df <- read.table("asv_table_ws4per5_not_normalized.txt",
                     header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
metadata <- read.table("ws4per5_design_BC.txt",
                       header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
taxonomy <- read.table("taxonomy_all.txt",
                       header = TRUE, sep = "\t", row.names = 1,
                       stringsAsFactors = FALSE, check.names = FALSE)
physeq <- phyloseq(
  otu_table(as.matrix(otu_df), taxa_are_rows = TRUE),
  sample_data(metadata),
  tax_table(as.matrix(taxonomy))
)

###anosim
ps_rel <- transform_sample_counts(physeq, function(x) x/sum(x) )

# --- ROOT SUBSET ---
ps_root <- subset_samples(ps_rel, compartment == "root")
# Remove taxa that are now zero in this specific subset
ps_root <- prune_taxa(taxa_sums(ps_root) > 0, ps_root)

# Distance matrix for Root only
d_root <- phyloseq::distance(ps_root, method = "bray")
meta_root <- as(sample_data(ps_root), "data.frame")

# ANOSIMs for Root
ano_root_genotype  <- vegan::anosim(d_root, grouping = meta_root$genotype, permutations = 999)
ano_root_genotype

# --- SHOOT SUBSET ---
ps_shoot <- subset_samples(ps_rel, compartment == "shoot")
ps_shoot <- prune_taxa(taxa_sums(ps_shoot) > 0, ps_shoot)

# Distance matrix for Shoot only
d_shoot <- phyloseq::distance(ps_shoot, method = "bray")
meta_shoot <- as(sample_data(ps_shoot), "data.frame")

# ANOSIMs for Shoot
ano_shoot_genotype  <- vegan::anosim(d_shoot, grouping = meta_shoot$genotype, permutations = 999)
ano_shoot_genotype

# --- MATRIX SUBSET ---
ps_matrix <- subset_samples(ps_rel, compartment == "matrix")
ps_matrix <- prune_taxa(taxa_sums(ps_matrix) > 0, ps_matrix)

# Distance matrix for matrix only
d_matrix <- phyloseq::distance(ps_matrix, method = "bray")
meta_matrix <- as(sample_data(ps_matrix), "data.frame")

# ANOSIMs for matrix
ano_matrix_genotype  <- vegan::anosim(d_matrix, grouping = meta_matrix$genotype, permutations = 999)
ano_matrix_genotype

### PERMANOVA (adonis2) + betadisper for spike-in normalized abundance (RBEC)
# 2) Transform + distance
#    (composition-only) Bray–Curtis on relative abundance
#PERMANOVA (adonis2)
#full factorial model ##root
ado_root <- vegan::adonis2(d_root ~  genotype, data = meta_root,permutations = 999, strata=meta_root$biological.replicate)
ado_root 
print(ado_root)
write.table(as.data.frame(ado_root), file = "adonis2_root_genotype_RA.txt", sep = "\t", quote = FALSE,row.names = TRUE,col.names = NA)

#full factorial model ##shoot
ado_shoot <- vegan::adonis2(d_shoot ~ genotype ,data = meta_shoot,permutations = 999, strata=meta_shoot$biological.replicate)
ado_shoot 
print(ado_shoot)
write.table(as.data.frame(ado_shoot), file = "adonis2_shoot_genotype_RA.txt", sep = "\t", quote = FALSE,row.names = TRUE,col.names = NA)

#full factorial model ##matrix
ado_matrix <- vegan::adonis2(d_matrix ~ genotype ,data = meta_matrix,permutations = 999, strata=meta_matrix$biological.replicate)
ado_matrix
print(ado_matrix)
write.table(as.data.frame(ado_matrix), file = "adonis2_matrix_genotype_RA.txt", sep = "\t", quote = FALSE,row.names = TRUE,col.names = NA)

#  betadisper (dispersion test)
#  Run for the SAME grouping used in PERMANOVA
# dispersion for genotype_compartment_treatment
# --- ROOT BETADISPER ---
# 1. Betadisper for Genotype_Treatment (most specific group)
disper_root <- vegan::betadisper(d_root, group = meta_root$genotype)
# Test if the dispersions (variances) are significantly different
disper_root_anova <- anova(disper_root)
disper_root_test<-permutest(disper_root, permutations = 999)
print(disper_root_test)

# --- SHOOT BETADISPER ---
disper_shoot <- vegan::betadisper(d_shoot, group = meta_shoot$genotype)
disper_shoot_anova <- anova(disper_shoot)
disper_shoot_test<-permutest(disper_shoot, permutations = 999)
print(disper_shoot_test)

# --- Matrix BETADISPER ---
disper_matrix <- vegan::betadisper(d_matrix, group = meta_matrix$genotype)
disper_matrix_anova <- anova(disper_matrix)
disper_matrix_test<-permutest(disper_matrix, permutations = 999)
print(disper_matrix_test)

##summary of betadisper
disper_root_test #genotype_treatment
disper_shoot_test #genotype_treatment
disper_matrix_test
