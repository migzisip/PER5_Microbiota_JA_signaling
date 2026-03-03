library(rstudioapi)
# Get the path of the current script
path <- rstudioapi::getActiveDocumentContext()$path
# Set the working directory to the folder containing that script
setwd(dirname(path))
# Verify the change
getwd()

library(vegan)
library(phyloseq)
library(pairwiseAdonis)
#Load data + build phyloseq
otu_df <- read.table("rbec_absolute_no_input.txt",
                     header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
metadata <- read.table("metadata_no_input.txt",
                       header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
taxonomy <- read.table("taxonomy.txt",
                       header = TRUE, sep = "\t", row.names = 1,
                       stringsAsFactors = FALSE, check.names = FALSE)
physeq <- phyloseq(
  otu_table(as.matrix(otu_df), taxa_are_rows = TRUE),
  sample_data(metadata),
  tax_table(as.matrix(taxonomy))
)

###anosim
ps_rel <- transform_sample_counts(physeq, function(x) x)
# --- ROOT SUBSET ---
ps_root <- subset_samples(ps_rel, compartment == "root")
# Remove taxa that are now zero in this specific subset
ps_root <- prune_taxa(taxa_sums(ps_root) > 0, ps_root)

# Distance matrix for Root only
d_root <- phyloseq::distance(ps_root, method = "bray")
meta_root <- as(sample_data(ps_root), "data.frame")

# ANOSIMs for Root
ano_root_inter<- vegan::anosim(d_root, grouping = meta_root$genotype_treatment, permutations = 999)
ano_root_inter 

##pairwise anosim-ws4_DCB/MeJA vs per5-mock 
root_group_to_compare_DCB <- c("Ws4_DCB", "per5_NA")
root_meta_subset_DCB <- meta_root[meta_root$genotype_treatment %in% root_group_to_compare_DCB, ]
root_d_subset_DCB <- as.dist(as.matrix(d_root)[rownames(root_meta_subset_DCB), rownames(root_meta_subset_DCB)])
root_ano_pairwise_DCB <- vegan::anosim(  root_d_subset_DCB, grouping = root_meta_subset_DCB$genotype_treatment,permutations = 999)
print(root_ano_pairwise_DCB)
summary(root_ano_pairwise_DCB)
ano_results_to_save <- data.frame(Comparison = paste(root_group_to_compare_DCB, collapse = "_vs_"),R_statistic = root_ano_pairwise_DCB$statistic,P_value = root_ano_pairwise_DCB$signif)
write.table(ano_results_to_save, "anosim_pairwise_Ws4-DCB vs per5-mock_root.txt", sep = "\t", row.names = FALSE)

root_group_to_compare_MeJA <- c("Ws4_MeJA", "per5_NA")
root_meta_subset_MeJA <- meta_root[meta_root$genotype_treatment %in% root_group_to_compare_MeJA, ]
root_d_subset_MeJA <- as.dist(as.matrix(d_root)[rownames(root_meta_subset_MeJA), rownames(root_meta_subset_MeJA)])
root_ano_pairwise_MeJA <- vegan::anosim(root_d_subset_MeJA, grouping = root_meta_subset_MeJA$genotype_treatment,permutations = 999)
print(root_ano_pairwise_MeJA)
summary(root_ano_pairwise_MeJA)
ano_results_to_save <- data.frame(Comparison = paste(root_group_to_compare_MeJA, collapse = "_vs_"),R_statistic = root_ano_pairwise_MeJA$statistic,P_value = root_ano_pairwise_MeJA$signif)
write.table(ano_results_to_save, "anosim_pairwise_Ws4-MeJA vs per5-mock_root.txt", sep = "\t", row.names = FALSE)

#pairwise anosim-ws4-mock vs ws4 DCB/MeJA
root_group_to_compare_NA_MeJA <- c("Ws4_NA","Ws4_MeJA")
root_meta_subset_NA_MeJA <- meta_root[meta_root$genotype_treatment %in% root_group_to_compare_NA_MeJA, ]
root_d_subset_NA_MeJA <- as.dist(as.matrix(d_root)[rownames(root_meta_subset_NA_MeJA), rownames(root_meta_subset_NA_MeJA)])
root_ano_pairwise_NA_MeJA <- vegan::anosim(root_d_subset_NA_MeJA, grouping = root_meta_subset_NA_MeJA$genotype_treatment,permutations = 999)
print(root_ano_pairwise_NA_MeJA)
summary(root_ano_pairwise_NA_MeJA)
ano_results_to_save <- data.frame(Comparison = paste(root_group_to_compare_NA_MeJA, collapse = "_vs_"),R_statistic = root_ano_pairwise_NA_MeJA$statistic,P_value = root_ano_pairwise_NA_MeJA$signif)
write.table(ano_results_to_save, "anosim_pairwise_Ws4-MeJA vs Ws4-mock_root.txt", sep = "\t", row.names = FALSE)

root_group_to_compare_NA_DCB <- c("Ws4_NA","Ws4_DCB")
root_meta_subset_NA_DCB <- meta_root[meta_root$genotype_treatment %in% root_group_to_compare_NA_DCB, ]
root_d_subset_NA_DCB <- as.dist(as.matrix(d_root)[rownames(root_meta_subset_NA_DCB), rownames(root_meta_subset_NA_DCB)])
root_ano_pairwise_NA_DCB <- vegan::anosim(root_d_subset_NA_DCB, grouping = root_meta_subset_NA_DCB$genotype_treatment,permutations = 999)
print(root_ano_pairwise_NA_DCB)
summary(root_ano_pairwise_NA_DCB)
ano_results_to_save <- data.frame(Comparison = paste(root_group_to_compare_NA_DCB, collapse = "_vs_"),R_statistic = root_ano_pairwise_NA_DCB$statistic,P_value = root_ano_pairwise_NA_DCB$signif)
write.table(ano_results_to_save, "anosim_pairwise_Ws4-DCB vs Ws4-mock_root.txt", sep = "\t", row.names = FALSE)

#pairwise anosim-ws4-mock vs per5-mock
root_group_to_compare_NA <- c("Ws4_NA","per5_NA")
root_meta_subset_NA <- meta_root[meta_root$genotype_treatment %in% root_group_to_compare_NA, ]
root_d_subset_NA <- as.dist(as.matrix(d_root)[rownames(root_meta_subset_NA), rownames(root_meta_subset_NA)])
root_ano_pairwise_NA <- vegan::anosim(root_d_subset_NA, grouping = root_meta_subset_NA$genotype_treatment,permutations = 999)
print(root_ano_pairwise_NA)
summary(root_ano_pairwise_NA)
ano_results_to_save <- data.frame(Comparison = paste(root_group_to_compare_NA, collapse = "_vs_"),R_statistic = root_ano_pairwise_NA$statistic,P_value = root_ano_pairwise_NA$signif)
write.table(ano_results_to_save, "anosim_pairwise_Ws4-mock vs per5-mock_root.txt", sep = "\t", row.names = FALSE)


# --- SHOOT SUBSET ---
ps_shoot <- subset_samples(ps_rel, compartment == "shoot")
ps_shoot <- prune_taxa(taxa_sums(ps_shoot) > 0, ps_shoot)

# Distance matrix for Shoot only
d_shoot <- phyloseq::distance(ps_shoot, method = "bray")
meta_shoot <- as(sample_data(ps_shoot), "data.frame")

# ANOSIMs for Shoot
ano_shoot_inter     <- vegan::anosim(d_shoot, grouping = meta_shoot$genotype_treatment, permutations = 999)
ano_shoot_inter  

##pairwise anosim-ws4_DCB/MeJA vs per5-mock 
shoot_group_to_compare_DCB <- c("Ws4_DCB", "per5_NA")
shoot_meta_subset_DCB <- meta_shoot[meta_shoot$genotype_treatment %in% shoot_group_to_compare_DCB, ]
shoot_d_subset_DCB <- as.dist(as.matrix(d_shoot)[rownames(shoot_meta_subset_DCB), rownames(shoot_meta_subset_DCB)])
shoot_ano_pairwise_DCB <- vegan::anosim(  shoot_d_subset_DCB, grouping = shoot_meta_subset_DCB$genotype_treatment,permutations = 999)
print(shoot_ano_pairwise_DCB)
summary(shoot_ano_pairwise_DCB)
ano_results_to_save <- data.frame(Comparison = paste(shoot_group_to_compare_DCB, collapse = "_vs_"),R_statistic = shoot_ano_pairwise_DCB$statistic,P_value = shoot_ano_pairwise_DCB$signif)
write.table(ano_results_to_save, "anosim_pairwise_Ws4-DCB vs per5-mock_shoot.txt", sep = "\t", row.names = FALSE)

shoot_group_to_compare_MeJA <- c("Ws4_MeJA", "per5_NA")
shoot_meta_subset_MeJA <- meta_shoot[meta_shoot$genotype_treatment %in% shoot_group_to_compare_MeJA, ]
shoot_d_subset_MeJA <- as.dist(as.matrix(d_shoot)[rownames(shoot_meta_subset_MeJA), rownames(shoot_meta_subset_MeJA)])
shoot_ano_pairwise_MeJA <- vegan::anosim(shoot_d_subset_MeJA, grouping = shoot_meta_subset_MeJA$genotype_treatment,permutations = 999)
print(shoot_ano_pairwise_MeJA)
summary(shoot_ano_pairwise_MeJA)
ano_results_to_save <- data.frame(Comparison = paste(shoot_group_to_compare_MeJA, collapse = "_vs_"),R_statistic = shoot_ano_pairwise_MeJA$statistic,P_value = shoot_ano_pairwise_MeJA$signif)
write.table(ano_results_to_save, "anosim_pairwise_Ws4-MeJA vs per5-mock_shoot.txt", sep = "\t", row.names = FALSE)

#pairwise anosim-ws4-mock vs ws4 DCB/MeJA
shoot_group_to_compare_NA_MeJA <- c("Ws4_NA","Ws4_MeJA")
shoot_meta_subset_NA_MeJA <- meta_shoot[meta_shoot$genotype_treatment %in% shoot_group_to_compare_NA_MeJA, ]
shoot_d_subset_NA_MeJA <- as.dist(as.matrix(d_shoot)[rownames(shoot_meta_subset_NA_MeJA), rownames(shoot_meta_subset_NA_MeJA)])
shoot_ano_pairwise_NA_MeJA <- vegan::anosim(shoot_d_subset_NA_MeJA, grouping = shoot_meta_subset_NA_MeJA$genotype_treatment,permutations = 999)
print(shoot_ano_pairwise_NA_MeJA)
summary(shoot_ano_pairwise_NA_MeJA)
ano_results_to_save <- data.frame(Comparison = paste(shoot_group_to_compare_NA_MeJA, collapse = "_vs_"),R_statistic = shoot_ano_pairwise_NA_MeJA$statistic,P_value = shoot_ano_pairwise_NA_MeJA$signif)
write.table(ano_results_to_save, "anosim_pairwise_Ws4-MeJA vs Ws4-mock_shoot.txt", sep = "\t", row.names = FALSE)

shoot_group_to_compare_NA_DCB <- c("Ws4_NA","Ws4_DCB")
shoot_meta_subset_NA_DCB <- meta_shoot[meta_shoot$genotype_treatment %in% shoot_group_to_compare_NA_DCB, ]
shoot_d_subset_NA_DCB <- as.dist(as.matrix(d_shoot)[rownames(shoot_meta_subset_NA_DCB), rownames(shoot_meta_subset_NA_DCB)])
shoot_ano_pairwise_NA_DCB <- vegan::anosim(shoot_d_subset_NA_DCB, grouping = shoot_meta_subset_NA_DCB$genotype_treatment,permutations = 999)
print(shoot_ano_pairwise_NA_DCB)
summary(shoot_ano_pairwise_NA_DCB)
ano_results_to_save <- data.frame(Comparison = paste(shoot_group_to_compare_NA_DCB, collapse = "_vs_"),R_statistic = shoot_ano_pairwise_NA_DCB$statistic,P_value = shoot_ano_pairwise_NA_DCB$signif)
write.table(ano_results_to_save, "anosim_pairwise_Ws4-DCB vs Ws4-mock_shoot.txt", sep = "\t", row.names = FALSE)

#pairwise anosim-ws4-mock vs per5-mock
shoot_group_to_compare_NA <- c("Ws4_NA","per5_NA")
shoot_meta_subset_NA <- meta_shoot[meta_shoot$genotype_treatment %in% shoot_group_to_compare_NA, ]
shoot_d_subset_NA <- as.dist(as.matrix(d_shoot)[rownames(shoot_meta_subset_NA), rownames(shoot_meta_subset_NA)])
shoot_ano_pairwise_NA <- vegan::anosim(shoot_d_subset_NA, grouping = shoot_meta_subset_NA$genotype_treatment,permutations = 999)
print(shoot_ano_pairwise_NA)
summary(shoot_ano_pairwise_NA)
ano_results_to_save <- data.frame(Comparison = paste(shoot_group_to_compare_NA, collapse = "_vs_"),R_statistic = shoot_ano_pairwise_NA$statistic,P_value = shoot_ano_pairwise_NA$signif)
write.table(ano_results_to_save, "anosim_pairwise_Ws4-mock vs per5-mock_shoot.txt", sep = "\t", row.names = FALSE)


### PERMANOVA (adonis2) + betadisper for spike-in normalized abundance (RBEC)
# 2) Transform + distance
#    (composition-only) Bray–Curtis on relative abundance
#PERMANOVA (adonis2)
# global test (n=1)
ado_root_genotype_treatment <- vegan::adonis2(d_root ~ genotype_treatment,data = meta_root,permutations = 999, strata= meta_root$genotype_treatment)
ado_root_genotype_treatment
write.table(as.data.frame(ado_root_genotype_treatment), file = "adonis2_root_genotype_treatment.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

##BH correction for p
# Apply Benjamini-Hochberg correction
corrected_p <- p.adjust(ado_root_genotype_treatment$`Pr(>F)`, method = "BH")
data.frame(ado_root_genotype_treatment, corrected_p) ##not neccessary for the global one

###pairwise adonis2
pairwise_results_root <- pairwise.adonis(
  x = d_root, 
  factors = meta_root$genotype_treatment, 
  sim.method = "bray",         
  p.adjust.m = "BH",           
  perm = 999                   
  )
pairwise_results_root
write.table(as.data.frame(pairwise_results_root), file = "pairwise_results_root.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

####manually run BH correction to double check
# 1. Define your data and grouping factor
dat <- d_root
groups <- as.factor(meta_root$genotype_treatment)
dat
# 2. Get all 6 combinations for the 4 groups
group_levels <- levels(groups)
pairs <- combn(group_levels, 2)

# 3. Initialize a table to store results
# We'll store the Pair name, the R2 (Effect Size), and the raw P-value
pairwise_df <- data.frame(Pair = character(),
                          R2 = numeric(),
                          P_value = numeric(),
                          stringsAsFactors = FALSE)

# 4. Loop through each of the 6 pairs
for(i in 1:ncol(pairs)) {
  p1 <- pairs[1, i]
  p2 <- pairs[2, i]
  
  # Create the index
  idx <- groups %in% c(p1, p2)
  
  # --- THE FIX STARTS HERE ---
  # If d_root is a distance matrix, we convert to matrix to subset rows AND columns
  if(inherits(dat, "dist")) {
    sub_dist <- as.matrix(dat)[idx, idx]
    sub_dist <- as.dist(sub_dist)
  } else {
    # If d_root is a community table, we subset rows normally
    sub_dist <- vegdist(dat[idx, ], method = "bray")
  }
  sub_groups <- groups[idx]
  # --- THE FIX ENDS HERE ---
  
  # Run adonis2 on the subsetted distance
  ad_res <- adonis2(sub_dist ~ sub_groups)
  
  pairwise_df[i, "Pair"] <- paste(p1, "vs", p2)
  pairwise_df[i, "R2"] <- ad_res$R2[1]
  pairwise_df[i, "P_value"] <- ad_res$`Pr(>F)`[1]
}

# 3. Apply Correction
pairwise_df$P_adj <- p.adjust(pairwise_df$P_value, method = "BH")
print(pairwise_df)
###end of double check

#full factorial model ##shoot
#global test (n=1)
ado_shoot_genotype_treatment <- vegan::adonis2(d_shoot ~ genotype_treatment,data = meta_shoot,permutations = 999, by = "terms")
ado_shoot_genotype_treatment
write.table(as.data.frame(ado_shoot_genotype_treatment), file = "adonis2_shoot_genotype_treatment.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
##BH correction for p
# Apply Benjamini-Hochberg correction
corrected_p <- p.adjust(ado_shoot_genotype_treatment$`Pr(>F)`, method = "BH")
data.frame(ado_shoot_genotype_treatment, corrected_p) ##not neccessary for the global one

###pairwise adonis2
pairwise_results_shoot <- pairwise.adonis(
  x = d_shoot, 
  factors = meta_shoot$genotype_treatment, 
  sim.method = "bray",         
  p.adjust.m = "BH",           
  perm = 999                   
)
pairwise_results_shoot
write.table(as.data.frame(pairwise_results_shoot), file = "pairwise_results_shoot.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

####manually run BH correction to double check
# 1. Setup
dat2 <- d_shoot
groups <- as.factor(meta_shoot$genotype_treatment)

# 2. Pairs
group_levels <- levels(groups)
pairs <- combn(group_levels, 2)

# 3. Results table
pairwise_df2 <- data.frame(Pair = character(),
                           R2 = numeric(),
                           P_value = numeric(),
                           stringsAsFactors = FALSE)

# 4. Loop
for(i in 1:ncol(pairs)) {
  p1 <- pairs[1, i]
  p2 <- pairs[2, i]
  
  idx <- groups %in% c(p1, p2)
  
  # --- UPDATED TO USE dat2 ---
  if(inherits(dat2, "dist")) {
    # If it's already a distance matrix
    sub_dist <- as.matrix(dat2)[idx, idx]
    sub_dist <- as.dist(sub_dist)
  } else {
    # If it's a community table, calculate distance now
    sub_dist <- vegdist(dat2[idx, ], method = "bray")
  }
  
  sub_groups <- groups[idx]
  
  # Run adonis2
  ad_res <- adonis2(sub_dist ~ sub_groups)
  
  pairwise_df2[i, "Pair"] <- paste(p1, "vs", p2)
  pairwise_df2[i, "R2"] <- ad_res$R2[1]
  pairwise_df2[i, "P_value"] <- ad_res$`Pr(>F)`[1]
}

# 5. Apply BH Correction
pairwise_df2$P_adj <- p.adjust(pairwise_df2$P_value, method = "BH")

print(pairwise_df2)
###end of double check

#  betadisper (dispersion test)
#  Run for the SAME grouping used in PERMANOVA
# dispersion for genotype_compartment_treatment
# --- ROOT BETADISPER ---
# 1. Betadisper for Genotype_Treatment (most specific group)
disper_root <- vegan::betadisper(d_root, group = meta_root$genotype_treatment)
# Test if the dispersions (variances) are significantly different
disper_root_anova <- anova(disper_root)
disper_root_test<-permutest(disper_root, permutations = 999)
print(disper_root_test)

# --- ROOT BETADISPER ---
disper_shoot <- vegan::betadisper(d_shoot, group = meta_shoot$genotype_treatment)
disper_shoot_anova <- anova(disper_shoot)
disper_shoot_test<-permutest(disper_shoot, permutations = 999)
print(disper_shoot_test)
# --- PAIRWISE BETADISPER for Ws4_DCB vs Ws4_NA (root)
disper_pairwise_DCB_NA_root <- vegan::betadisper(root_d_subset_NA_DCB, group = root_meta_subset_NA_DCB$genotype_treatment)
anova(disper_pairwise_DCB_NA_root)
disper_pairwise_DCB_NA_root_test<-permutest(disper_pairwise_DCB_NA_root, permutations = 999)
disper_pairwise_DCB_NA_root_test
tukey_results_ws4_DCB_NA_root <- TukeyHSD(disper_pairwise_DCB_NA_root)
print(tukey_results_ws4_DCB_NA_root)

# --- PAIRWISE BETADISPER for Ws4_MeJA vs Ws4_NA (root)
disper_pairwise_MeJA_NA_root <- vegan::betadisper(root_d_subset_NA_MeJA, group = root_meta_subset_NA_MeJA$genotype_treatment)
anova(disper_pairwise_MeJA_NA_root)
disper_pairwise_MeJA_NA_root_test<-permutest(disper_pairwise_MeJA_NA_root, permutations = 999)
disper_pairwise_MeJA_NA_root_test
tukey_results_ws4_MeJA_per5_root <- TukeyHSD(disper_pairwise_MeJA_NA_root)
print(tukey_results_ws4_MeJA_per5_root)

# --- PAIRWISE BETADISPER for Ws4_NA vs Ws4_NA (root)
disper_pairwise_NA_root <- vegan::betadisper(root_d_subset_NA, group = root_meta_subset_NA$genotype_treatment)
anova(disper_pairwise_NA_root)
disper_pairwise_root_test<-permutest(disper_pairwise_NA_root, permutations = 999)
disper_pairwise_root_test
tukey_results_ws4_per5_root <- TukeyHSD(disper_pairwise_NA_root)
print(tukey_results_ws4_per5_root)
# --- PAIRWISE BETADISPER for Ws4_DCB vs per5_NA (root)
disper_pairwise_DCB_root <- vegan::betadisper(root_d_subset_DCB, group = root_meta_subset_DCB$genotype_treatment)
anova(disper_pairwise_DCB_root)
disper_pairwise_DCB_root_test<-permutest(disper_pairwise_DCB_root, permutations = 999)
disper_pairwise_DCB_root_test
tukey_results_ws4_DCB_per5 <- TukeyHSD(disper_pairwise_DCB_root)
print(tukey_results_ws4_DCB_per5)
# --- PAIRWISE BETADISPER for Ws4_MeJA vs per5_NA (root)
disper_pairwise_MeJA_root <- vegan::betadisper(root_d_subset_MeJA, group = root_meta_subset_MeJA$genotype_treatment)
anova(disper_pairwise_MeJA_root)
disper_pairwise_MeJA_root_test<-permutest(disper_pairwise_MeJA_root, permutations = 999)
disper_pairwise_MeJA_root_test
tukey_results_ws4_MeJA_per5 <- TukeyHSD(disper_pairwise_MeJA_root)
print(tukey_results_ws4_MeJA_per5)

# Plot the groups and their distance to centroids ##just to check
##plot(disper_pairwise_MeJA_root, main = "Multivariate Dispersion", hull = FALSE, ellipse = TRUE)
# --- PAIRWISE BETADISPER for Ws4_DCB vs Ws4_NA (shoot)
disper_pairwise_DCB_NA_shoot <- vegan::betadisper(shoot_d_subset_NA_DCB, group = shoot_meta_subset_NA_DCB$genotype_treatment)
anova(disper_pairwise_DCB_NA_shoot)
disper_pairwise_DCB_NA_shoot_test<-permutest(disper_pairwise_DCB_NA_shoot, permutations = 999)
disper_pairwise_DCB_NA_shoot_test
tukey_results_ws4_DCB_NA_shoot <- TukeyHSD(disper_pairwise_DCB_NA_shoot)
print(tukey_results_ws4_DCB_NA_shoot)

# --- PAIRWISE BETADISPER for Ws4_MeJA vs Ws4_NA (shoot)
disper_pairwise_MeJA_NA_shoot <- vegan::betadisper(shoot_d_subset_NA_MeJA, group = shoot_meta_subset_NA_MeJA$genotype_treatment)
anova(disper_pairwise_MeJA_NA_shoot)
disper_pairwise_MeJA_NA_shoot_test<-permutest(disper_pairwise_MeJA_NA_shoot, permutations = 999)
disper_pairwise_MeJA_NA_shoot_test
tukey_results_ws4_MeJA_per5_shoot <- TukeyHSD(disper_pairwise_MeJA_NA_shoot)
print(tukey_results_ws4_MeJA_per5_shoot)

# --- PAIRWISE BETADISPER for Ws4_NA vs Ws4_NA (shoot)
disper_pairwise_NA_shoot <- vegan::betadisper(shoot_d_subset_NA, group = shoot_meta_subset_NA$genotype_treatment)
anova(disper_pairwise_NA_shoot)
disper_pairwise_shoot_test<-permutest(disper_pairwise_NA_shoot, permutations = 999)
disper_pairwise_shoot_test
tukey_results_ws4_per5_shoot <- TukeyHSD(disper_pairwise_NA_shoot)
print(tukey_results_ws4_per5_shoot)

# --- PAIRWISE BETADISPER for Ws4_DCB vs per5_NA (shoot)
disper_pairwise_DCB_shoot <- vegan::betadisper(shoot_d_subset_DCB, group = shoot_meta_subset_DCB$genotype_treatment)
anova(disper_pairwise_DCB_shoot)
disper_pairwise_DCB_shoot_test<-permutest(disper_pairwise_DCB_shoot, permutations = 999)
disper_pairwise_DCB_shoot_test
tukey_results_ws4_DCB_per5_shoot <- TukeyHSD(disper_pairwise_DCB_shoot)
print(tukey_results_ws4_DCB_per5_shoot)

# --- PAIRWISE BETADISPER for Ws4_MeJA vs per5_NA (shoot)
disper_pairwise_MeJA_shoot <- vegan::betadisper(shoot_d_subset_MeJA, group = shoot_meta_subset_MeJA$genotype_treatment)
anova(disper_pairwise_MeJA_shoot)
disper_pairwise_MeJA_shoot_test<-permutest(disper_pairwise_MeJA_shoot, permutations = 999)
disper_pairwise_MeJA_shoot_test
tukey_results_ws4_MeJA_per5_shoot <- TukeyHSD(disper_pairwise_MeJA_shoot)
print(tukey_results_ws4_MeJA_per5_shoot)

shoot_meta_subset_NA <- meta_shoot[meta_shoot$genotype_treatment %in% shoot_group_to_compare_NA, ]
shoot_d_subset_NA <- as.dist(as.matrix(d_shoot)[rownames(shoot_meta_subset_NA), rownames(shoot_meta_subset_NA)])

