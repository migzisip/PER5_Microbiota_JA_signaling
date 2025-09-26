#20250909 updated by Chiao-Jung Han
##microbial profiling for per5 project
##bubble plot for spike-in normalized abudance
##rbec
setwd("C://Users//user//Desktop//per5 microbial profiling//Rbec//bubble plot")
library(phyloseq)
library(dplyr)
library(ggplot2)
library(biomformat)
library(ggtext)
library(forcats)

# Read the table and set OTU IDs as row names
otu_df <- read.table("rbec_absolute_no_input.txt", header = TRUE, sep = "\t", row.names = 1, check.names = TRUE)
otu_table_obj <-otu_table(otu_df,taxa_are_rows = TRUE)
metadata <- read.table("metadata_no_input.txt", header = TRUE, sep = "\t", row.names = 1)
sample_data_obj <- sample_data(metadata)  # converts to phyloseq sample_data object
taxonomy <- read.table("taxonomy.tsv", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
tax_table_obj <- tax_table(as.matrix(taxonomy))
physeq <- phyloseq(otu_table_obj, sample_data_obj, tax_table_obj)

# Make sure sample IDs match OTU table
sample_names(physeq)
rownames(metadata)

# Extract absolute abundances as matrix
# Calculate scaling factor per sample
spike <- read.csv("spike_no_input.csv",  row.names = 1, header= TRUE)
spike_counts <- spike 
spike_normalized_abundance <- otu_table_obj / spike_counts
abs_matrix<-spike_normalized_abundance ##spike normalized counts
write.csv(abs_matrix, file = "abs_matrix.csv", row.names = TRUE) ##check the spike-normalized read counts

# Create new phyloseq object with absolute abundances
abs_otu <- otu_table(abs_matrix, taxa_are_rows = TRUE)
abs_physeq <- phyloseq(abs_otu, sample_data(physeq), tax_table(physeq))

abs_genus <- tax_glom(abs_physeq, taxrank = "Genus")
plot_bar(abs_genus, fill = "Genus") +
  ylab("Spike-in normalized abundance") + theme(legend.position = "right")+
  theme_minimal()


library(RColorBrewer)

# Display available palettes (optional)
# RColorBrewer::display.brewer.all()

# Use Set3 (which supports up to 12 colors) + combine two palettes
my_colors <- c(brewer.pal(12, "Set3"), brewer.pal(6, "Paired"))  # Total = 18
plot_bar(abs_genus, fill = "Genus") +
  ylab("Spike-in normalized abundance") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = my_colors)


##draw the bubble plot
df <- psmelt(abs_genus)  # or abs_physeq if you're not aggregating to genus


##to remove the input/NA column
abs_family <- tax_glom(abs_physeq, taxrank = "Family")

#Create a named vector for mapping
unique(df$GroupLabel)
df$GroupLabel <- paste(df$compartment, df$genotype, df$treatment, sep = "\n")
df$GroupLabel <- factor(df$GroupLabel, levels = unique(df$GroupLabel[order(df$compartment, df$genotype, df$treatment)]))
df <- df %>% filter(GroupLabel != "input\nNA\nNA")
#make a customized order sequence
df$GroupLabel <-factor(df$GroupLabel, 
                       levels = c("root\nWs4\nNA","root\nWs4\nDCB","root\nWs4\nMeJA","root\nper5\nNA",
                                  "shoot\nWs4\nNA","shoot\nWs4\nDCB","shoot\nWs4\nMeJA", "shoot\nper5\nNA", 
                                  "agar\nWs4\nNA","agar\nWs4\nDCB","agar\nWs4\nMeJA","agar\nper5\nNA"
                                 ))

df <- df %>%
  filter(!is.na(Family) & Family != "" & Family != "unknown")

# Arrange the data frame by Phylum first, then class, and lastly by Family
df <- df %>%
  arrange(Phylum,Class, Family)

# Re-factor Family using that order
df$Family <- factor(df$Family, levels = unique(df$Family))
df$GroupLabel <- df %>% 
  arrange(Phylum, Family) %>%
  mutate(Family = factor(Family, levels = unique(Family)))

df$GroupLabel <- df$GroupLabel$GroupLabel
##re level the class sequence
df$Phylum <- factor(df$Phylum, levels = rev(levels(factor(df$Phylum))))



ggplot(df, aes(x = GroupLabel, y = Family, size = Abundance, fill = Phylum)) +
  geom_point(shape = 21, color = "black", alpha = 0.8) +
  scale_size_continuous(name = "Spike in-normalized abundance", range = c(1, 10)) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Spike-in normalized abundance (Rbec)",
       x = "condition",
       y = "Family") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##above mentioned is by sample, let's do a plot only present their avg
df_avg <- df %>%
  group_by(GroupLabel, Family, Class, Phylum) %>%
  arrange(Phylum, Class,Family) %>%
  summarise(AvgAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()

##re level the class sequence
df_avg$Class <- factor(df_avg$Class, levels = unique(df_avg$Class))
df_avg$Family <- factor(df_avg$Family, levels = rev(unique(df_avg$Family)))

##my color set
my_color<- c("Actinobacteria" = "#F6AEB0","Bacteroidetes"="#A4E0E2","Firmicutes"="#F69230","Alphaproteobacteria"="#AFDFC9","Betaproteobacteria"="#16946D","Gammaproteobacteria"="#708A63")
##make a hierachical color sets for text
family_class_df<- data.frame(taxonomy$Family, taxonomy$Class)
class_colors <-my_color
family_class_df$color<- class_colors[family_class_df$taxonomy.Class]
family_class_df<-family_class_df[-c(1,2), , drop = FALSE]

#plot
ggplot(df_avg, aes(x = GroupLabel, y = Family, size = AvgAbundance, fill = Class)) +
  geom_point(shape = 21, color = "#A9A9A9",alpha = 1.2, stroke = 0) +
  scale_colour_manual(values = my_color) +
  scale_fill_manual(values = my_color) +
  labs(title = " Spike-in normalized abundance (Rbec)",
       x = "Genotype + Compartment + Treatment",
       y = "Family") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.text = element_text(size = 9))+ 
  scale_size_continuous(name = "Normalized Abundance", trans = "sqrt", range = c(1, 10))

##make sure the group label id
unique(df_avg$GroupLabel)
##run statistics (Wilcoxon test)
library(dplyr)
library(rstatix)
library(tidyr)
library(tibble)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("group_by", "dplyr")
conflict_prefer("select", "dplyr")

#define the statistical matrix
my_data<-df

# Define groups
ref1 <- "agar\nWs4\nNA"
grp1 <- c("agar\nper5\nNA", "agar\nWs4\nDCB", "agar\nWs4\nMeJA")

ref2 <- "root\nWs4\nNA"
grp2 <- c("root\nper5\nNA", "root\nWs4\nDCB", "root\nWs4\nMeJA")

ref3 <- "shoot\nWs4\nNA"
grp3 <- c("shoot\nper5\nNA", "shoot\nWs4\nDCB", "shoot\nWs4\nMeJA")


# Helper function to run wilcoxon test between each group and its reference
run_group_test <- function(df_input, test_groups, ref_group, set_id) {
  all_results <- list()
  families <- unique(df_input$Family)
  
  for (fam in families) {
    fam_data <- df_input %>% filter(Family == fam)
    
    for (grp in test_groups) {
      x_vals <- fam_data %>% filter(GroupLabel == grp) %>% pull(Abundance)
      y_vals <- fam_data %>% filter(GroupLabel == ref_group) %>% pull(Abundance)
      
      if (length(x_vals) > 0 && length(y_vals) > 0) {
        # Turn OFF confidence interval to avoid edge-case crash
        wt <- wilcox.test(x_vals, y_vals, conf.int = FALSE)
        
        result_row <- tibble::tibble(
          Family = fam,
          group1 = grp,
          group2 = ref_group,
          p = wt$p.value,
          statistic = wt$statistic,
          estimate = median(x_vals) - median(y_vals),
          ComparisonSet = set_id
        )
        
        all_results[[length(all_results) + 1]] <- result_row
      }
    }
  }
  
  bind_rows(all_results)
}

# Run tests for all 3 groups
stats1 <- run_group_test(my_data, grp1, ref1, "Set1")
stats2 <- run_group_test(my_data, grp2, ref2, "Set2")
stats3 <- run_group_test(my_data, grp3, ref3, "Set3")
warnings()


# Combine and adjust p-values within each set
pairwise_stats <- bind_rows(stats1, stats2, stats3) %>%
  group_by(Family, ComparisonSet) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  ungroup()

# Assign significance categories
sig_labels <- pairwise_stats %>%
  transmute(
    label = paste(Family, group1, sep = "_"),
    Significance = case_when(
      p.adj <= 0.05 ~ "p ≤ 0.05",
      p.adj <= 0.1 ~ "0.05 < p ≤ 0.1",
      TRUE ~ "ns"
    )
  ) %>%
  distinct()

# Add label to df_avg for joining
df_avg <- df_avg %>%
  mutate(label = paste(Family, GroupLabel, sep = "_"))

# Join significance labels
df_avg <- df_avg %>%
  left_join(sig_labels, by = "label")



# Merge and annotate df_avg
df_avg <-df_avg %>%
  mutate(
    Significance = coalesce(Significance, "ns"),
    stroke_width = case_when(
      Significance == "p ≤ 0.05" ~ 3.5,
      Significance == "0.05 < p ≤ 0.1" ~ 1,
      TRUE ~ 0
    ))

df_avg <- df_avg %>%
  mutate(
    Significance = coalesce(Significance, "ns"),
    stroke_color = case_when(
      Significance == "p ≤ 0.05" ~ "black",
      Significance == "0.05 < p ≤ 0.1" ~ "#A9A9A9",
      TRUE ~ "white"
    ))
##re arrange the level
df_avg$Significance <-factor(df_avg$Significance, levels = c("p ≤ 0.05","0.05 < p ≤ 0.1", "ns")) 

# Annotate df_avg with significance
df_avg <- df_avg %>%
  mutate(
    Significance = coalesce(Significance, "ns"),
    Significance = factor(Significance, levels = c("p ≤ 0.05", "0.05 < p ≤ 0.1", "ns"))
  )

##plot with sig-label
ggplot(df_avg, aes(x = GroupLabel, y = Family,
                   size = AvgAbundance, fill = Class)) +
  geom_point(aes(
    shape = Significance,
    color = Significance,
    stroke = 1.5
  ), alpha = 0.9) +
  scale_size_continuous(
    name = "Normalized Abundance",
    trans = "sqrt", range = c(1, 10),
    breaks = c(0.01,0.1, 0.5,1,2.2)
  ) +
  scale_fill_manual(values = my_color) +
  scale_shape_manual(
    values = c("p ≤ 0.05" = 21, "0.05 < p ≤ 0.1" = 21, "ns" = 21),
    guide = guide_legend()
  ) +
  scale_color_manual(
    values = c("p ≤ 0.05" = "black", "0.05 < p ≤ 0.1" = "#A9A9A9", "ns" = "white"),
    guide = guide_legend()
  ) +
  labs(
    title = "Spike-in Normalized Abundance (Rbec)",
    x = "Genotype + Compartment + Treatment",
    y = "Family"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

