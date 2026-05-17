#-------------------------------------------------------------------------------
# Title: M. tuberculosis Variant Analysis (combined VCF)
# Data: all_variants.csv from bcftools joint calling (haploid, filtered)
# Author: (your name)
#-------------------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(RColorBrewer)
library(flextable)

#-------------------------------------------------------------------------------
# 1. Load the combined variant data
#-------------------------------------------------------------------------------
variants <- read.csv("all_variants.csv", stringsAsFactors = FALSE)

# Quick look
dim(variants)          # number of variants × number of columns
head(variants, 3)
glimpse(variants)

#-------------------------------------------------------------------------------
# 2. Classify variants: SNP vs INDEL
#-------------------------------------------------------------------------------
variants <- variants %>%
  mutate(
    # SNP: both REF and ALT are single bases
    variant_type = ifelse(nchar(REF) == 1 & nchar(ALT) == 1, "SNP", "INDEL")
  )
view(variants)
# Count and proportion
counts_summary <- table(variants$variant_type)
counts_prop <- prop.table(table(variants$variant_type))

# exporting output to tables in word
variant_counts <- as.data.frame(counts_summary)
variant_counts %>%
  flextable() %>%
  autofit() %>%
  set_caption("Summary of variants") %>% 
  save_as_docx(path = "variant_count_summary.docx")

variant_prop <- as.data.frame(counts_prop)
variant_prop %>%
  flextable() %>%
  autofit() %>%
  set_caption("Variants proportions") %>% 
  save_as_docx(path = "variant_proportion.docx")

#------------------counts#-------------------------------------------------------------------------------
# 3. Summary statistics of QUAL and DP
#-------------------------------------------------------------------------------
summary(variants$QUAL)
summary(variants$DP)

# Quality distribution by variant type
qual_dist <- variants %>%
  mutate(
    QUAL = as.numeric(QUAL),
    DP = as.numeric(DP)
  ) %>%
  group_by(variant_type) %>%
  summarise(
    mean_QUAL   = mean(QUAL, na.rm = TRUE),
    median_QUAL = median(QUAL, na.rm = TRUE),
    mean_DP     = mean(DP, na.rm = TRUE),
    median_DP   = median(DP, na.rm = TRUE)
  )
qual_dist
#-------------------------------------------------------------------------------
# 4. Reshape genotype columns into long format (one row per variant per sample)
#-------------------------------------------------------------------------------
# Genotype columns start with "GT_"
gt_cols <- grep("^GT_", names(variants), value = TRUE)

variants_long <- variants %>%
  pivot_longer(cols = all_of(gt_cols),
               names_to = "sample",
               values_to = "genotype") %>%
  mutate(sample = str_remove(sample, "^GT_"))   # clean sample names

# Encode genotype: 0 = reference, 1 = alternate, . = missing
variants_long <- variants_long %>%
  mutate(genotype = case_when(
    genotype == "0" ~ "ref",
    genotype == "1" ~ "alt",
    genotype == "." ~ "missing",
    TRUE ~ "unknown"
  ))

view(variants_long)

#-------------------------------------------------------------------------------
# 5. Per‑sample genotype counts
#-------------------------------------------------------------------------------
sample_summary <- variants_long %>%
  group_by(sample, genotype) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = count, values_fill = 0) %>%
  mutate(total_variants = ref + alt + missing)

sample_summary

#-------------------------------------------------------------------------------
# 6. Variant sharing patterns
#-------------------------------------------------------------------------------
# For each variant, how many samples carry the alternate allele?
variants_sharing <- variants_long %>%
  filter(genotype == "alt") %>%
  group_by(CHROM, POS, REF, ALT) %>%
  summarise(num_alt_samples = n(),
            alt_samples = paste(sample, collapse = ","), .groups = "drop")

# Count of variants by sharing category
variants_sharing %>%
  count(num_alt_samples) %>%
  arrange(num_alt_samples)

# Variants unique to a single sample
unique_variants <- variants_sharing %>% filter(num_alt_samples == 1)
head(unique_variants)

# Variants present in all 4 samples
fixed_variants <- variants_sharing %>% filter(num_alt_samples == 4)
head(fixed_variants)

#-------------------------------------------------------------------------------
# 7. Visualizations
#-------------------------------------------------------------------------------
# (a) Quality histogram by variant type
ggplot(variants, aes(x = QUAL, fill = variant_type)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Variant Quality Distribution",
       x = "Phred-scaled QUAL score", y = "Count") +
  theme_minimal()

# (b) Depth histogram
ggplot(variants, aes(x = DP)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Read Depth Distribution across Variants",
       x = "Depth (DP)", y = "Count") +
  theme_minimal()

# (c) Boxplot of DP per variant type
ggplot(variants, aes(x = variant_type, y = DP, fill = variant_type)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Depth by Variant Type", x = "", y = "Depth") +
  theme_minimal()

# (d) Number of alternate alleles per sample (bar plot)
alt_counts <- variants_long %>%
  filter(genotype == "alt") %>%
  count(sample, name = "alt_count")

ggplot(alt_counts, aes(x = sample, y = alt_count, fill = sample)) +
  geom_col() +
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "Number of Alternate Alleles per Sample",
       x = "Sample", y = "Alternate allele count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# (e) Heatmap of variant sharing intensity
# Create a presence/absence matrix (variants × samples)
presence_mat <- variants_long %>%
  mutate(value = ifelse(genotype == "alt", 1, 0)) %>%
  select(CHROM, POS, sample, value) %>%
  pivot_wider(names_from = sample, values_from = value, values_fill = 0) %>%
  column_to_rownames(var = "POS")  # using POS as row name (CHROM is same for all)

# Only first 100 variants for legibility (remove if you want all)
presence_mat_sub <- presence_mat[1:min(100, nrow(presence_mat)), ]

# Heatmap with pheatmap (install if needed)
if (!require(pheatmap)) install.packages("pheatmap")
library(pheatmap)
pheatmap(presence_mat_sub,
         main = "Variant Presence/Absence (first 100 variants)",
         color = c("white", "steelblue"),
         legend = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         angle_col = 45)

#-------------------------------------------------------------------------------
# 8. Save key results to files
#-------------------------------------------------------------------------------
write.csv(sample_summary, "sample_variant_summary.csv", row.names = FALSE)
write.csv(variants_sharing, "variant_sharing_summary.csv", row.names = FALSE)

# Save plots if desired
ggsave("QUAL_histogram.png", width = 8, height = 5, dpi = 300)
ggsave("DP_histogram.png", width = 8, height = 5, dpi = 300)
ggsave("alt_alleles_per_sample.png", width = 6, height = 5, dpi = 300)

#-------------------------------------------------------------------------------
# End of script
#-------------------------------------------------------------------------------