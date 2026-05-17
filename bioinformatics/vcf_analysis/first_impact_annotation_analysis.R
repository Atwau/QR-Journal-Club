# ==============================================================================
# Single Sample MTB Variant Analysis Script
# ==============================================================================

# 1. Load the required library
# ------------------------------------------------------------------------------
# The tidyverse package makes reading and filtering data much easier.
library(tidyverse)
library(flextable)

# 2. Read the specific CSV file
# ------------------------------------------------------------------------------
# Read the data into R
srr285_primary_impact <- read_tsv("mtb_assignment/SRR28515580_extracted_1st_impact.tsv")
srr345_primary_impact <- read_tsv("mtb_assignment/SRR34502717_extracted_1st_impact.tsv")
srr354_primary_impact <- read_tsv("mtb_assignment/SRR35497378_extracted_1st_impact.tsv")
srr358_primary_impact <- read_tsv("mtb_assignment/SRR35891567_extracted_1st_impact.tsv")

spec(srr285_primary_impact)
view(srr285_primary_impact)


# ------------------------------------------------------------------------------
# Clean the Column Names
# Remove the messy "ANN[*]." prefix from SnpSift to make columns easier to type
srr285_primary_impact <- srr285_primary_impact %>%
  rename(
    EFFECT = `ANN[0].EFFECT`,
    IMPACT = `ANN[0].IMPACT`,
    GENE   = `ANN[0].GENE`,
    HGVS_P = `ANN[0].HGVS_P`
  )
glimpse(srr285_primary_impact)

# SRR345
# Remove the messy "ANN[*]." prefix from SnpSift to make columns easier to type
srr345_primary_impact <- srr345_primary_impact %>%
  rename(
    EFFECT = `ANN[0].EFFECT`,
    IMPACT = `ANN[0].IMPACT`,
    GENE   = `ANN[0].GENE`,
    HGVS_P = `ANN[0].HGVS_P`
  )
glimpse(srr345_primary_impact)

# SRR354
# Remove the messy "ANN[*]." prefix from SnpSift to make columns easier to type
srr354_primary_impact <- srr354_primary_impact %>%
  rename(
    EFFECT = `ANN[0].EFFECT`,
    IMPACT = `ANN[0].IMPACT`,
    GENE   = `ANN[0].GENE`,
    HGVS_P = `ANN[0].HGVS_P`
  )
glimpse(srr354_primary_impact)

# SRR358
# Remove the messy "ANN[*]." prefix from SnpSift to make columns easier to type
srr358_primary_impact <- srr358_primary_impact %>%
  rename(
    EFFECT = `ANN[0].EFFECT`,
    IMPACT = `ANN[0].IMPACT`,
    GENE   = `ANN[0].GENE`,
    HGVS_P = `ANN[0].HGVS_P`
  )
glimpse(srr358_primary_impact)

#===============================================================================

# Analysing each sample individually


# 4. Basic Overview of the total number of variants in each sample

#Total number of variants in SRR285
print(nrow(srr285_primary_impact))

# Variant impact count (i.e.HIGH, MODERATE, LOW, or MODIFIER)
srr285_impact_summary <- srr285_primary_impact %>% count(IMPACT, name = "Count")
print(srr285_impact_summary)

# exporting output to tables in word
srr285_impact_summary <- as.data.frame(srr285_impact_summary)
srr285_impact_summary %>%
  flextable() %>%
  autofit() %>%
  set_caption("Summary of variant impact for SRR288* Sample") %>% 
  save_as_docx(path = "mtb_assignment/srr285_impact_summary.docx")

#-------------------------------------------------------------------------------
#Total number of variants in SRR345
print(nrow(srr345_primary_impact))

# Variant impact count (i.e.HIGH, MODERATE, LOW, or MODIFIER)
srr345_impact_summary <- srr345_primary_impact %>% count(IMPACT, name = "Count")
print(srr345_impact_summary)

# exporting output to tables in word
srr345_impact_summary <- as.data.frame(srr345_impact_summary)
srr345_impact_summary %>%
  flextable() %>%
  autofit() %>%
  set_caption("Summary of variant impact for SRR345* Sample") %>% 
  save_as_docx(path = "mtb_assignment/srr345_impact_summary.docx")

#-------------------------------------------------------------------------------
#Total number of variants in SRR354
print(nrow(srr354_primary_impact))

# Variant impact count (i.e.HIGH, MODERATE, LOW, or MODIFIER)
srr354_impact_summary <- srr354_primary_impact %>% count(IMPACT, name = "Count")
print(srr354_impact_summary)

# exporting output to tables in word
srr354_impact_summary <- as.data.frame(srr354_impact_summary)
srr354_impact_summary %>%
  flextable() %>%
  autofit() %>%
  set_caption("Summary of variant impact for SRR354* Sample") %>% 
  save_as_docx(path = "mtb_assignment/srr354_impact_summary.docx")

#-------------------------------------------------------------------------------
#Total number of variants in SRR358
print(nrow(srr358_primary_impact))

# Variant impact count (i.e.HIGH, MODERATE, LOW, or MODIFIER)
srr358_impact_summary <- srr358_primary_impact %>% count(IMPACT, name = "Count")
print(srr358_impact_summary)

# exporting output to tables in word
srr358_impact_summary <- as.data.frame(srr358_impact_summary)
srr358_impact_summary %>%
  flextable() %>%
  autofit() %>%
  set_caption("Summary of variant impact for SRR358* Sample") %>% 
  save_as_docx(path = "mtb_assignment/srr358_impact_summary.docx")


#===============================================================================

# Combining Data into One Large Dataset
# ------------------------------------------------------------------------------
all_samples <- bind_rows(
  "SRR285" = srr285_primary_impact,
  "SRR345" = srr345_primary_impact,
  "SRR354" = srr354_primary_impact,
  "SRR358" = srr358_primary_impact,
  .id = "Sample"
)

all_samples
view(all_samples)

# ------------------------------------------------------------------------------
# Generate Combined Overview Table and Export to Word
# Print total variants per sample to the console
total_variants <- all_samples %>% count(Sample, name = "Total_Variants")
print(total_variants)

# Create the combined cross-tabulated impact summary table
combined_impact_summary <- all_samples %>%
  count(Sample, IMPACT) %>%
  pivot_wider(names_from = Sample, values_from = n, values_fill = 0)

# View the table in the console
print(combined_impact_summary)

# Export the combined output table to a single Word document
combined_impact_summary %>%
  as.data.frame() %>%
  flextable() %>%
  autofit() %>%
  set_caption("Combined Summary of Variant Impacts Across MTB Samples") %>% 
  save_as_docx(path = "mtb_assignment/combined_impact_summary.docx")


# ------------------------------------------------------------------------------
# Filter for Significant Mutations

# Keep only the mutations that are likely to change how the protein works
# (HIGH = Frameshifts/Stop codons, MODERATE = Missense amino acid changes)
all_significant_variants <- all_samples %>%
  filter(IMPACT %in% c("HIGH", "MODERATE"))

# Number of significant (High/Moderate) variants per sample
# Instead of a single total row count, this counts significant variants for each sample
all_significant_counts <- all_significant_variants %>% 
  count(Sample, name = "Significant_Variants")
print(all_significant_counts)


# ------------------------------------------------------------------------------
# Checking for Drug Resistance Mutations
# Defining the most important MTB resistance genes
dr_genes <- c("rpoB", "katG", "inhA", "gyrA", "gyrB", "embB", "pncA", "rpsL")

# Search the significant variants for hits in the drug resistance genes
all_dr_mutations <- all_significant_variants %>%
  filter(GENE %in% dr_genes) %>%
  # Included 'Sample' in the selection to keep track of the source sample
  select(Sample, CHROM, POS, GENE, EFFECT, IMPACT, HGVS_P) %>%
  # Sorted by Sample first, then by Gene name
  arrange(Sample, GENE)

print(all_dr_mutations)


# ------------------------------------------------------------------------------
# Saving the Combined Results
# Save the combined drug resistance hits to a single CSV file
write_csv(all_dr_mutations, "mtb_assignment/combined_drug_resistance.csv")

# Exporting combined output to a single Word document table
all_dr_mutations_df <- as.data.frame(all_dr_mutations)
all_dr_mutations_df %>%
  flextable() %>%
  autofit() %>%
  set_caption("Summary of Drug Resistance Mutation Genes Across All MTB Samples") %>% 
  save_as_docx(path = "mtb_assignment/combined_drug_resistance.docx")

# ==============================================================================

# Summary of variant types in the most clinically relevant isolate 
# Define the drug resistance genes
resistance_genes <- c("rpoB", "katG", "inhA", "gyrA", "gyrB", "embB", "pncA", "rpsL")

# 2. Identify which sample has the highest number of drug resistance mutations
top_isolate_info <- all_samples %>%
  filter(IMPACT %in% c("HIGH", "MODERATE"), GENE %in% resistance_genes) %>%
  count(Sample, name = "DR_Mutation_Count") %>%
  arrange(desc(DR_Mutation_Count)) %>%
  slice_head(n=1) # To extract the sample with the maximum count (1 for 1 sample)
# to select top two samples, then we adjust slice_head(n=1) to slice_head(n=2)

top_isolate_info

# Store the name of the top sample (e.g., "SRR285")
most_relevant_sample <- top_isolate_info$Sample

cat("\n--- Most Clinically Relevant Isolate Identified: --- \n", most_relevant_sample, "\n")

# Filtering the complete dataset for only this specific isolate
relevant_isolate_data <- all_samples %>%
  filter(Sample == most_relevant_sample)

# Summarizing the distinct variant types (EFFECT) within this isolate
variant_type_summary <- relevant_isolate_data %>%
  count(EFFECT, IMPACT, name = "Count") %>%
  arrange(desc(Count))

# Printing the summary table to the console
print(variant_type_summary)

# Exporting the summary table to a Word document
variant_type_summary %>%
  as.data.frame() %>%
  flextable() %>%
  autofit() %>%
  set_caption(paste("Summary of Variant Functional Types in the Most Clinically Relevant Isolate:", most_relevant_sample)) %>% 
  save_as_docx(path = "mtb_assignment/most_relevant_isolate_variant_summary.docx")

