# ==============================================================================
# Single Sample MTB Variant Analysis Script
# ==============================================================================

# 1. Load the required library
# ------------------------------------------------------------------------------
# The tidyverse package makes reading and filtering data much easier.
library(tidyverse)


# ------------------------------------------------------------------------------
# Read the data into R
srr28_data <- read_tsv("mtb_assignment/SRR28515580_extracted.tsv")
spec(srr28_data)
view(srr28_data)

# 3. Clean the Column Names
# ------------------------------------------------------------------------------
# Remove the messy "ANN[*]." prefix from SnpSift to make columns easier to type
clean_srr28 <- srr28_data %>%
  rename(
    EFFECT = `ANN[*].EFFECT`,
    IMPACT = `ANN[*].IMPACT`,
    GENE   = `ANN[*].GENE`,
    HGVS_P = `ANN[*].HGVS_P`
  )

glimpse(clean_srr28)
view(clean_srr28)

unique(clean_srr28$EFFECT)
unique(clean_srr28$IMPACT)

# 4. Basic Overview of the Sample
# ------------------------------------------------------------------------------
cat("\n--- Total number of variants in srr28 ---\n")
print(nrow(clean_srr28))

cat("\n--- Breakdown of variants by their IMPACT ---\n")
# This shows how many mutations are HIGH, MODERATE, LOW, or MODIFIER
impact_summary <- clean_srr28 %>% count(IMPACT, name = "Count")
print(impact_summary)
view(impact_summary)

# 5. Filter for Significant Mutations
# ------------------------------------------------------------------------------
# Keep only the mutations that are likely to change how the protein works
# (HIGH = Frameshifts/Stop codons, MODERATE = Missense amino acid changes)
unique(clean_srr28$IMPACT)
significant_variants <- clean_srr28 %>%
  filter(str_detect(IMPACT, "HIGH|MODERATE"))
view(significant_variants)

cat("\n--- Number of significant (High/Moderate) variants ---\n")
print(nrow(significant_variants))

# 6. Check for Drug Resistance Mutations
# ------------------------------------------------------------------------------
# 1. Define the resistance genes as a single search pattern separated by OR (|)
dr_genes_pattern <- "rpoB|katG|inhA|gyrA|gyrB|embB|pncA|rpsL"

# 2. Use str_detect to find any of these genes inside the comma-separated lists
dr_mutations <- significant_variants %>%
  filter(str_detect(GENE, dr_genes_pattern)) %>%
  select(CHROM, POS, GENE, EFFECT, IMPACT, HGVS_P) %>%
  arrange(GENE)

view(dr_mutations)
# 3. Print the results
print(dr_mutations)

# 7. Save the Results
# ------------------------------------------------------------------------------
# Save the drug resistance hits to a new file for later use
write_csv(dr_mutations, "mtb_assignment/SRR28515580_drug_resistance.csv")

cat("\nDone! Drug resistance summary saved to r_analysis/SRR28515580_drug_resistance.csv\n")