# Install vcfR if not already installed
# install.packages("vcfR")

# Load the required library
library(vcfR)

# 1. Read the VCF file into the R environment
# This parses the meta-information, fixed fields, and genotype fields into a vcfR object
vcf_file <- "D:/teaching/QR Journal Club/bioinformatics/vcf_analysis/variants.raw.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)

# 2. Inspect the vcfR object structure
show(vcf)

# 3. Extract the Fixed Fields (CHROM, POS, REF, ALT, QUAL, FILTER)
# getFIX() creates a matrix of the fixed information
fixed_data <- getFIX(vcf)
head(fixed_data)

# 4. Extract Genotype (GT) data
# extract.gt() extracts specific FORMAT fields. Here, extracting the genotype calls.
gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
head(gt_matrix)

# 5. Extract Read Depth (DP) per sample
# Converts the DP field into a numeric matrix for quantitative analysis
dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
head(dp_matrix)

# 6. Convert to a chromR object for quality control visualization
# The chromR object aggregates the data for plotting
chrom <- create.chromR(name="My_Variants", vcf=vcf, verbose=FALSE)

# Calculate summary statistics for the variants
chrom <- proc.chromR(chrom, verbose=FALSE)


# Generate standard quality plots (Read Depth, Mapping Quality, etc.)
plot(chrom)

# Plot the distribution of variants across the genome position
chromoqc(chrom)
