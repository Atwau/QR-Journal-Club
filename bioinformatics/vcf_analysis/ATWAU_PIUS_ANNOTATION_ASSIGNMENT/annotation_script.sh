#!/bin/bash

# after installing snpeff package
# We need SnpEff database for which we shall use as a reference to annotate our own varinats
# you can make a custom snpeff database
# but that is not necessary if there is already a well annotated snpeff database
# In the case of MTB, there are already well annotated snpeff databases for reference
# So we proceed to search for existing mtb snpeff databases
snpEff databases | grep -i "tuberculosis_h37rv"

# Look through and select which you will download
# SnpEff will download the database and place it in a location of its choice i.e.
# where the Conda installation of SnpEff expects to find it by default, you do not need to specify
# the data directory when running your annotations.
snpEff download Mycobacterium_tuberculosis_h37rv

#=======================================================================================

# Next is you write a script to do the actual annotation
# Ensure the output directory exists
mkdir -p annotation

# Loop through all filtered VCF files
for vcf in vcf/filtered/*.filtered.vcf; do
    # Extract the base filename
    basename=$(basename "$vcf" .filtered.vcf)

    echo "Annotating $basename..."

    # Replace the chromosome name and pipe to snpEff
    sed 's/NC_000962.3/Chromosome/g' "$vcf" | \
    snpEff eff Mycobacterium_tuberculosis_h37rv - \
    > "annotation/${basename}.ann.vcf"

    echo "Finished $basename"
done

echo "All files successfully annotated."

#=======================================================================================

# After annotation, you can view the annotated vcf files
# you can also decide to extract only a few columns from the annotated vcf files
# and to be more precise, you can decide to extract only primary impact variants

# 1. Extracting only relevant columns with annotation information
#---------------------------------------------------------------------------------------

mkdir -p r_analysis/annotated

# Loop through all annotated VCF files in the annotation directory
for vcf in annotation/*.ann.vcf; do
    # Extract the base filename (e.g., SRR28515580)
    basename=$(basename "$vcf" .ann.vcf)

    echo "Extracting fields for $basename..."

    # Run SnpSift with -s directive
    # The -s ',' flag forces multiple values to stay in their assigned column
    # This is very useful as earlier on I had faced challenge of distorted output values
    # For example the affected genes were not visible anywhere in the output .tsv file
    SnpSift extractFields -s ',' "$vcf" \
        CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_P" \
        > "r_analysis/annotated/${basename}_extracted.tsv"

    echo "Finished $basename"
done

echo "All extractions complete. Files are in r_analysis/annotated/"


# 2. Extracting only primary impact variants
#----------------------------------------------------------------------------------------

# Create the output directory
mkdir -p r_analysis/annotated

# Loop through all annotated VCF files in the annotation directory
for vcf in annotation/*.ann.vcf; do
    # Extract the base filename (e.g., SRR28515580)
    basename=$(basename "$vcf" .ann.vcf)

    echo "Extracting primary fields for $basename..."

    # Changed ANN[*] to ANN[0] to isolate the single highest-impact consequence.
    # The -s flag is removed because we are no longer processing multi-value lists.
    SnpSift extractFields "$vcf" \
        CHROM POS REF ALT "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].HGVS_P" \
        > "r_analysis/annotated/${basename}_extracted_1st_impact.tsv"

    echo "Finished $basename"
done

echo "All extractions complete. Files are in r_analysis/annotated/"

#=======================================================================================

# The rest of the analysis can be done in R

#=======================================================================================

# MTB PROFILINING 

#!/bin/bash

# Loop through the rest of samples with *.fastq.gz
for sample in SRR28515580 SRR35497378 SRR35891567; do
    echo "Processing lineage for sample: $sample"
    tb-profiler profile \
        -1 data/trimmed/${sample}_1_paired.fastq.gz \
        -2 data/trimmed/${sample}_2_paired.fastq.gz \
        --dir tb_profiling_results/trimmed \
        -p $sample \
        -t 4 \
        --txt
done

# Explicitly process SRR34502717
echo "Processing lineage for sample: SRR34502717"
tb-profiler profile \
    -1 data/trimmed/SRR34502717_1_paired.fq.gz \
    -2 data/trimmed/SRR34502717_2_paired.fq.gz \
    --dir tb_profiling_results/trimmed \
    -p SRR34502717 \
    -t 4 \
    --txt


#======================================================================================

# PHYLOGENETIC TREE CONSTRUCTION

