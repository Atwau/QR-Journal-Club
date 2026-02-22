#===============================================================================
# DESeq analysis for RNA-seq data
#===============================================================================

# Reference materials
# Original study:https://www.science.org/doi/10.1126/scisignal.adf1947
# Raw data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203159

#-------------------------------------------------------------------------------


# 1. Load necessary libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

# 2. Set working directory
# Replace with the path to the directory containing your input files
getwd()
setwd("/home/.../") # not necessary since I am working in R project

# 3. Read in the count data and sample information
counts_table <- read.csv("bioinformatics/raw_counts.tsv", sep="\t", row.names=1)
sample_info <- read.csv("bioinformatics/design.tsv", sep="\t", row.names=1)

# Inspecting raw counts and sample info imports
glimpse(counts_table)
head(counts_table)
glimpse(sample_info)

view(counts_table)
view(sample_info)


# 4. Set factor levels for the groups
sample_info$group <- factor(sample_info$group, levels = c("control", "tgf-beta"))

# 5. Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_table, 
                              colData = sample_info, 
                              design = ~ group)


#-------------------------------------------------------------------------------

# 6. Filter out genes with low counts 
# Keep genes with counts >= 10 in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

dds

# 8. Run the DESeq2 analysis pipeline
dds <- DESeq(dds, test="Wald", sfType="poscounts")
?DESeq

# 9. Extract and format the results
deseq_results <- results(dds)
deseq_results <- as.data.frame(deseq_results)

view(deseq_results)

# Add a column for gene names using the row names
deseq_results$gene_name <- rownames(deseq_results)

view(deseq_results)

# Reorder the columns to bring gene_name first
deseq_results <- subset(deseq_results, select = c("gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))

view(deseq_results)

# Save the complete results to a file
write.table(deseq_results, file="bioinformatics/deseq_results.tsv", sep="\t", quote=FALSE, row.names=FALSE)

#-------------------------------------------------------------------------------

# 10. Filter for significantly differentially expressed genes
# Criteria: adjusted p-value < 0.05 and absolute log2FoldChange >= 1
sig_genes <- subset(deseq_results, padj < 0.05 & abs(log2FoldChange) >= 1)

# Order the significant genes by their adjusted p-value
sig_genes <- sig_genes[order(sig_genes$padj), ]

# Save the significant genes to an output file
write.table(sig_genes, file="bioinformatics/significant_genes.tsv", sep="\t", quote=FALSE, row.names=FALSE)


#===============================================================================
# RESULTS VISUALIZATION
#===============================================================================

# 11. Histogram of adjusted P-values
hist(deseq_results$padj, breaks=seq(0, 1, length=21), 
     main="Histogram of padj", xlab="padj")


#-------------------------------------------------------------------------------

# 12. Volcano Plot
# Set clean, standard margins (Bottom, Left, Top, Right)
par(mar=c(5, 5, 4, 2), cex=1.0) 


# Create the base plot
plot(x = deseq_results$log2FoldChange, 
     y = -log10(deseq_results$padj), 
     main = "Volcano Plot", 
     xlab = "log2 Fold Change", 
     ylab = "-log10(padj)", 
     pch = 20, 
     col = "gray40", # Light gray for non-significant points
     cex = 0.5)

# Add colored points for significance
points(x = sig_genes$log2FoldChange, 
       y = -log10(sig_genes$padj), 
       pch = 20, 
       col = ifelse(sig_genes$log2FoldChange > 0, "red", "blue"), 
       cex = 0.7)

# Add the legend (using bty="n" ensures no box covers the dots)
legend("topright", 
       legend = c("Downregulated", "Upregulated"), 
       pch = 20, 
       col = c("blue", "red"), 
       bty = "n", 
       cex = 0.4,
       pt.cex = 1.2)

#-------------------------------------------------------------------------------

# 13. Extract Normalized Counts and Plot Heatmap
# Extract normalized counts
normalized_counts <- counts(dds, normalized=TRUE)

# Perform Variance Stabilizing Transformation (VST) for visualization
vsd <- vst(dds, blind=FALSE)

# Select the top 10 genes with the smallest adjusted p-values
top10_genes <- head(sig_genes$gene_name, 10)

# Extract log-transformed values for the top 10 genes from the VST object
mat <- assay(vsd)[top10_genes, ]

# Plot the heatmap with row and column clustering
pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE)













