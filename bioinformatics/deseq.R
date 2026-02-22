# 1. Load necessary libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

# 2. Set working directory
# Replace with the path to the directory containing your input files
getwd()
setwd("/home/.../") 

# 3. Read in the count data and sample information
counts_table <- read.csv("bioinformatics/raw_counts.tsv", sep="\t", row.names=1)
sample_info <- read.csv("bioinformatics/design.tsv", sep="\t", row.names=1)

# 4. Set factor levels for the groups
factors <- factor(sample_info$group)
groups <- unique(sample_info$group)
groups <- rev(groups) # Reverse to ensure the control group is set as the reference
sample_info$group <- factors

# 5. Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_table, 
                              colData = sample_info, 
                              design = ~ group)

# 6. Set the reference level for the group factor to 'control'
dds$group <- relevel(dds$group, ref = "control")

# 7. Filter out genes with low counts 
# Keep genes with counts >= 10 in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# 8. Run the DESeq2 analysis pipeline
dds <- DESeq(dds, test="Wald", sfType="poscounts")

# 9. Extract and format the results
deseq_results <- results(dds)
deseq_results <- as.data.frame(deseq_results)

# Add a column for gene names using the row names
deseq_results$gene_name <- rownames(deseq_results)

# Reorder the columns to bring gene_name first
deseq_results <- subset(deseq_results, select = c("gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))

# Save the complete results to a file
write.table(deseq_results, file="deseq_results.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# 10. Filter for significantly differentially expressed genes
# Criteria: adjusted p-value < 0.05 and absolute log2FoldChange >= 1
sig_genes <- subset(deseq_results, padj < 0.05 & abs(log2FoldChange) >= 1)

# Order the significant genes by their adjusted p-value
sig_genes <- sig_genes[order(sig_genes$padj), ]

# Save the significant genes to an output file
write.table(sig_genes, file="significant_genes.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# --- Plotting Visualizations ---

# 11. Histogram of adjusted P-values
hist(deseq_results$padj, breaks=seq(0, 1, length=21), main="Histogram of padj", xlab="padj")

# 12. Volcano Plot
# Set plot margins and symbol size
par(mar=c(4, 4, 4, 1), cex=1.5) 

# Plot all genes
plot(x = deseq_results$log2FoldChange, 
     y = -log10(deseq_results$padj), 
     main = "Volcano Plot", 
     xlab = "log2 Fold Change", 
     ylab = "-log10(padj)", 
     pch = 20, 
     cex = 0.5)

# Highlight significant genes (Downregulated in blue, Upregulated in red)
points(x = sig_genes$log2FoldChange, 
       y = -log10(sig_genes$padj), 
       pch = 20, 
       col = ifelse(sig_genes$log2FoldChange > 0, "red", "blue"), 
       cex = 1)

# Add a legend for the highlighted genes
legend("topleft", legend=c("Down", "Up"), pch=20, col=c("blue", "red"))

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