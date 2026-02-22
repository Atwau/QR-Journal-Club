---
title: "pca_plot"
output: html_document
---

# PCA Plots in DESeq2 Analysis

In the context of **DESeq2** (and transcriptomics generally), a **Principal Component Analysis (PCA)** plot is a fundamental quality control (QC) tool used to visualize the overall similarity between samples. It reduces the complexity of high-dimensional gene expression data into a 2D or 3D coordinate system.


## Primary Uses of PCA in DESeq2

### 1. Assessment of Biological Replicates
The most immediate use is to verify that your **biological replicates** cluster together. If samples from the "Treatment" group are scattered across the plot or overlap significantly with the "Control" group, it suggests high intra-group variability, which reduces the power to detect differentially expressed genes (DEGs).

### 2. Identification of Batch Effects
PCA is highly effective at revealing systematic biases. If samples cluster by technical factors—such as the date of sequencing, the technician, or the flow cell—rather than by biological condition, a **batch effect** is present.
* **Correction:** If a batch effect is visible, you should include the batch as a covariate in your DESeq2 design formula: `~ batch + condition`.

### 3. Outlier Detection
A PCA plot allows you to identify "rogue" samples that sit far away from their designated group. Common causes for outliers include:
* Library preparation failure.
* Sample contamination.
* Mislabeling (e.g., a "Control" sample clustering perfectly with the "Treatment" group).


## Technical Mechanics

In a standard DESeq2 workflow, PCA is performed on **transformed data** (using `vst` or `rlog`) rather than raw counts to stabilize variance across the dynamic range of expression.

* **PC1 (Principal Component 1):** The axis that explains the highest percentage of variance in the dataset. Ideally, this represents your primary experimental variable (e.g., Treatment vs. Control).
* **PC2 (Principal Component 2):** The axis explaining the second-highest variance, orthogonal to PC1. This often captures secondary biological effects or technical noise.


## Limitations

While useful, PCA has specific boundaries:

* **Linearity:** PCA assumes linear relationships between variables. If the underlying biology is highly non-linear, other methods like **t-SNE** or **UMAP** might be used for visualization, though PCA remains the gold standard for QC.
* **Top Variance Only:** By default, the `plotPCA` function in DESeq2 uses only the **top 500 most variable genes**. This is usually sufficient to capture the primary "signal," but it does not represent the entire transcriptome.

