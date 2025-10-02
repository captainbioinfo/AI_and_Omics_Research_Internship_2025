#-------------
# Assignment 
#-------------

# Work through the full preprocessing workflow with your own dataset

# 1. Perform quality control before and after normalization and 
# check whether any arrays are flagged as outliers. 
# note down how many you found before and after normalization

# 2. Normalize the data and then apply filtering to remove low-intensity probes 
# and note how many transcripts remain. 

# 3. Use the phenotype information to define your target groups and re-label them (e.g normal vs cancer)
##############################################################################################################


#######################################################################
#### 0. Install and Load Required Packages ####
#######################################################################
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("GEOquery","lumi","arrayQualityMetrics"))

library(GEOquery)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(lumi)
library(Biobase)
library(arrayQualityMetrics)

#########################################################################################
# 1. Download expression matrix + metadata via series matrix
gse <- getGEO("GSE156451", GSEMatrix = TRUE)
eset <- gse[[1]]

expr_mat <- exprs(eset)          # processed / normalized values (series matrix)
pheno <- pData(eset)
head(eset)
# 2. If you have raw counts (from supplementary files)
# counts = raw counts matrix (genes × samples)

counts <- read.delim("GSE156451_raw_counts_GRCh38.p13_NCBI.tsv", 
                     row.names = 1, check.names = FALSE)

dim(counts)        # genes × samples
head(counts[,1:5]) # preview first 5 samples


#  Check sample overlap
colnames(counts)[1:5]       # sample names from counts
rownames(pheno)[1:5]        # sample names from pheno

#  Subset pheno to only include samples in counts
common_samples <- intersect(colnames(counts), rownames(pheno))

counts <- counts[, common_samples]
pheno  <- pheno[common_samples, ]

#  Rebuild ExpressionSet
exprs_eset <- ExpressionSet(
  assayData = as.matrix(log2(counts + 1)),
  phenoData = AnnotatedDataFrame(pheno)
)

# Run QC
arrayQualityMetrics(expressionset = exprs_eset,
                    outdir = "QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = FALSE)  # already log2

#  Build ExpressionSet from RNA-seq counts
exprs_eset <- ExpressionSet(assayData = as.matrix(log2(counts + 1)), 
                            phenoData = AnnotatedDataFrame(pheno))

#  Run arrayQualityMetrics QC report
arrayQualityMetrics(expressionset = exprs_eset,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = FALSE)  # already log2-transformed

#########################################################################
# 2. Normalize the data (DESeq2 variance stabilizing transformation)
#########################################################################

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = pheno,
                              design    = ~ 1)  # no condition yet


# Filter low counts (keep genes with at least 10 reads in ≥ 3 samples)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]   ### I just dont want to loose feature for ML step that is why I keep it low.
cat("Genes after filtering:", nrow(dds), "\n")
# Genes after filtering: 23370
# Variance stabilizing transformation (VST normalization)
vsd <- vst(dds, blind = TRUE)
norm_counts <- assay(vsd)   # normalized expression matrix

#########################################################################
# Define phenotype groups (Normal vs Cancer, etc.)
#########################################################################

# Inspect phenotype metadata
colnames(pheno)

# Example: use tissue:ch1 column (you’ll need to check which column has normal/tumor info)
table(pheno$`tissue:ch1`)

# Correctly create condition based on actual labels
pheno$condition <- ifelse(pheno$`tissue:ch1` == "Native tissue",
                          "Normal",
                          "Cancer")

# Check
table(pheno$condition)

# Relabel in DESeq2 object
colData(dds)$condition <- pheno$condition

#########################################################################
# 4. Quick plots
#########################################################################
library(ggplot2)
library(DESeq2)
library(pheatmap)

# ------------------------
# PCA Plot
# ------------------------
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  scale_color_manual(values = c("Normal" = "blue", "Cancer" = "red")) +
  ggtitle("PCA of VST-normalized counts")

# Save PCA plot
ggsave("PCA_plot.png", pca_plot, width = 6, height = 5)
ggsave("PCA_plot.pdf", pca_plot, width = 6, height = 5)

# ------------------------
# Sample Distance Heatmap
# ------------------------
# Calculate Euclidean distances between samples
sample_dists <- dist(t(norm_counts))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(norm_counts)
colnames(sample_dist_matrix) <- colnames(norm_counts)

# Heatmap of sample distances
pheatmap(sample_dist_matrix,
         annotation_col = pheno["condition"],
         main = "Sample-to-sample distances",
         filename = "Sample_Distance_Heatmap.png",
         width = 6, height = 5)

# Save as PDF too
pheatmap(sample_dist_matrix,
         annotation_col = pheno["condition"],
         main = "Sample-to-sample distances",
         filename = "Sample_Distance_Heatmap.pdf",
         width = 6, height = 5)
#####################################################################
# Perform PCA on VST-normalized counts
library(ggplot2)
library(pheatmap)
library(DESeq2)

# ------------------------
# 1. PCA Plot
# ------------------------
# Perform PCA on VST-normalized counts
vsd_mat <- assay(vsd)
pca <- prcomp(t(vsd_mat), scale. = TRUE)

# Percentage of variance explained
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

# Create PCA dataframe
pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     condition = pheno$condition,
                     sample = colnames(vsd_mat))

# PCA plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  scale_color_manual(values = c("Normal" = "blue", "Cancer" = "red")) +
  ggtitle("PCA of VST-normalized counts")

# Save PCA plot
ggsave("PCA_plot.png", pca_plot, width = 6, height = 5)
ggsave("PCA_plot.pdf", pca_plot, width = 6, height = 5)

# ------------------------
# 2. Sample Distance Heatmap
# ------------------------
# Calculate Euclidean distances between samples
sample_dists <- dist(t(vsd_mat))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(vsd_mat)
colnames(sample_dist_matrix) <- colnames(vsd_mat)

# Heatmap of sample distances
pheatmap(sample_dist_matrix,
         annotation_col = pheno["condition"],
         main = "Sample-to-sample distances",
         filename = "Sample_Distance_Heatmap.png",
         width = 6, height = 5)

# Save as PDF too
pheatmap(sample_dist_matrix,
         annotation_col = pheno["condition"],
         main = "Sample-to-sample distances",
         filename = "Sample_Distance_Heatmap.pdf",
         width = 6, height = 5)
