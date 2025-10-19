#--------------------
#### Assignment #### 
#--------------------

# Map probe IDs to gene symbols using AnnotationDbi
#        i.  Find the appropriate annotation package for your platform and use it 
#        ii. Check how many probes map to the same gene and handle duplicates

# Perform differential gene expression analysis using the Limma package
# Create a volcano plot showing upregulated and downregulated genes
# Generate a heatmap of the top 25 DEGs
# Save DEG results (complete, upregulated, downregulated) as CSV files
# Export both plots as PNG images in the Results folder

# Write a short result summary (4–5 lines) explaining 
#        i.   how multiple probes can map to the same gene, how you handled duplicate probes,
#        ii.  Which contrast or comparison did you perform (e.g cancer_vs_normal, diabetes_vs_normal etc)
#        iii. and summarize how many genes were upregulated and downregulated based on your DEG results.
################################# Setup ####################################

# Install required packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("GEOquery", "DESeq2", "org.Hs.eg.db", "AnnotationDbi", 
                       "pheatmap", "ggplot2", "Biobase", "arrayQualityMetrics", 
                       "apeglm", "limma")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# Load libraries
library(GEOquery)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(ggplot2)
library(Biobase)
library(arrayQualityMetrics)
library(apeglm)
library(limma)

################################# Data Loading ####################################

# 1. Download expression matrix + metadata via series matrix
message("Downloading GEO dataset...")
gse <- getGEO("GSE156451", GSEMatrix = TRUE)
eset <- gse[[1]]

expr_mat <- exprs(eset)
pheno <- pData(eset)

# 2. Load raw counts
message("Loading raw counts...")
counts <- read.delim("GSE156451_raw_counts_GRCh38.p13_NCBI.tsv",
                     row.names = 1, check.names = FALSE)

# Align samples
common_samples <- intersect(colnames(counts), rownames(pheno))
counts <- counts[, common_samples]
pheno <- pheno[common_samples, ]

message(paste("Total samples:", ncol(counts)))
message(paste("Total genes before filtering:", nrow(counts)))

################################# Quality Control ####################################

# 3. Build ExpressionSet for QC
message("Running quality control...")
exprs_eset <- ExpressionSet(
  assayData = as.matrix(log2(counts + 1)),
  phenoData = AnnotatedDataFrame(pheno)
)

if (!dir.exists("Results")) dir.create("Results")
if (!dir.exists("Results/QC_Raw_Data")) dir.create("Results/QC_Raw_Data", recursive = TRUE)

arrayQualityMetrics(expressionset = exprs_eset,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = FALSE)

################################# Gene ID Mapping ####################################

# 4. Map Ensembl IDs to Gene Symbols
message("Mapping gene IDs to symbols...")

# Clean Ensembl IDs: remove version numbers (.1, .2, etc)
ensembl_ids <- gsub("\\.\\d+$", "", rownames(counts))

# Map to Gene Symbols using org.Hs.eg.db
# Note: For RNA-seq data with Ensembl IDs, we use ENSEMBL as keytype
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",  # Here we can use "ENSEMBL" also
                       multiVals = "first")

# Create mapping data frame
gene_mapping <- data.frame(
  EnsemblID = ensembl_ids,
  GeneSymbol = gene_symbols,
  stringsAsFactors = FALSE
)

# Check for duplicate mappings (multiple Ensembl IDs to same gene symbol)
duplicate_genes <- gene_mapping$GeneSymbol[duplicated(gene_mapping$GeneSymbol) & !is.na(gene_mapping$GeneSymbol)]
n_duplicates <- length(unique(duplicate_genes))

message(paste("Number of genes with multiple Ensembl IDs:", n_duplicates))

# Handle duplicates: Keep the one with highest mean expression
if (n_duplicates > 0) {
  message("Handling duplicate gene mappings by selecting highest mean expression...")
  
  # Calculate mean expression for each gene
  mean_expr <- rowMeans(counts)
  gene_mapping$MeanExpr <- mean_expr
  
  # For each duplicated gene symbol, keep only the Ensembl ID with highest expression
  gene_mapping <- gene_mapping[order(gene_mapping$GeneSymbol, -gene_mapping$MeanExpr), ]
  gene_mapping <- gene_mapping[!duplicated(gene_mapping$GeneSymbol) | is.na(gene_mapping$GeneSymbol), ]
  
  # Update counts matrix
  counts <- counts[gene_mapping$EnsemblID, ]
  rownames(counts) <- ifelse(is.na(gene_mapping$GeneSymbol), 
                             gene_mapping$EnsemblID, 
                             gene_mapping$GeneSymbol)
}

message(paste("Genes after handling duplicates:", nrow(counts)))

################################# Normalization ####################################

# 5. Define groups (Normal vs Cancer)
message("Defining experimental groups...")
pheno$condition <- ifelse(pheno$`tissue:ch1` == "Native tissue", "Normal", "Cancer")
pheno$condition <- factor(pheno$condition, levels = c("Normal", "Cancer"))

message(table(pheno$condition))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = pheno,
                              design = ~ condition)

# Filter low counts (genes with at least 10 reads in at least 3 samples)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
message(paste("Genes after filtering low counts:", nrow(dds)))

# Normalize using variance stabilizing transformation
message("Normalizing data...")
vsd <- vst(dds, blind = TRUE)
norm_counts <- assay(vsd)

################################# Differential Expression Analysis ####################################

# 6. Run DESeq2 analysis
message("Running differential expression analysis...")
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", "Cancer", "Normal"))

# Apply log fold change shrinkage
res <- lfcShrink(dds, contrast = c("condition", "Cancer", "Normal"), type = "apeglm")

# Add gene symbols to results
res$GeneSymbol <- rownames(res)

# Remove NA padj values
res <- res[!is.na(res$padj), ]

message(paste("Total genes in results:", nrow(res)))

################################# Categorize DEGs ####################################

# 7. Define thresholds and categorize genes
logfc_cutoff <- 1
adjp_cutoff <- 0.05

res$Regulation <- ifelse(res$log2FoldChange > logfc_cutoff & res$padj <= adjp_cutoff, "Upregulated",
                         ifelse(res$log2FoldChange < -logfc_cutoff & res$padj <= adjp_cutoff, "Downregulated", 
                                "Not Significant"))

# Extract significant DEGs
up_genes <- subset(res, Regulation == "Upregulated")
down_genes <- subset(res, Regulation == "Downregulated")

message(paste("Upregulated genes:", nrow(up_genes)))
message(paste("Downregulated genes:", nrow(down_genes)))

################################# Save Results ####################################

# 8. Save DEG results as CSV files
message("Saving results...")
write.csv(as.data.frame(res), file = "Results/DEG_all_Cancer_vs_Normal.csv", row.names = TRUE)
write.csv(as.data.frame(up_genes), file = "Results/DEG_upregulated_Cancer_vs_Normal.csv", row.names = TRUE)
write.csv(as.data.frame(down_genes), file = "Results/DEG_downregulated_Cancer_vs_Normal.csv", row.names = TRUE)

################################# Volcano Plot ####################################

# 9. Create volcano plot
message("Creating volcano plot...")
res_df <- as.data.frame(res)
res_df$Regulation <- factor(res_df$Regulation, 
                            levels = c("Upregulated", "Downregulated", "Not Significant"))

colors <- c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Regulation)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 12) +
  labs(title = "Volcano Plot: Cancer vs Normal",
       subtitle = paste0("Upregulated: ", nrow(up_genes), " | Downregulated: ", nrow(down_genes)),
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("Results/volcano_Cancer_vs_Normal.png", plot = volcano_plot, 
       width = 10, height = 8, dpi = 300)

################################# Heatmap ####################################

# 10. Generate heatmap of top 25 DEGs
message("Creating heatmap...")

# Get top 25 DEGs by adjusted p-value
top25_genes <- rownames(head(res[order(res$padj), ], 25))

# Extract normalized counts for top genes
mat <- norm_counts[top25_genes, ]

# Z-score normalization (across samples for each gene)
mat_z <- t(scale(t(mat)))

# Create annotation for samples
annotation_col <- data.frame(Condition = pheno$condition)
rownames(annotation_col) <- colnames(mat)

# Define colors for annotation
annotation_colors <- list(Condition = c(Normal = "#4DAF4A", Cancer = "#E41A1C"))

# Create heatmap
png("Results/heatmap_top25_Cancer_vs_Normal.png", width = 1000, height = 1200, res = 150)
pheatmap(mat_z, 
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         main = "Top 25 DEGs: Cancer vs Normal",
         show_colnames = FALSE,
         fontsize_row = 8,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

################################# Summary Report ####################################

# 11. Write summary report
message("Writing summary report...")

summary_lines <- c(
  "===== DIFFERENTIAL GENE EXPRESSION ANALYSIS SUMMARY =====",
  "",
  "Dataset: GSE156451 (RNA-seq data)",
  "",
  "1. GENE ID MAPPING & DUPLICATE HANDLING:",
  paste0("   - Data type: RNA-seq with Ensembl gene IDs (not probe-based microarray data)"),
  paste0("   - Annotation package used: org.Hs.eg.db"),
  paste0("   - Multiple Ensembl IDs can map to the same gene symbol due to:"),
  paste0("     * Alternative transcripts (isoforms)"),
  paste0("     * Gene annotation updates/versions"),
  paste0("     * Paralogous genes with similar sequences"),
  paste0("   - Duplicate handling strategy: For genes with multiple Ensembl IDs, retained the ID with highest mean expression"),
  paste0("   - Number of genes with multiple mappings: ", n_duplicates),
  "",
  "2. DIFFERENTIAL EXPRESSION COMPARISON:",
  paste0("   - Contrast performed: Cancer vs Normal"),
  paste0("   - Method: DESeq2 with apeglm log fold change shrinkage"),
  paste0("   - Significance thresholds: |log2FC| > ", logfc_cutoff, " and adjusted p-value < ", adjp_cutoff),
  "",
  "3. DIFFERENTIAL EXPRESSION RESULTS:",
  paste0("   - Total genes analyzed: ", nrow(res)),
  paste0("   - Upregulated genes (Cancer vs Normal): ", nrow(up_genes)),
  paste0("   - Downregulated genes (Cancer vs Normal): ", nrow(down_genes)),
  paste0("   - Total significant DEGs: ", nrow(up_genes) + nrow(down_genes)),
  "",
  "4. OUTPUT FILES GENERATED:",
  "   - DEG_all_Cancer_vs_Normal.csv (all genes with statistics)",
  "   - DEG_upregulated_Cancer_vs_Normal.csv (upregulated genes only)",
  "   - DEG_downregulated_Cancer_vs_Normal.csv (downregulated genes only)",
  "   - volcano_Cancer_vs_Normal.png (volcano plot visualization)",
  "   - heatmap_top25_Cancer_vs_Normal.png (heatmap of top 25 DEGs)",
  "",
  paste0("Analysis completed: ", Sys.time()),
  "=============================================="
)

writeLines(summary_lines, con = "Results/summary_Cancer_vs_Normal.txt")

# Print summary to console
cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\n\n")


