# --------------------------
# Assignment 2
# --------------------------
# In this assignment you will work with the results of differential gene expression (DGE) analysis. 
#The analysis produces two key measures for each gene:

# log2FoldChange (log2FC): 
# Indicates the magnitude and direction of change in gene expression. 
# Positive values suggest higher expression(upregulated gene) in the experimental condition compared to control. 
# Negative values suggest lower expression (downregulated gene). 
# The absolute value reflects the strength of the change.

# Adjusted p-value (padj): 
# Represents the statistical significance of the observed difference, corrected for multiple testing. 
# A smaller value indicates stronger evidence that the observed difference is not due to chance.

# Write a function classify_gene() 

# that takes:
#   - logFC (log2FoldChange)
#   - padj  (adjusted p-value)

# and returns:
#   - "Upregulated" if log2FC > 1 and padj < 0.05
#   - "Downregulated" if log2FC < -1 and padj < 0.05
#   - "Not_Significant" otherwise


# Then:
#   - Apply it in a for-loop to process both datasets (DEGs_data_1.csv, DEGs_data_2.csv)
#   - Replace missing padj values with 1
#   - Add a new column 'status'
#   - Save processed files into Results folder
#   - Print summary counts of significant, upregulated, and downregulated genes
#   - Use table() for summaries

# Data Availability
# The input files are available in the GitHub repository:
#      DEGs_Data_1.csv
#      DEGs_Data_2.csv

# Each file contains three columns: 
# Gene_Id	
# padj	
# logFC
############################### @@@@@@@@@@@@@@@ ########################
# Step 1: Define classification function
# --------------------------
classify_gene <- function(logFC, padj) {
  if (is.na(padj)) padj <- 1   # Treat missing padj as 1
  if (logFC > 1 & padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

# --------------------------
# Step 2: Prepare output folder
# --------------------------
if (!dir.exists("Results")) {
  dir.create("Results")
}
### Load the data ##
# --------------------------
# Step 3: List input files
# --------------------------
files <- c("raw_data/DEGs_Data_1.csv", "raw_data/DEGs_Data_2.csv")

# --------------------------
# Step 4: Process each file
# --------------------------
for (file in files) {
  cat("\nProcessing file:", file, "\n")
  
  # Load dataset
  df <- read.csv(file)
  
  # Replace missing padj values with 1
  df$padj[is.na(df$padj)] <- 1
  
  # Classify each gene
  df$status <- mapply(classify_gene, df$logFC, df$padj)
  
  # Save processed data
  out_file <- paste0("Results/Processed_", basename(file))
  write.csv(df, out_file, row.names = FALSE)
  
  # --------------------------
  # Summaries
  # --------------------------
  cat("Overall classification:\n")
  print(table(df$status))
  
  cat("Significant genes (Up or Down):\n")
  print(table(df$status[df$status != "Not_Significant"]))
  
  cat("Upregulated only:\n")
  print(table(df$status == "Upregulated"))
  
  cat("Downregulated only:\n")
  print(table(df$status == "Downregulated"))
}
#################################################### END ######################################################