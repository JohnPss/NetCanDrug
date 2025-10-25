# ===================================================================
# NETCANDRUG PROJECT - SCRIPT 02: DIFFERENTIAL EXPRESSION ANALYSIS
# ===================================================================
#
# OBJECTIVE:
# 1. Load the SummarizedExperiment (SE) object prepared in script 01.
# 2. Perform differential expression analysis (DEA) to compare
#    Tumor vs. Normal samples using the edgeR package.
# 3. Filter results to find Differentially Expressed Genes (DEGs)
#    based on project criteria.
# 4. Save the final DEG list in a .csv file for use in the next
#    phase in Python.
#
# ===================================================================

# --- 1. Load libraries ---
library(SummarizedExperiment)
library(edgeR)
library(dplyr)

# --- 2. Define directories and load data ---
## Robust project root detection helper (reuses the same approach as other scripts)
find_project_root <- function(start_path = NULL) {
  if (is.null(start_path)) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- "--file="
    idx <- grep(file_arg, args)
    if (length(idx) > 0) {
      start_path <- normalizePath(sub(file_arg, "", args[idx[1]]))
    } else if (!is.null(sys.frames()[[1]]$ofile)) {
      start_path <- normalizePath(sys.frames()[[1]]$ofile)
    } else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      path <- rstudioapi::getActiveDocumentContext()$path
      if (nzchar(path)) start_path <- normalizePath(path)
    } else {
      start_path <- normalizePath(getwd())
    }
  }

  current <- normalizePath(dirname(start_path))
  root_markers <- c("config.py", ".git", "requirements.txt", "src")
  while (TRUE) {
    files <- list.files(current, all.files = TRUE)
    if (any(root_markers %in% files)) return(normalizePath(current))
    parent <- dirname(current)
    if (parent == current) break
    current <- parent
  }
  return(normalizePath(dirname(start_path)))
}

project_root <- find_project_root()
data_file <- file.path(project_root, "data", "processed", "tcga_brca_prepared_counts.RData")

if (!file.exists(data_file)) stop("Arquivo .RData não encontrado! Rode primeiro o Script 01.")
load(data_file)

cat("TCGA data loaded successfully.\n")
print(tcga_brca_se)


cat("TCGA data loaded successfully.\n")
print(tcga_brca_se)

# --- 3. Preparation for edgeR ---
# Create DGEList object, which is edgeR's main data structure
# colData(tcga_brca_se)$sample_type contains "Primary Tumor" or "Solid Tissue Normal" information
dge <- DGEList(
  counts = assay(tcga_brca_se),
  group = colData(tcga_brca_se)$sample_type,
  genes = rowData(tcga_brca_se)
)

cat("\nDGEList object created.\n")

# Filter genes with low counts
# Keep only genes with at least 10 counts in at least 70% of samples
# This removes noise and increases statistical power
keep <- filterByExpr(dge, group = dge$samples$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat(paste("Genes filtered. Remaining:", nrow(dge), "genes for analysis.\n"))

# Normalization
# Calculate TMM normalization factors to correct sequencing biases
dge <- calcNormFactors(dge, method = "TMM")
cat("Normalization factors calculated (TMM).\n")

# --- 4. Differential Expression Analysis ---
# Create design matrix: specifies the statistical model to be tested
# We're comparing the two groups: Tumor vs. Normal
design <- model.matrix(~ 0 + group, data = dge$samples)
colnames(design) <- levels(dge$samples$group)
colnames(design) <- make.names(colnames(design)) # Ensure valid column names

cat("Design matrix created.\n")
print(head(design))

# Estimate dispersion (variability) of the data
dge <- estimateDisp(dge, design)
cat("Dispersion estimated.\n")

# Fit generalized linear model (GLM) and perform statistical test
fit <- glmQLFit(dge, design)
contrast_matrix <- makeContrasts(Tumor_vs_Normal = Primary.Tumor - Solid.Tissue.Normal, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast_matrix)

cat("Model fitted and statistical test performed.\n")

# Extract results
results_table <- topTags(qlf, n = Inf)$table

cat("Results table generated. First rows:\n")
print(head(results_table))


# --- 5. Filter and Save DEGs ---
# Apply criteria defined in the project
fdr_threshold <- 0.01
logfc_threshold <- 1.5

# Add gene_symbol as column
# rowData(tcga_brca_se) contains 'gene_name' column with symbols
gene_symbols <- rowData(tcga_brca_se)$gene_name
names(gene_symbols) <- rownames(results_table) # Align with rownames

# Add gene_symbol to results table
results_table$gene_symbol <- gene_symbols[rownames(results_table)]

degs <- results_table %>%
  filter(FDR < fdr_threshold, abs(logFC) > logfc_threshold) %>%
  filter(!is.na(gene_symbol)) %>% # Remove genes without symbol
  arrange(desc(abs(logFC)))

cat(paste("\nTotal DEGs found:", nrow(degs), "\n"))
cat(paste("Up-regulated DEGs (overexpressed):", sum(degs$logFC > 0), "\n"))
cat(paste("Down-regulated DEGs (underexpressed):", sum(degs$logFC < 0), "\n"))

# Save including gene_symbol
output_filename <- "data/processed/deg_list.csv"
write.csv(degs, file = output_filename, row.names = FALSE)

cat("========================================================\n")
cat("DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED!\n")
cat(paste("The file 'deg_list.csv' was saved at:", "data/processed/", "\n"))
cat("This file will be the input for network construction in Python.\n")
cat("========================================================\n")
