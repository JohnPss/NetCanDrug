# ===================================================================
# NETCANDRUG PROJECT - SCRIPT 01: TCGA DATA PREPARATION PIPELINE
# ===================================================================
#
# OBJECTIVE:
# 1. Define a query to fetch Breast Cancer data (TCGA-BRCA).
# 2. Download the data (TCGAbiolinks will skip existing files).
# 3. Prepare the data in a SummarizedExperiment object.
# 4. Save the final object for future use.
#
# This workflow (Query -> Download -> Prepare) is more robust and avoids
# missing metadata errors.
#
# ===================================================================

# --- 1. Load required libraries ---
library(TCGAbiolinks)
library(SummarizedExperiment)

# --- 2. Define working directory ---
# Get the script directory and set project root
script_dir <- getwd()
work_dir <- dirname(script_dir)
setwd(work_dir)

# --- 3. QUERY Step: Define exactly which data we want ---
cat("STEP 1: Creating GDC query...\n")

query_brca <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

cat("Query created successfully. Summary:\n")
print(query_brca)


# --- 4. DOWNLOAD Step: Download data defined by the query ---
cat("\nSTEP 2: Downloading data (this may take a while)... \n")
# TCGAbiolinks is smart: it will check files in the directory
# and only download what's missing
GDCdownload(
  query = query_brca,
  method = "api",
  directory = "data/raw/GDCdata" # Output directory
)
cat("Download completed.\n")


# --- 5. PREPARE Step: Organize downloaded data ---
cat("\nSTEP 3: Preparing data in a SummarizedExperiment object...\n")
# GDCprepare uses the 'query' which has all metadata,
# avoiding previous errors
tcga_brca_se <- GDCprepare(
  query = query_brca,
  directory = "data/raw/GDCdata"
)

cat("Preparation completed. Object structure:\n")
print(tcga_brca_se)


# --- 6. Save processed object ---
cat("\nSTEP 4: Saving processed object for future use...\n")
output_path <- file.path(work_dir, "data/processed/")

# Create directory if it doesn't exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

output_filename <- file.path(output_path, "tcga_brca_prepared_counts.RData")
save(tcga_brca_se, file = output_filename)


cat("========================================================\n")
cat("TCGA PIPELINE COMPLETED!\n")
cat(paste("The file 'tcga_brca_prepared_counts.RData' was saved at:", output_path, "\n"))
cat("The next script (02_differential_expression.R) will load this file.\n")
cat("========================================================\n")
