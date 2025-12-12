# ===================================================================
# NETCANDRUG PROJECT - SCRIPT 06: CREATE ID TRANSLATION MAP
# ===================================================================
# OBJECTIVE: Use the biomaRt package to create a reliable translation
#            map from UniProt ID to Gene Symbol.
# ===================================================================

# Install biomaRt if not already installed
if (!requireNamespace("biomaRt", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("biomaRt")
}

library(biomaRt)

# Auto-detect project root and set working directory (works with Rscript / RStudio)
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
setwd(project_root)

print("--- Connecting to Ensembl via biomaRt ---")
# Use Ensembl database for human genes
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

print("--- Fetching ID mapping (may take a minute) ---")
# Fetch the attributes we want: UniProt ID and gene name (HGNC symbol)
id_map <- getBM(
    attributes = c("uniprot_gn_id", "hgnc_symbol"),
    mart = ensembl
)

# Clean the result
# Remove rows where one of the IDs is missing
id_map_clean <- id_map[id_map$uniprot_gn_id != "" & id_map$hgnc_symbol != "", ]
# Rename columns for consistency with Python
colnames(id_map_clean) <- c("uniprot_id", "gene_symbol")

print(paste("Mapping created with", nrow(id_map_clean), "entries."))

# Save the mapping file
output_file <- "data/processed/uniprot_to_genesymbol.csv"
write.csv(id_map_clean, file = output_file, row.names = FALSE)

cat("\n========================================================\n")
cat("TRANSLATION MAP CREATED SUCCESSFULLY!\n")
cat(paste("File saved at:", output_file, "\n"))
cat("Now run the Python script 07.\n")
cat("========================================================\n")
