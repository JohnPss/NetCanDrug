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

# Get the script directory and set project root
setwd("/home/johnpss/Desktop/NetCanDrugTest1/NetCanDrugTest3")

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
