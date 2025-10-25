# NetCanDrug: Network-Based Cancer Drug Repositioning Pipeline

## Overview

NetCanDrug is a computational pipeline that identifies potential drug candidates for cancer treatment by integrating:
- Differential gene expression analysis
- Protein-protein interaction networks
- Drug-gene interaction databases
- Pathway enrichment analysis

## Project Structure

```
NetCanDrugTest3/
├── config.py                  # Configuration file with all parameters
├── src/                       # Source code directory
│   ├── 01_prepare_tcga_data.R          # Download and prepare TCGA data
│   ├── 02_differential_expression.R    # Identify differentially expressed genes
│   ├── 03_build_tumor_network.py       # Build PPI network from DEGs
│   ├── 04_analyze_tumor_topology.py    # Calculate network topology metrics
│   ├── 05_process_dgidb.py             # Map drugs to network targets
│   ├── 06_create_id_map.R              # Create UniProt-to-gene mapping
│   ├── 07_pathway_enrichment.py        # Identify dysregulated pathways
│   └── 08_calculate_netcandrug_score.py # Calculate final drug scores
└── data/                      # Data directory (not in repository)
    ├── raw/                   # Raw input data
    └── processed/             # Processed output data
```

## Pipeline Steps

### 1. Prepare TCGA Data (R)
Downloads and prepares breast cancer gene expression data from TCGA.

**Input:** None (downloads from TCGA)  
**Output:** `data/processed/tcga_brca_prepared_counts.RData`

### 2. Differential Expression Analysis (R)
Identifies genes that are differentially expressed between tumor and normal samples.

**Input:** `data/processed/tcga_brca_prepared_counts.RData`  
**Output:** `data/processed/deg_list.csv`

### 3. Build Tumor Network (Python)
Constructs a protein-protein interaction network using DEGs and STRING database.

**Input:**
- `data/processed/deg_list.csv`
- `data/raw/string_links.txt`
- `data/raw/9606.protein.aliases.v12.0.txt`

**Output:** `data/processed/ppi_tumor_network.graphml`

### 4. Analyze Tumor Topology (Python)
Calculates network centrality metrics (degree, betweenness, closeness).

**Input:** `data/processed/ppi_tumor_network.graphml`  
**Output:** `data/processed/tumor_network_topology.csv`

### 5. Process DGIdb (Python)
Maps FDA-approved drugs to targets in the tumor network.

**Input:**
- `data/raw/interactions.tsv` (from DGIdb)
- `data/processed/tumor_network_topology.csv`

**Output:** `data/processed/drug_target_mapping.csv`

### 6. Create ID Map (R)
Creates a mapping between UniProt IDs and gene symbols using biomaRt.

**Input:** None (queries Ensembl)  
**Output:** `data/processed/uniprot_to_genesymbol.csv`

### 7. Pathway Enrichment (Python)
Identifies biological pathways enriched with DEGs.

**Input:**
- `data/processed/deg_list.csv`
- `data/processed/reactome_gene_pathway_map.csv`
- `data/processed/uniprot_to_genesymbol.csv`

**Output:** `data/processed/dysregulated_pathways.csv`

### 8. Calculate NetCanDrug Score (Python)
Integrates all previous results to rank drug candidates.

**Input:**
- `data/processed/tumor_network_topology.csv`
- `data/processed/drug_target_mapping.csv`
- `data/processed/dysregulated_pathways.csv`
- `data/processed/deg_list.csv`

**Output:** `data/processed/final_drug_ranking.csv`

## Requirements

### R Packages
```r
install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "edgeR", "biomaRt"))
install.packages("dplyr")
```

### Python Packages
```bash
pip install pandas numpy networkx python-igraph scikit-learn scipy statsmodels tqdm
```

## Required Data Files

Before running the pipeline, download these files:

1. **STRING Database** (for script 03):
   - Download from: https://string-db.org/cgi/download
   - Files needed:
     - `9606.protein.links.v12.0.txt` → save as `data/raw/string_links.txt`
     - `9606.protein.aliases.v12.0.txt` → save as `data/raw/9606.protein.aliases.v12.0.txt`

2. **DGIdb Interactions** (for script 05):
   - Download from: https://www.dgidb.org/downloads
   - File: `interactions.tsv` → save as `data/raw/interactions.tsv`

3. **Reactome Pathway Map** (for script 07):
   - Download from: https://reactome.org/download-data
   - File: Gene-to-pathway mapping → save as `data/processed/reactome_gene_pathway_map.csv`

## Configuration

All key parameters are centralized in `config.py`. Modify this file to:
- Change file paths
- Adjust analysis thresholds
- Modify scoring weights

Key parameters:
- `STRING_SCORE_THRESHOLD`: Minimum STRING confidence score (default: 700)
- `PATHWAY_FDR_THRESHOLD`: FDR threshold for pathway significance (default: 0.05)
- `SCORE_WEIGHTS`: Weights for final drug scoring (topology: 0.4, pathway: 0.3, expression: 0.2, clinical: 0.1)

## Usage

Run scripts in order:

```bash
# R scripts
Rscript src/01_prepare_tcga_data.R
Rscript src/02_differential_expression.R

# Python scripts
python src/03_build_tumor_network.py
python src/04_analyze_tumor_topology.py
python src/05_process_dgidb.py

# R script
Rscript src/06_create_id_map.R

# Python scripts
python src/07_pathway_enrichment.py
python src/08_calculate_netcandrug_score.py
```

## Output

The final output is `data/processed/final_drug_ranking.csv`, which contains:
- Drug names
- NetCanDrug Score (composite score)
- Number of targets
- Individual scores: topology, pathway, expression, clinical

Drugs are ranked by NetCanDrug Score in descending order.

## Citation

If you use this pipeline, please cite the relevant databases:
- TCGA: https://www.cancer.gov/tcga
- STRING: https://string-db.org
- DGIdb: https://www.dgidb.org
- Reactome: https://reactome.org

## License

This project is for research purposes.
