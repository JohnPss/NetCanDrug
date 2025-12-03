"""
NetCanDrug Project Configuration File

This file contains all the key parameters and paths used throughout the pipeline.
Modify these values according to your environment and requirements.
"""

import os

# ===================================================================
# PROJECT DIRECTORY CONFIGURATION
# ===================================================================
# Base directory for the project (automatically set to project root)
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

# Data directories
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
RAW_DATA_DIR = os.path.join(DATA_DIR, "raw")
PROCESSED_DATA_DIR = os.path.join(DATA_DIR, "processed")

# ===================================================================
# INPUT FILE PATHS
# ===================================================================
# DEG list from differential expression analysis
DEG_FILE = os.path.join(PROCESSED_DATA_DIR, "deg_list.csv")

# STRING database files
STRING_LINKS_FILE = os.path.join(RAW_DATA_DIR, "string_links.txt")
STRING_ALIASES_FILE = os.path.join(RAW_DATA_DIR, "9606.protein.aliases.v12.0.txt")

# DGIdb interactions file
DGIDB_FILE = os.path.join(RAW_DATA_DIR, "interactions.tsv")

# Reactome pathway mapping
PATHWAY_MAP_FILE = os.path.join(RAW_DATA_DIR, "reactome_gene_pathway_map.csv")

# ID mapping file
ID_MAP_FILE = os.path.join(PROCESSED_DATA_DIR, "uniprot_to_genesymbol.csv")

# ===================================================================
# OUTPUT FILE PATHS
# ===================================================================
# Tumor network graph
TUMOR_NETWORK_FILE = os.path.join(PROCESSED_DATA_DIR, "ppi_tumor_network.graphml")

# Topology analysis results
TOPOLOGY_FILE = os.path.join(PROCESSED_DATA_DIR, "tumor_network_topology.csv")

# Drug-target mapping
DRUG_TARGET_FILE = os.path.join(PROCESSED_DATA_DIR, "drug_target_mapping.csv")

# Dysregulated pathways
PATHWAYS_FILE = os.path.join(PROCESSED_DATA_DIR, "dysregulated_pathways.csv")

# Final drug ranking
FINAL_RANKING_FILE = os.path.join(PROCESSED_DATA_DIR, "final_drug_ranking.csv")

# ===================================================================
# ANALYSIS PARAMETERS
# ===================================================================

# Script 03: Build Tumor Network
# STRING confidence score threshold (700 = high confidence)
STRING_SCORE_THRESHOLD = 700

# Preferred sources for gene name mapping
PREFERRED_GENE_SOURCES = [
    "Ensembl_HGNC",
    "BioMart_HUGO",
    "Ensembl_EntrezGene"
]

# Chunk size for reading large STRING files
STRING_CHUNK_SIZE = 1000000

# Script 07b: Pathway Enrichment
# FDR threshold for pathway significance
PATHWAY_FDR_THRESHOLD = 0.05

# Script 08: NetCanDrug Score Calculation
# Score weights for final ranking
SCORE_WEIGHTS = {
    'topology': 0.45,
    'pathway': 0.40,
    'expression': 0.05,
    'clinical': 0.00     # AUMENTAR de 0.1
}
# Topology score components weights
TOPOLOGY_DEGREE_WEIGHT = 0.6
TOPOLOGY_BETWEENNESS_WEIGHT = 0.4

# ===================================================================
# DATA FILTERING PARAMETERS
# ===================================================================

# Drug filtering patterns
# Database ID patterns to remove from drug names
DATABASE_ID_PATTERNS = r'IUPHAR\.LIGAND|CHEMBL|CHEBI|DRUGBANK|PUBCHEM'

# Supplements and natural products to exclude
SUPPLEMENTS_TO_EXCLUDE = [
    'ECHINACEA',
    'GINSENG',
    'GREEN TEA',
    'TURMERIC',
    'GARLIC',
    'VITAMIN'
]

# ===================================================================
# INTERACTION TYPE KEYWORDS
# ===================================================================

# Keywords for inhibitor interactions
INHIBITOR_KEYWORDS = ['inhibitor', 'antagonist', 'blocker', 'suppressor']

# Keywords for activator interactions
ACTIVATOR_KEYWORDS = ['activator', 'agonist', 'inducer', 'potentiator']

# ===================================================================
# FUNCTION TO CREATE DIRECTORIES
# ===================================================================

def create_directories():
    """Create necessary directories if they don't exist"""
    directories = [DATA_DIR, RAW_DATA_DIR, PROCESSED_DATA_DIR]
    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)
            print(f"Created directory: {directory}")

# ===================================================================
# VALIDATION FUNCTION
# ===================================================================

def validate_input_files(*files):
    """
    Validate that required input files exist.
    
    Args:
        *files: Variable number of file paths to validate
        
    Returns:
        bool: True if all files exist, False otherwise
    """
    missing_files = []
    for file_path in files:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    if missing_files:
        print("ERROR: The following required files are missing:")
        for file_path in missing_files:
            print(f"  - {file_path}")
        return False
    return True
