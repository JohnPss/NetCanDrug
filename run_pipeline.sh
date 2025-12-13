#!/bin/bash

# ===================================================================
# NETCANDRUG PIPELINE RUNNER
# ===================================================================
# This script executes the pipeline steps 01 through 08 sequentially.
# It assumes that 'python3' and 'Rscript' are available in the PATH.

# Stop execution if any command fails
set -e

echo "==================================================================="
echo "STARTING NETCANDRUG PIPELINE"
echo "==================================================================="

# --- Step 01: Prepare TCGA Data ---
echo ""
echo ">>> Running Step 01: Prepare TCGA Data..."
Rscript src/pipeline/01_prepare_tcga_data.R

# --- Step 02: Differential Expression Analysis ---
echo ""
echo ">>> Running Step 02: Differential Expression Analysis..."
Rscript src/pipeline/02_differential_expression.R

# --- Step 03: Build Tumor PPI Network ---
echo ""
echo ">>> Running Step 03: Build Tumor PPI Network..."
python3 src/pipeline/03_build_tumor_network.py

# --- Step 04: Analyze Tumor Topology ---
echo ""
echo ">>> Running Step 04: Analyze Tumor Topology..."
python3 src/pipeline/04_analyze_tumor_topology.py

# --- Step 05: Process DGIdb Data ---
echo ""
echo ">>> Running Step 05: Process DGIdb Data..."
python3 src/pipeline/05_process_dgidb.py

# --- Step 06: Create ID Map ---
echo ""
echo ">>> Running Step 06: Create ID Map..."
Rscript src/pipeline/06_create_id_map.R

# --- Step 07: Pathway Enrichment ---
echo ""
echo ">>> Running Step 07: Pathway Enrichment..."
python3 src/pipeline/07_pathway_enrichment.py

# --- Step 08: Calculate NetCanDrug Score ---
echo ""
echo ">>> Running Step 08: Calculate NetCanDrug Score..."
python3 src/pipeline/08_calculate_netcandrug_score.py

echo ""
echo "==================================================================="
echo "PIPELINE COMPLETED SUCCESSFULLY!"
echo "==================================================================="
