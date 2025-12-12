# ===================================================================
# NETCANDRUG PROJECT - SCRIPT 04: TUMOR NETWORK TOPOLOGY ANALYSIS
# ===================================================================
#
# METHODOLOGICAL APPROACH:
# Differential (Δ) analysis between networks of different sizes is methodologically
# flawed. The correct and more robust approach is to analyze the topology
# of the disease (tumor) network directly. The hubs of this network are the targets.
#
# OBJECTIVE OF THIS SCRIPT:
# 1. Load ONLY the tumor network.
# 2. Calculate centrality metrics (Degree, Betweenness, Closeness).
# 3. Save this table, which will be the input for the Topology Score.
#
# ===================================================================

import pandas as pd
import networkx as nx
import igraph as ig
from tqdm import tqdm
import os
import sys

# Add project root to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

try:
    from src.config import config
except ImportError:
    # Fallback for when running from root
    sys.path.append(os.getcwd())
    from src.config import config

print("=" * 80)
print("SCRIPT 04: TUMOR NETWORK TOPOLOGY ANALYSIS")
print("=" * 80)

# --- Configuration ---
tumor_network_file = config.TUMOR_NETWORK_FILE
output_file = config.TOPOLOGY_FILE

# --- STEP 1: Load Tumor Network ---
print("\n--- STEP 1: Loading tumor network ---")
if not os.path.exists(tumor_network_file):
    print(f"❌ ERROR: Tumor network file not found at '{tumor_network_file}'")
    print("   Please run Script 03 first.")
    exit()

G_tumor_nx = nx.read_graphml(tumor_network_file)
print(f"  ✅ Tumor network loaded: {G_tumor_nx.number_of_nodes()} nodes, {G_tumor_nx.number_of_edges()} edges")

# --- STEP 2: Calculate Centrality Metrics (Fast Version with iGraph) ---
print("\n--- STEP 2: Calculating centrality metrics (may take 1-2 minutes) ---")

# Convert to iGraph for performance
try:
    G_tumor_ig = ig.Graph.from_networkx(G_tumor_nx)
    node_names = G_tumor_ig.vs["_nx_name"] # In newer igraph versions
except KeyError:
    node_names = G_tumor_ig.vs["id"] # In older versions

print("  - Calculating Degree...")
degree = G_tumor_ig.degree()

print("  - Calculating Betweenness...")
betweenness = G_tumor_ig.betweenness(directed=False)

print("  - Calculating Closeness...")
closeness = G_tumor_ig.closeness()

# Normalize betweenness to [0, 1] range for consistency
n = len(node_names)
if n > 2:
    # Standard networkx normalization is 1/((n-1)*(n-2)), but igraph doesn't normalize by default.
    # We'll apply simple normalization by dividing by maximum, which is sufficient for ranking.
    max_betweenness = max(betweenness)
    if max_betweenness > 0:
        betweenness = [b / max_betweenness for b in betweenness]

print("  ✅ Metrics calculated successfully.")

# --- STEP 3: Assemble and Save Results Table ---
print("\n--- STEP 3: Assembling results table ---")

# Map STRING ID to Gene Symbol
gene_map = nx.get_node_attributes(G_tumor_nx, 'gene_symbol')

# Create DataFrame
df_topology = pd.DataFrame({
    'string_id': node_names,
    'degree': degree,
    'betweenness_centrality': betweenness,
    'closeness_centrality': closeness
})

# Add gene name (gene_symbol)
df_topology['gene_symbol'] = df_topology['string_id'].map(gene_map)

# Reorder columns for better visualization
df_topology = df_topology[['gene_symbol', 'string_id', 'degree', 'betweenness_centrality', 'closeness_centrality']]

# Sort by most important proteins (highest degree, then highest betweenness)
df_topology = df_topology.sort_values(by=['degree', 'betweenness_centrality'], ascending=False)

# Ensure output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Save to file
df_topology.to_csv(output_file, index=False)

print(f"  ✅ Topology table saved at: {output_file}")

# --- STEP 4: Final Report ---
print("\n" + "=" * 80)
print("TUMOR NETWORK TOPOLOGY REPORT")
print("=" * 80)

print(f"\n📈 DEGREE STATISTICS (Connectivity):")
print(f"  - Mean:    {df_topology['degree'].mean():.2f}")
print(f"  - Maximum: {df_topology['degree'].max()} (Most connected protein)")
print(f"  - Minimum: {df_topology['degree'].min()}")

print(f"\n🔝 TOP 10 MOST IMPORTANT PROTEINS (HUBS) IN THE CANCER NETWORK:")
print(df_topology[['gene_symbol', 'degree', 'betweenness_centrality']].head(10).to_string(index=False))

print("\n" + "=" * 80)
print("SCRIPT 04 COMPLETED SUCCESSFULLY!")
print("=" * 80)
print("\n🎯 NEXT STEPS:")
print("  - The file 'tumor_network_topology.csv' is the final basis for the Topology Score.")
print("  - Now we can proceed to Script 05: Process DGIdb.")
print("=" * 80)

