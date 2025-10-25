# ===================================================================
# NETCANDRUG PROJECT - SCRIPT 03: BUILD TUMOR PPI NETWORK
# ===================================================================
#
# OBJECTIVE:
# 1. Load the list of Differentially Expressed Genes (DEGs) generated in R.
# 2. Map gene names to STRING database IDs.
# 3. Filter the STRING interaction database to keep only high-confidence
#    interactions between DEGs.
# 4. Build the corresponding graph (network) using the NetworkX library.
# 5. Add node attributes (logFC, FDR) for future analyses.
# 6. Save the final graph in GraphML format.
#
# ===================================================================

import pandas as pd
import networkx as nx
import os
import sys

# Add parent directory to path to import config
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

# --- 1. Define paths and parameters ---
print("--- STEP 1: Defining paths and parameters ---")

# Input file paths
deg_file = config.DEG_FILE
string_links_file = config.STRING_LINKS_FILE
string_aliases_file = config.STRING_ALIASES_FILE

# Output file path
output_graph_file = config.TUMOR_NETWORK_FILE

# Parameters
SCORE_THRESHOLD = config.STRING_SCORE_THRESHOLD

# --- 2. Load DEGs and create a set for fast lookup ---
print("\n--- STEP 2: Loading DEG list ---")
degs_df = pd.read_csv(deg_file)

# Use the 'gene_name' column which already contains the correct symbols
# Remove null values in case any gene has no name
deg_genes = set(degs_df['gene_name'].dropna())

print(f"Total of {len(deg_genes)} differentially expressed genes loaded.")

# --- 3. Map gene names to STRING IDs ---
print("\n--- STEP 3: Mapping gene names to STRING IDs ---")
if not os.path.exists(string_aliases_file):
    print("\n!!! ERROR: STRING aliases file not found !!!")
    print(f"Please download '9606.protein.aliases.v12.0.txt.gz' from the STRING website,")
    print(f"uncompress it and place it in '{os.path.dirname(string_aliases_file)}'")
    exit()

# Create a dictionary to map gene_name -> string_id
gene_to_string_id = {}
# Reliable sources for official gene names
preferred_sources = config.PREFERRED_GENE_SOURCES

print("Reading aliases file and searching for reliable sources (this may take a minute)...")
with open(string_aliases_file, 'r') as f:
    next(f) # Skip header
    for line in f:
        try:
            string_id, alias, source = line.strip().split('\t')
            # If the gene name is in our DEG list...
            if alias in deg_genes:
                # ...and if it comes from one of our preferred sources...
                if any(ps in source for ps in preferred_sources):
                    # ...add it to our translation dictionary
                    if alias not in gene_to_string_id:
                        gene_to_string_id[alias] = string_id
        except ValueError:
            # Ignore malformed lines in the aliases file
            continue

# Create a set with the STRING IDs of our DEGs for filtering
deg_string_ids = set(gene_to_string_id.values())

print(f"Mapping complete. {len(deg_string_ids)} of {len(deg_genes)} DEGs were mapped to STRING IDs.")

if len(deg_string_ids) == 0:
    print("\n!!! WARNING: No genes were mapped. Please check if gene names in 'deg_list.csv'")
    print("and the sources in 'preferred_sources' are correct. Exiting script.")
    exit()

# --- 4. Build network from STRING interactions ---
print("\n--- STEP 4: Building network based on interactions ---")
print("Reading STRING interactions file and filtering (this may take a while)...")

# Read interactions file in chunks
chunk_iter = pd.read_csv(
    string_links_file,
    sep=' ',
    usecols=['protein1', 'protein2', 'combined_score'],
    chunksize=config.STRING_CHUNK_SIZE
)

filtered_interactions = []
for chunk in chunk_iter:
    chunk_filtered = chunk[chunk['combined_score'] >= SCORE_THRESHOLD]
    chunk_filtered = chunk_filtered[
        chunk_filtered['protein1'].isin(deg_string_ids) &
        chunk_filtered['protein2'].isin(deg_string_ids)
    ]
    filtered_interactions.append(chunk_filtered)

interactions_df = pd.concat(filtered_interactions, ignore_index=True)

print(f"Filtering complete. {len(interactions_df)} interactions retained.")

# Create the graph
G = nx.from_pandas_edgelist(
    interactions_df,
    source='protein1',
    target='protein2',
    edge_attr='combined_score'
)

print("Initial graph constructed.")

# --- 5. Add attributes to nodes ---
print("\n--- STEP 5: Adding attributes (logFC, FDR) to network nodes ---")

# Create a reverse dictionary to map string_id -> gene_name
string_id_to_gene = {v: k for k, v in gene_to_string_id.items()}

# Use 'gene_name' as index for attribute lookup
degs_df.set_index('gene_name', inplace=True)
node_attributes = {}

for node in G.nodes():
    gene_symbol = string_id_to_gene.get(node)
    if gene_symbol and gene_symbol in degs_df.index:
        node_attributes[node] = {
            'gene_symbol': gene_symbol,
            'logFC': degs_df.loc[gene_symbol, 'logFC'],
            'FDR': degs_df.loc[gene_symbol, 'FDR']
        }

nx.set_node_attributes(G, node_attributes)

print("Attributes added successfully.")

# --- 6. Final report and graph saving ---
print("\n--- STEP 6: Final report and saving ---")
print("\n=== TUMOR NETWORK STATISTICS ===")
print(f"Number of nodes (proteins): {G.number_of_nodes()}")
print(f"Number of edges (interactions): {G.number_of_edges()}")
print("======================================")

# Ensure output directory exists
os.makedirs(os.path.dirname(output_graph_file), exist_ok=True)

nx.write_graphml(G, output_graph_file)

print(f"\nGraph successfully saved at: '{output_graph_file}'")
print("\n========================================================\n")

