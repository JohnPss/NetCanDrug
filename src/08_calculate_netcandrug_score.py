# ==================================================================
# NETCANDRUG PROJECT - SCRIPT 08: CALCULATE FINAL NETCANDRUG SCORE
# ==================================================================
#
# OBJECTIVE:
# 1. Integrate all processed results.
# 2. Calculate the four component scores for each drug.
# 3. Generate the final ranking of candidate drugs.
#
# SCORING APPROACH:
# - Expression Score is GRADUATED by the absolute value of logFC.
# - Topology Score includes a soft PENALTY for drugs with very few targets.
#
# ==================================================================

import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm
import os
import sys

# Add parent directory to path to import config
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

print("=" * 80)
print("SCRIPT 08: CALCULATE FINAL NETCANDRUG SCORE")
print("=" * 80)

# --- 1. Path Configuration ---
topology_file = config.TOPOLOGY_FILE
drug_target_file = config.DRUG_TARGET_FILE
pathway_file = config.PATHWAYS_FILE
deg_file = config.DEG_FILE
full_pathway_map_file = config.PATHWAY_MAP_FILE
id_map_file = config.ID_MAP_FILE
output_file = config.FINAL_RANKING_FILE

WEIGHTS = config.SCORE_WEIGHTS

# --- 2. Load all data files ---
print("\n--- STEP 1: Loading all data files ---")
try:
    df_topology = pd.read_csv(topology_file)
    df_drug_target = pd.read_csv(drug_target_file)
    df_pathways = pd.read_csv(pathway_file)
    df_degs = pd.read_csv(deg_file)
    df_full_pathway_map = pd.read_csv(full_pathway_map_file)
    df_id_map = pd.read_csv(id_map_file)
    print("  ✅ All files loaded successfully.")
except FileNotFoundError as e:
    print(f"❌ ERROR: File not found: {e.filename}")
    exit()

# --- 3. Prepare Data Structures ---
print("\n--- STEP 2: Preparing data structures for calculations ---")

scaler = MinMaxScaler()
df_topology['degree_norm'] = scaler.fit_transform(df_topology[['degree']])
df_topology['betweenness_norm'] = scaler.fit_transform(df_topology[['betweenness_centrality']])
df_topology['topology_score_gene'] = (config.TOPOLOGY_DEGREE_WEIGHT * df_topology['degree_norm'] + 
                                       config.TOPOLOGY_BETWEENNESS_WEIGHT * df_topology['betweenness_norm'])
gene_to_topology = df_topology.set_index('gene_symbol')['topology_score_gene'].to_dict()

drug_to_targets = df_drug_target.groupby('drug_name')['target_gene'].apply(list).to_dict()

uniprot_to_gene = pd.Series(df_id_map.gene_symbol.values, index=df_id_map.uniprot_id).to_dict()
df_full_pathway_map['gene_symbol'] = df_full_pathway_map['uniprot_id'].map(uniprot_to_gene)
df_full_pathway_map.dropna(subset=['gene_symbol'], inplace=True)

dysregulated_pathways_set = set(df_pathways['pathway_name'])
pathway_map_df_dys = df_full_pathway_map[df_full_pathway_map['pathway_name'].isin(dysregulated_pathways_set)]
gene_to_dys_pathways = pathway_map_df_dys.groupby('gene_symbol')['pathway_name'].apply(set).to_dict()
total_dysregulated_pathways = len(dysregulated_pathways_set)

gene_to_logfc = df_degs.set_index('gene_name')['logFC'].to_dict()
print("  ✅ Auxiliary data structures created.")

# --- 4. Calculate Scores ---
print("\n--- STEP 3: Calculating the 4 scores for each drug ---")
results = []
unique_drugs = df_drug_target['drug_name'].unique()

for drug in tqdm(unique_drugs, desc="Calculating scores"):
    targets = drug_to_targets.get(drug, [])
    if not targets: continue

    # --- Calculate Topology Score (TS) ---
    ts = np.mean([gene_to_topology.get(t, 0) for t in targets]) if targets else 0
    
    # PENALTY FOR FEW TARGETS
    # Penalize drugs with only 1 target (may be too specific/fragile)
    # The sigmoid function maps target count to penalty between ~0.27 (1 target) and ~1.0 (many targets)
    target_diversity_penalty = 1 / (1 + np.exp(-len(set(targets)) + 2))
    ts = ts * target_diversity_penalty
    
    # --- Calculate Pathway Score (PS) ---
    pathways_hit = set().union(*[gene_to_dys_pathways.get(t, set()) for t in targets])
    ps = len(pathways_hit) / total_dysregulated_pathways if total_dysregulated_pathways > 0 else 0

    # --- Calculate Expression Score (ES) ---
    expression_scores = []
    drug_interactions = df_drug_target[df_drug_target['drug_name'] == drug]
    for _, row in drug_interactions.iterrows():
        target, interaction = row['target_gene'], str(row['interaction_type']).lower()
        logfc = gene_to_logfc.get(target)
        score = 0
        if logfc is not None:
            is_inhibitor = any(i in interaction for i in config.INHIBITOR_KEYWORDS)
            is_activator = any(a in interaction for a in config.ACTIVATOR_KEYWORDS)
            
            # EXPRESSION SCORE WITH GRADATION
            if is_inhibitor and logfc > 0:
                score = abs(logfc)  # Higher logFC, better score
            elif is_activator and logfc < 0:
                score = abs(logfc)
            # If action is not consistent, score remains 0
            
        expression_scores.append(score)
    es = np.mean(expression_scores) if expression_scores else 0
    cs = 1.0

    results.append({
        'drug_name': drug, 'topology_score_raw': ts, 'pathway_score_raw': ps,
        'expression_score_raw': es, 'clinical_score_raw': cs, 'num_targets': len(set(targets))
    })

df_scores = pd.DataFrame(results)
print("  ✅ Raw scores calculated.")

# --- 5. Normalize and Calculate Final Score ---
print("\n--- STEP 4: Normalizing and calculating final score ---")
for score_col in ['topology_score_raw', 'pathway_score_raw', 'expression_score_raw']:
    norm_col = score_col.replace('_raw', '_norm')
    if df_scores[score_col].nunique() > 1:
        df_scores[norm_col] = MinMaxScaler().fit_transform(df_scores[[score_col]])
    else:
        df_scores[norm_col] = df_scores[score_col]
df_scores['cs_norm'] = df_scores['clinical_score_raw']

df_scores['NetCanDrug_Score'] = (
    WEIGHTS['topology'] * df_scores['topology_score_norm'] +
    WEIGHTS['pathway'] * df_scores['pathway_score_norm'] +
    WEIGHTS['expression'] * df_scores['expression_score_norm'] +
    WEIGHTS['clinical'] * df_scores['cs_norm']
)

final_ranking = df_scores.sort_values('NetCanDrug_Score', ascending=False)
final_ranking = final_ranking[['drug_name', 'NetCanDrug_Score', 'num_targets', 'topology_score_norm', 'pathway_score_norm', 'expression_score_norm', 'cs_norm']]
print("  ✅ Final ranking generated.")

# --- 6. Save and Present Results ---
print("\n--- STEP 5: Saving and presenting final ranking ---")

# Ensure output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

final_ranking.to_csv(output_file, index=False)
print(f"  ✅ Final ranking saved at: {output_file}")

print("\n" + "=" * 80)
print("FINAL CANDIDATE RANKING - NETCANDRUG PROJECT")
print("=" * 80)
print("\n🔝 Top 20 Candidate Drugs:")
display_ranking = final_ranking.rename(columns={
    'topology_score_norm': 'ts_norm',
    'pathway_score_norm': 'ps_norm',
    'expression_score_norm': 'es_norm'
})
print(display_ranking.head(20).to_string(index=False, formatters={'NetCanDrug_Score': '{:.4f}'.format, 'ts_norm': '{:.3f}'.format, 'ps_norm': '{:.3f}'.format, 'es_norm': '{:.3f}'.format, 'cs_norm': '{:.1f}'.format}))

print("\n" + "=" * 80)
print("PROJECT COMPLETED SUCCESSFULLY!")
print("=" * 80)

