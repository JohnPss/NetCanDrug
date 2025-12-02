# src/script_08_final.py
"""
Final merged and improved version of SCRIPT 08: Calculate Final NetCanDrug Score
Replace existing SCRIPT 08 with this file after review.
"""

import os
import sys
import logging
from tqdm import tqdm

import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler

# Add parent directory to path to import config
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

# Setup logging based on config
logging.basicConfig(level=logging.DEBUG if getattr(config, "VERBOSE", False) else logging.INFO,
                    format="%(asctime)s %(levelname)s: %(message)s")

def safe_get(d, k, default=0):
    return d.get(k, default)

def calc_target_diversity_penalty(n_targets, shift):
    # sigmoid-based penalty, exposed param shift in config
    return 1.0 / (1.0 + np.exp(-n_targets + shift))

def normalize_df_column(df, col):
    # return normalized column (0-1) using MinMaxScaler
    scaler = MinMaxScaler()
    if df[col].nunique() > 1:
        return scaler.fit_transform(df[[col]]).flatten()
    else:
        return df[col].values

def main():
    logging.info("SCRIPT 08 FINAL: Calculate NetCanDrug Score (merged + improved)")

    # --- 1. Paths & config ---
    topology_file = config.TOPOLOGY_FILE
    drug_target_file = config.DRUG_TARGET_FILE
    pathway_file = config.PATHWAYS_FILE
    deg_file = config.DEG_FILE
    full_pathway_map_file = config.PATHWAY_MAP_FILE
    id_map_file = config.ID_MAP_FILE
    output_file = config.FINAL_RANKING_FILE

    WEIGHTS = config.SCORE_WEIGHTS
    SHIFT = getattr(config, "TARGET_PENALTY_SHIFT", 2.0)
    UNKNOWN_CREDIT = getattr(config, "UNKNOWN_INTERACTION_CREDIT", 0.6)
    OTHER_CREDIT = getattr(config, "OTHER_INTERACTION_CREDIT", 0.4)

    # --- 2. Load files ---
    try:
        df_topology = pd.read_csv(topology_file)
        df_drug_target = pd.read_csv(drug_target_file)
        df_pathways = pd.read_csv(pathway_file)
        df_degs = pd.read_csv(deg_file)
        df_full_pathway_map = pd.read_csv(full_pathway_map_file)
        df_id_map = pd.read_csv(id_map_file)
        logging.info("All files loaded successfully.")
    except FileNotFoundError as e:
        logging.error(f"File not found: {e.filename}")
        sys.exit(1)

    # --- 3. Standardize gene names (very important) ---
    # Ensure gene symbol columns exist
    for df_name, df in [("topology", df_topology), ("degs", df_degs), ("full_pathway_map", df_full_pathway_map)]:
        if 'gene_symbol' in df.columns:
            df['gene_symbol'] = df['gene_symbol'].astype(str).str.upper().str.strip()
    # Specific columns used:
    if 'gene_name' in df_degs.columns:
        df_degs['gene_name'] = df_degs['gene_name'].astype(str).str.upper().str.strip()
    if 'gene_symbol' in df_topology.columns:
        df_topology['gene_symbol'] = df_topology['gene_symbol'].astype(str).str.upper().str.strip()
    if 'uniprot_id' in df_full_pathway_map.columns:
        df_full_pathway_map['uniprot_id'] = df_full_pathway_map['uniprot_id'].astype(str).str.strip()

    # --- 4. Topology scores per gene (normalized separately) ---
    # Use separate scalers (explicit)
    scaler_degree = MinMaxScaler()
    scaler_between = MinMaxScaler()
    if 'degree' not in df_topology.columns or 'betweenness_centrality' not in df_topology.columns:
        logging.error("Expected 'degree' and 'betweenness_centrality' in topology file.")
        sys.exit(1)

    df_topology['degree_norm'] = scaler_degree.fit_transform(df_topology[['degree']])
    df_topology['betweenness_norm'] = scaler_between.fit_transform(df_topology[['betweenness_centrality']])
    df_topology['topology_score_gene'] = (config.TOPOLOGY_DEGREE_WEIGHT * df_topology['degree_norm'] +
                                          config.TOPOLOGY_BETWEENNESS_WEIGHT * df_topology['betweenness_norm'])
    gene_to_topology = df_topology.set_index('gene_symbol')['topology_score_gene'].to_dict()

    # --- 5. Drug-target mapping (list) ---
    if 'drug_name' not in df_drug_target.columns or 'target_gene' not in df_drug_target.columns:
        logging.error("Expected 'drug_name' and 'target_gene' in drug_target file.")
        sys.exit(1)

    # Normalize drug_name and target_gene columns
    df_drug_target['drug_name'] = df_drug_target['drug_name'].astype(str).str.upper().str.strip()
    df_drug_target['target_gene'] = df_drug_target['target_gene'].astype(str).str.upper().str.strip()
    drug_to_targets = df_drug_target.groupby('drug_name')['target_gene'].apply(list).to_dict()

    # --- 6. Pathway mapping ---
    uniprot_to_gene = pd.Series(df_id_map.gene_symbol.values, index=df_id_map.uniprot_id).to_dict()
    df_full_pathway_map['gene_symbol'] = df_full_pathway_map['uniprot_id'].map(uniprot_to_gene)
    df_full_pathway_map.dropna(subset=['gene_symbol'], inplace=True)
    df_full_pathway_map['gene_symbol'] = df_full_pathway_map['gene_symbol'].astype(str).str.upper().str.strip()

    dysregulated_pathways_set = set(df_pathways['pathway_name'].unique())
    pathway_map_df_dys = df_full_pathway_map[df_full_pathway_map['pathway_name'].isin(dysregulated_pathways_set)]
    gene_to_dys_pathways = pathway_map_df_dys.groupby('gene_symbol')['pathway_name'].apply(set).to_dict()
    total_dysregulated_pathways = len(dysregulated_pathways_set)

    # --- 7. DEGs map (logFC) ---
    gene_to_logfc = df_degs.set_index('gene_name')['logFC'].to_dict()

    # --- 8. Iterate drugs and compute scores ---
    results = []
    unique_drugs = list(df_drug_target['drug_name'].unique())
    logging.info(f"Calculating scores for {len(unique_drugs)} unique drugs...")

    for drug in tqdm(unique_drugs, desc="Calculating scores"):
        targets = drug_to_targets.get(drug, [])
        if not targets:
            continue

        # Topology score
        ts_raw = np.mean([safe_get(gene_to_topology, t, 0.0) for t in targets]) if targets else 0.0
        penalty = calc_target_diversity_penalty(len(set(targets)), SHIFT)
        ts_raw = ts_raw * penalty

        # Pathway score
        pathways_hit = set().union(*[gene_to_dys_pathways.get(t, set()) for t in targets]) if targets else set()
        ps_raw = len(pathways_hit) / total_dysregulated_pathways if total_dysregulated_pathways > 0 else 0.0

        # Expression score (improved)
        expression_scores = []
        interactions_for_drug = df_drug_target[df_drug_target['drug_name'] == drug]
        for _, row in interactions_for_drug.iterrows():
            target = str(row['target_gene']).upper().strip()
            interaction_raw = row.get('interaction_type', '')
            interaction = str(interaction_raw).lower() if not pd.isna(interaction_raw) else ''
            logfc = gene_to_logfc.get(target)
            score = 0.0
            if logfc is not None:
                # Determine inhibitor/activator
                is_inhibitor = any(k in interaction for k in getattr(config, 'INHIBITOR_KEYWORDS', []))
                is_activator = any(k in interaction for k in getattr(config, 'ACTIVATOR_KEYWORDS', []))

                if is_inhibitor:
                    if logfc > 0:
                        score = abs(logfc)
                elif is_activator:
                    if logfc < 0:
                        score = abs(logfc)
                else:
                    # Unknown or other interaction types: give partial credit
                    if interaction in ['', 'nan', 'n/a', 'none', 'unknown', 'other']:
                        score = abs(logfc) * UNKNOWN_CREDIT
                    else:
                        score = abs(logfc) * OTHER_CREDIT

            expression_scores.append(score)

        es_raw = float(np.mean(expression_scores)) if expression_scores else 0.0

        # Clinical score placeholder
        cs_raw = 0.0

        results.append({
            'drug_name': drug,
            'topology_score_raw': ts_raw,
            'pathway_score_raw': ps_raw,
            'expression_score_raw': es_raw,
            'clinical_score_raw': cs_raw,
            'num_targets': len(set(targets))
        })

    df_scores = pd.DataFrame(results)
    logging.info(f"Raw scores calculated for {len(df_scores)} drugs.")

    # Debug info (only if verbose)
    if getattr(config, "VERBOSE", False):
        logging.debug("Raw scores (summary):")
        logging.debug(df_scores[['expression_score_raw','topology_score_raw','pathway_score_raw']].describe())

    # --- 9. Normalize and final score ---
    for score_col in ['topology_score_raw', 'pathway_score_raw', 'expression_score_raw']:
        norm_col = score_col.replace('_raw', '_norm')
        df_scores[norm_col] = normalize_df_column(df_scores, score_col)

    df_scores['cs_norm'] = df_scores['clinical_score_raw']  # already normalized placeholder

    df_scores['NetCanDrug_Score'] = (
        WEIGHTS['topology'] * df_scores['topology_score_norm'] +
        WEIGHTS['pathway'] * df_scores['pathway_score_norm'] +
        WEIGHTS['expression'] * df_scores['expression_score_norm'] +
        WEIGHTS['clinical'] * df_scores['cs_norm']
    )

    final_ranking = df_scores.sort_values('NetCanDrug_Score', ascending=False)
    final_ranking = final_ranking[['drug_name', 'NetCanDrug_Score', 'num_targets',
                                   'topology_score_norm', 'pathway_score_norm', 'expression_score_norm', 'cs_norm']]

    # --- 10. Save ---
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    final_ranking.to_csv(output_file, index=False)
    logging.info(f"Final ranking saved at: {output_file}")

    # Present top20
    print("\n🔝 Top 20 candidate drugs (NetCanDrug):")
    print(final_ranking.head(20).to_string(index=False,
                                          formatters={'NetCanDrug_Score':'{:.4f}'.format,
                                                      'topology_score_norm':'{:.3f}'.format,
                                                      'pathway_score_norm':'{:.3f}'.format,
                                                      'expression_score_norm':'{:.3f}'.format,
                                                      'cs_norm':'{:.1f}'.format}))

if __name__ == "__main__":
    main()
