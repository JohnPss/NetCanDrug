#!/usr/bin/env python3


import os
import sys
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Basic scientific stack
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import hypergeom
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, auc, precision_recall_curve
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.utils import resample

# Reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Make sure results directory exists
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True, parents=True)

# Add project root to path and import config
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

try:
    from src.config import config
except Exception as e:
    print("Warning: could not import config module. Make sure src/config/config.py exists.")
    raise

# ---------- UTILITIES & CONSTANTS ----------

def normalize_drug_names(series: pd.Series) -> pd.Series:
    return (series.astype(str).str.upper().str.strip()
            .str.replace(r'[^A-Z0-9]', '', regex=True))
            

def load_fda_breast_drugs():
    """Load FDA breast cancer drugs from CSV file"""
    fda_file_path = Path(config.DATA_DIR) / "reference" / "fda_breast_cancer_drugs.csv"

    if not fda_file_path.exists():
        raise FileNotFoundError(
            f"FDA drug list not found at: {fda_file_path}"
        )

    fda_df = pd.read_csv(fda_file_path)

    if "drug_name" not in fda_df.columns:
        raise ValueError("CSV must contain column 'drug_name'")

    drug_names = normalize_drug_names(fda_df["drug_name"]).tolist()
    return drug_names, set(drug_names)


FDA_BREAST_DRUGS_LIST, FDA_BREAST_DRUGS = load_fda_breast_drugs()



def normalize_drug_names(series: pd.Series) -> pd.Series:
    return (series.astype(str).str.upper().str.strip()
            .str.replace(r'[^A-Z0-9]', '', regex=True))


def ensure_df_column(df: pd.DataFrame, col: str, default=np.nan):
    if col not in df.columns:
        df[col] = default
    return df


# ---------- SECTION 1: LOAD DATA & PREPARE ----------

# ---------- SECTION 1: LOAD DATA & PREPARE ----------

def recalculate_score_without_clinical(df: pd.DataFrame):
    """
    Recalculates the final score using ONLY biological components
    (Topology, Pathway, Expression) and EXCLUDING Clinical Score/FDA status
    to prevent data leakage during validation.
    """
    print("[PREP] Recalculating scores explicitly excluding Clinical Component...")
    
    # Ensure columns exist and fill NaNs
    cols = ['topology_score_norm', 'pathway_score_norm', 'expression_score_norm']
    for c in cols:
        ensure_df_column(df, c, 0.0)
        
    # Validation Weights (Must sum to 1.0 or be consistent with Ablation 'Full Model')
    # Ajuste estes pesos conforme sua configuração 'Full Model' na Ablation
    W_TOPO = 0.45
    W_PATH = 0.40
    W_EXPR = 0.05
    
    df['Validation_Score'] = (
        W_TOPO * df['topology_score_norm'] +
        W_PATH * df['pathway_score_norm'] +
        W_EXPR * df['expression_score_norm']
    )
    return df

def load_data():
    print("[LOAD] Reading input files...")
    df_ranking = pd.read_csv(config.FINAL_RANKING_FILE)
    df_drug_targets = pd.read_csv(config.DRUG_TARGET_FILE)
    df_deg = pd.read_csv(config.DEG_FILE)
    
    # Normalize drug names
    for df in [df_ranking, df_drug_targets]:
        if 'drug_name' in df.columns:
            df['drug_name'] = normalize_drug_names(df['drug_name'])
            
    # CRITICAL: Create the Clean Validation Score
    df_ranking = recalculate_score_without_clinical(df_ranking)
    
    return df_ranking, df_drug_targets, df_deg


# ---------- SECTION 2: VALIDATION (TRAIN/TEST, Precision@K) ----------

def train_test_validation(df_ranking: pd.DataFrame, fda_list=FDA_BREAST_DRUGS_LIST):
    print("\n[VALIDATION] Train/test split & Precision@K on HELD-OUT FDA drugs")
    
    # create train/test split
    train_drugs, test_drugs = train_test_split(fda_list, test_size=0.4, random_state=RANDOM_SEED)
    
    # Normalização para garantir match
    test_drugs_norm = [d.upper().replace(r'[^A-Z0-9]', '') for d in test_drugs]
    
    df = df_ranking.copy()
    # Usar Validation_Score se existir (criado na etapa anterior), senão usa NetCanDrug_Score
    score_col = 'Validation_Score' if 'Validation_Score' in df.columns else 'NetCanDrug_Score'
    df = df.sort_values(score_col, ascending=False)
    
    df['is_hit'] = df['drug_name'].isin(test_drugs_norm)
    
    test_results = []
    for k in [10, 20, 30, 50, 100]:
        top_k = df.head(k)
        hits = int(top_k['is_hit'].sum())
        precision = hits / k
        
        test_results.append({
            'Metric_Type': 'Held-Out (Test Set)', # <--- NOVA COLUNA
            'K': k, 
            'Precision': precision, 
            'Hits': hits,
            'Total_Positives_In_Set': len(test_drugs)
        })
        
    test_results_df = pd.DataFrame(test_results)
    test_results_df.to_csv(RESULTS_DIR / "validation_heldout.csv", index=False)
    print(test_results_df.to_string(index=False))
    
    return train_drugs, test_drugs, test_results_df


# ---------- SECTION 3: ABLATION STUDY ----------

def ablation_study(df_ranking: pd.DataFrame):
    print("\n[ABLATION] Running ablation study...")
    df = df_ranking.copy()
    # Ensure normalized component columns exist
    ensure_df_column(df, 'topology_score_norm', 0.0)
    ensure_df_column(df, 'pathway_score_norm', 0.0)
    ensure_df_column(df, 'expression_score_norm', 0.0)
    configurations = {
        'Full Model': {'w_topo': 0.45, 'w_path': 0.40, 'w_expr': 0.05},
        'No Topology': {'w_topo': 0.00, 'w_path': 0.60, 'w_expr': 0.40},
        'No Pathway': {'w_topo': 0.65, 'w_path': 0.00, 'w_expr': 0.35},
        'No Expression': {'w_topo': 0.55, 'w_path': 0.45, 'w_expr': 0.00},
        'Topology Only': {'w_topo': 1.00, 'w_path': 0.00, 'w_expr': 0.00},
        'Pathway Only': {'w_topo': 0.00, 'w_path': 1.00, 'w_expr': 0.00},
        'Expression Only': {'w_topo': 0.00, 'w_path': 0.00, 'w_expr': 1.00}
    }
    ablation_results = []
    for name, weights in configurations.items():
        df['score_ablation'] = (
            weights['w_topo'] * df['topology_score_norm'] +
            weights['w_path'] * df['pathway_score_norm'] +
            weights['w_expr'] * df['expression_score_norm']
        )
        df_sorted = df.sort_values('score_ablation', ascending=False)
        precision_20 = sum(df_sorted.head(20)['is_test_fda'].fillna(False)) / 20
        ablation_results.append({'Configuration': name, 'Precision@20': precision_20})
    ablation_df = pd.DataFrame(ablation_results).sort_values('Precision@20', ascending=False)
    ablation_df.to_csv(RESULTS_DIR / "ablation_study.csv", index=False)
    print(ablation_df.to_string(index=False))
    return ablation_df


# ---------- SECTION 4: BASELINE (EXPRESSION ONLY) ----------

def compute_baseline_expression(df_drug_targets: pd.DataFrame, df_deg: pd.DataFrame, test_drugs=None):
    print("\n[BASELINE] Computing expression-only baseline ranking...")
    drug_targets = df_drug_targets.copy()
    deg = df_deg.copy()
    # normalize columns
    drug_targets['drug_name'] = normalize_drug_names(drug_targets['drug_name'])
    deg['gene_name'] = deg['gene_name'].astype(str).str.upper()
    # merge on target gene <-> gene_name
    merged = drug_targets.merge(deg, left_on='target_gene', right_on='gene_name', how='inner')
    merged = merged.dropna(subset=['logFC'])
    def calculate_expression_score(row):
        logfc = row['logFC']
        interaction = str(row.get('interaction_type', '')).lower()
        is_inhibitor = any(k in interaction for k in ['inhibitor', 'antagonist', 'blocker'])
        is_activator = any(k in interaction for k in ['activator', 'agonist', 'inducer'])
        if is_inhibitor and logfc > 0:
            return abs(logfc)
        elif is_activator and logfc < 0:
            return abs(logfc)
        else:
            return abs(logfc) * 0.5
    merged['expr_score'] = merged.apply(calculate_expression_score, axis=1)
    baseline_expr = merged.groupby('drug_name')['expr_score'].mean().reset_index().rename(columns={'expr_score':'baseline_expression_score'})
    baseline_expr = baseline_expr.sort_values('baseline_expression_score', ascending=False)
    baseline_expr.to_csv(RESULTS_DIR / "baseline_expression_ranking.csv", index=False)
    # Precision@20 if test_drugs provided
    baseline_precision_20 = None
    if test_drugs is not None:
        baseline_expr['is_test_fda'] = baseline_expr['drug_name'].isin(test_drugs)
        baseline_precision_20 = sum(baseline_expr.head(20)['is_test_fda']) / 20
        print(f"Baseline Precision@20 (test set): {baseline_precision_20:.3f}")
    print("Saved baseline_expression_ranking.csv")
    return baseline_expr, baseline_precision_20


# ---------- SECTION 5 & 11: UNIFIED PERMUTATION TEST ----------

def run_permutation_test(df_ranking: pd.DataFrame, 
                         target_drug_list: list, 
                         mode_label: str,
                         n_iter=10000, 
                         top_k=20):
    """
    Generic permutation test.
    mode_label: 'Held-Out' or 'Full_FDA' to label the output.
    target_drug_list: The specific list of drugs to count as hits.
    """
    print(f"\n[PERMUTATION] Running test for: {mode_label} (Top {top_k}, {n_iter} iter)...")
    
    # Preparar dados do modelo
    df = df_ranking.copy()
    score_col = 'Validation_Score' if 'Validation_Score' in df.columns else 'NetCanDrug_Score'
    df = df.sort_values(score_col, ascending=False)
    
    # Normalizar lista de alvo
    target_set = set([d.upper().replace(r'[^A-Z0-9]', '') for d in target_drug_list])
    all_drugs = df['drug_name'].dropna().unique()
    
    # 1. Calcular precisão real do modelo
    top_k_drugs = df.head(top_k)['drug_name'].tolist()
    real_hits = len(set(top_k_drugs).intersection(target_set))
    model_precision = real_hits / top_k
    
    # 2. Definir seed local para garantir que rodar 2x dê o mesmo resultado
    rng = np.random.RandomState(RANDOM_SEED)
    
    # 3. Permutações
    random_precisions = []
    for i in range(n_iter):
        # Use 'rng.choice' ao invés de 'np.random.choice'
        shuffled = rng.choice(all_drugs, top_k, replace=False)
        hits = len(set(shuffled).intersection(target_set))
        random_precisions.append(hits / top_k)
        
    random_precisions = np.array(random_precisions)
    mean_random = random_precisions.mean()
    std_random = random_precisions.std()
    
    # CORREÇÃO CRÍTICA DO P-VALOR:
    # Conta quantas vezes o aleatório foi IGUAL ou MELHOR que o modelo
    # Adiciona +1 no numerador e denominador para evitar p=0 (correção de Laplace/conservadora)
    n_better_or_equal = (random_precisions >= model_precision).sum()
    p_val = (n_better_or_equal + 1) / (n_iter + 1)
    
    results = {
        'Metric_Type': mode_label,
        'Top_K': top_k,
        'Model_Precision': model_precision,
        'Random_Mean': mean_random,
        'Random_Std': std_random,
        'Fold_Improvement': model_precision / mean_random if mean_random > 0 else np.inf,
        'Empirical_p_value': p_val
    }
    
    # Print rápido
    print(f"   Model: {model_precision:.3f} | Random: {mean_random:.3f} | p-val: {p_val:.4f}")
    
    return pd.DataFrame([results])

# ---------- SECTION 6: FINAL REPORT & COMPARISONS ----------

def generate_final_comparison(netcandrug_precision20, baseline_precision20, mean_random):
    comparison_df = pd.DataFrame({
        'Method': ['NetCanDrug (Multi-component)', 'Baseline (Expression only)', 'Random (Mean of 1000)'],
        'Precision@20': [netcandrug_precision20, baseline_precision20, mean_random]
    })
    def fold(x): return f"{x/mean_random:.2f}x" if mean_random>0 else "inf"
    comparison_df['Fold vs Random'] = [fold(x) for x in comparison_df['Precision@20']]
    comparison_df.to_csv(RESULTS_DIR / "baseline_comparisons.csv", index=False)
    print("\nMethod comparison:")
    print(comparison_df.to_string(index=False))
    return comparison_df


# ---------- SECTION 7: Siltuximab CASE STUDY (CORRIGIDO) ----------

def siltuximab_case_study(ppi_graph_path=config.TUMOR_NETWORK_FILE,
                         drug_targets_path=config.DRUG_TARGET_FILE,
                         pathway_map_path=config.PATHWAY_MAP_FILE,
                         id_map_path=config.ID_MAP_FILE):
    print("\n[CASE STUDY] Siltuximab (IL6) - subgraph & pathways")
    
    if not os.path.exists(ppi_graph_path):
        print("Graph file not found. Skipping case study.")
        return None

    # Load network
    G = nx.read_graphml(ppi_graph_path)
    
    # --- CORREÇÃO: Mapear Gene Symbol -> Node ID ---
    # O grafo usa IDs do STRING (9606.ENSP...), mas o alvo é "IL6".
    # Criamos um mapa reverso para encontrar o ID correto.
    symbol_to_id = {}
    for node, data in G.nodes(data=True):
        # Tenta pegar 'gene_symbol', se falhar tenta 'name', se falhar usa o próprio ID
        symbol = data.get('gene_symbol')
        if symbol:
            symbol_to_id[str(symbol).upper()] = node
    
    # Load targets
    drug_targets = pd.read_csv(drug_targets_path)
    drug_targets['drug_name'] = normalize_drug_names(drug_targets['drug_name'])
    
    # get siltuximab targets
    siltuximab_targets = drug_targets[drug_targets['drug_name'] == 'SILTUXIMAB']['target_gene'].tolist()
    if len(siltuximab_targets) == 0:
        print("No Siltuximab targets found in mapping. Skipping case study.")
        return None
        
    target_symbol = siltuximab_targets[0].upper() # Ex: "IL6"
    
    # Busca o ID real do nó usando o mapa
    target_node_id = symbol_to_id.get(target_symbol)
    
    if target_node_id and target_node_id in G:
        print(f"Target found: {target_symbol} -> Node ID: {target_node_id}")
        
        # Build 2-hop subgraph usando o ID correto
        subgraph_nodes = set([target_node_id])
        for neighbor in G.neighbors(target_node_id):
            subgraph_nodes.add(neighbor)
        for neighbor in list(subgraph_nodes):
            for second_neighbor in G.neighbors(neighbor):
                subgraph_nodes.add(second_neighbor)
        subgraph = G.subgraph(subgraph_nodes).copy()
        
        print(f"Subgraph size: {subgraph.number_of_nodes()} nodes, {subgraph.number_of_edges()} edges")
        
        # Visualize
        plt.figure(figsize=(12, 10))
        pos = nx.spring_layout(subgraph, k=0.5, iterations=50, seed=RANDOM_SEED)
        
        # Cores: Destacar o alvo
        node_colors = ['red' if node == target_node_id else 'lightblue' for node in subgraph.nodes()]
        node_sizes = [3000 if node == target_node_id else 500 for node in subgraph.nodes()]
        
        # Labels: Usar o gene_symbol para plotar, não o ID feio do STRING
        labels = {}
        for node in subgraph.nodes():
            labels[node] = G.nodes[node].get('gene_symbol', str(node))

        nx.draw_networkx_nodes(subgraph, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8)
        nx.draw_networkx_edges(subgraph, pos, alpha=0.3)
        nx.draw_networkx_labels(subgraph, pos, labels=labels, font_size=8, font_weight='bold')
        
        plt.title(f"Subgrafo 2-hop ao redor do alvo ({target_symbol})", fontsize=14)
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(RESULTS_DIR / "siltuximab_subgraph.png", dpi=300)
        plt.close()
        
        # Stats básicas
        deg = nx.degree_centrality(G).get(target_node_id, 0)
        bet = nx.betweenness_centrality(G).get(target_node_id, 0) # Cuidado: isso pode demorar no grafo inteiro
        # Se betweenness demorar muito, remova a linha acima ou calcule apenas no subgrafo (embora seja metricamente diferente)
        
        print(f"Stats for {target_symbol}: Degree Centrality={deg:.4f}")
        return {'target': target_symbol, 'nodes': subgraph.number_of_nodes()}
        
    else:
        print(f"Target node {target_symbol} (ID: {target_node_id}) not found in graph keys.")
        return None


# ---------- SECTION 8: ROC & PRECISION-RECALL CURVES ----------

def roc_pr_analysis(netcandrug_path="data/processed/final_drug_ranking.csv",
                    baseline_expr_path=RESULTS_DIR / "baseline_expression_ranking.csv"):
    print("\n[ROC/PR] Computing ROC and Precision-Recall curves...")
    
    # Load and clean NetCanDrug
    netcandrug = pd.read_csv(netcandrug_path)
    netcandrug['drug_name'] = normalize_drug_names(netcandrug['drug_name'])
    netcandrug = recalculate_score_without_clinical(netcandrug)
    baseline_expr = pd.read_csv(baseline_expr_path)
    nc_score_col = 'Validation_Score'
    for df in [netcandrug, baseline_expr]:
        if 'drug_name' in df.columns:
            df['drug_name'] = normalize_drug_names(df['drug_name'])
    # Prepare ground truth
    fda_set = FDA_BREAST_DRUGS
    def prepare_df(df, score_col):
        df = df.copy()
        df['y_true'] = df['drug_name'].isin(fda_set).astype(int)
        df['y_score'] = df[score_col]
        return df[['drug_name','y_true','y_score']].dropna()
    # Select score column names (attempt common variants)
    nc_score_col = None
    for candidate in ['NetCanDrug_Score', 'netcandrug_score', 'score']:
        if candidate in netcandrug.columns:
            nc_score_col = candidate
            break
    bl_score_col = None
    for candidate in ['baseline_expression_score', 'baseline_expression_score', 'expr_score']:
        if candidate in baseline_expr.columns:
            bl_score_col = candidate
            break
    if nc_score_col is None or bl_score_col is None:
        print("Could not locate expected score columns in the provided ranking files. Aborting ROC/PR analysis.")
        return None
    nc_data = prepare_df(netcandrug, nc_score_col)
    bl_data = prepare_df(baseline_expr, bl_score_col)
    # Align on common drugs
    common_drugs = set(nc_data['drug_name']).intersection(set(bl_data['drug_name']))
    nc_al = nc_data[nc_data['drug_name'].isin(common_drugs)]
    bl_al = bl_data[bl_data['drug_name'].isin(common_drugs)]
    # ROC and PR curves
    fpr_nc, tpr_nc, _ = roc_curve(nc_al['y_true'], nc_al['y_score'])
    auc_nc = auc(fpr_nc, tpr_nc)
    fpr_bl, tpr_bl, _ = roc_curve(bl_al['y_true'], bl_al['y_score'])
    auc_bl = auc(fpr_bl, tpr_bl)
    precision_nc, recall_nc, _ = precision_recall_curve(nc_al['y_true'], nc_al['y_score'])
    ap_nc = average_precision_score(nc_al['y_true'], nc_al['y_score'])
    precision_bl, recall_bl, _ = precision_recall_curve(bl_al['y_true'], bl_al['y_score'])
    ap_bl = average_precision_score(bl_al['y_true'], bl_al['y_score'])
    
    # --- INICIO DA MODIFICAÇÃO BOOTSTRAP ---
    n_boot = 2000
    boot_aucs = []
    boot_aps = []
    
    # Converter para array numpy evita erros no resample
    y_true_arr = nc_al['y_true'].values
    y_score_arr = nc_al['y_score'].values
    
    rng = np.random.RandomState(42) # Seed fixa para o bootstrap
    
    for i in range(n_boot):
        # Reamostragem com reposição (mantendo proporção de classes se possível, ou simples)
        y_true_b, y_score_b = resample(y_true_arr, y_score_arr, replace=True, random_state=rng)
        
        # Só calcula se houver pelo menos uma classe positiva e uma negativa
        if len(np.unique(y_true_b)) < 2:
            continue
            
        # Calcula métricas da amostra
        boot_aucs.append(roc_auc_score(y_true_b, y_score_b))
        boot_aps.append(average_precision_score(y_true_b, y_score_b))
        
    # Calcular Intervalos de Confiança (2.5% e 97.5%)
    auc_lower, auc_upper = np.percentile(boot_aucs, [2.5, 97.5])
    ap_lower, ap_upper = np.percentile(boot_aps, [2.5, 97.5])
    # --- FIM DA MODIFICAÇÃO ---
    
    # Plot
    # Plot ROC
    plt.figure(figsize=(7, 6))
    ax = plt.gca()
    ax.plot(fpr_nc, tpr_nc, label=f'NetCanDrug (AUC={auc_nc:.3f})')
    ax.plot(fpr_bl, tpr_bl, label=f'Baseline (AUC={auc_bl:.3f})')
    ax.plot([0,1],[0,1],'k--', label='Random')
    ax.set_xlabel('False Positive Rate'); ax.set_ylabel('True Positive Rate'); ax.set_title('ROC Curve')
    ax.legend(loc='lower right'); ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "roc_curve.png", dpi=300)
    plt.close()

    # Plot PR
    plt.figure(figsize=(7, 6))
    ax = plt.gca()
    baseline_precision = nc_al['y_true'].sum() / len(nc_al) if len(nc_al)>0 else 0.0
    ax.plot(recall_nc, precision_nc, label=f'NetCanDrug (AP={ap_nc:.3f})')
    ax.plot(recall_bl, precision_bl, label=f'Baseline (AP={ap_bl:.3f})')
    ax.plot([0,1],[baseline_precision, baseline_precision],'k--', label=f'Random (AP={baseline_precision:.3f})')
    ax.set_xlabel('Recall'); ax.set_ylabel('Precision'); ax.set_title('Precision-Recall Curve')
    ax.legend(loc='upper right'); ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "pr_curve.png", dpi=300)
    plt.close()
    # Print summary metrics
    print("\nROC / PR summary:")
    print(f"NetCanDrug AUC: {auc_nc:.3f} (95% CI: {auc_lower:.3f}-{auc_upper:.3f})")
    print(f"NetCanDrug AP:  {ap_nc:.3f} (95% CI: {ap_lower:.3f}-{ap_upper:.3f})")
    print(f"Baseline:    AUC={auc_bl:.4f}, AP={ap_bl:.4f}")
    print(f"Improvement AUC: {((auc_nc/auc_bl - 1)*100) if auc_bl>0 else np.nan:.1f}%")
    print(f"Improvement AP:  {((ap_nc/ap_bl - 1)*100) if ap_bl>0 else np.nan:.1f}%")
    return {'auc_nc': auc_nc, 'ap_nc': ap_nc, 'auc_bl': auc_bl, 'ap_bl': ap_bl}


# ---------- SECTION 9: DISTRIBUTIONS & STATISTICAL TESTS ----------
def distribution_and_statistics(final_ranking_path=config.FINAL_RANKING_FILE):
    print("\n[STATS] Distribution analysis (Original Colors + Fixed Text)")
    df = pd.read_csv(final_ranking_path)
    df['drug_name'] = normalize_drug_names(df.get('drug_name', pd.Series(dtype=str)))
    df['is_fda_approved'] = df['drug_name'].isin(FDA_BREAST_DRUGS)
    
    # Recalcular score limpo para plotagem
    df = recalculate_score_without_clinical(df)
    score_col = 'Validation_Score' # Usar esta variável no plot
    
    # Ensure component columns
    components = ['topology_score_norm', 'pathway_score_norm', 'expression_score_norm']
    for comp in components:
        ensure_df_column(df, comp, np.nan)
        
    # --- STYLE: Original Colors (Blue/Orange) ---
    c_other = 'tab:blue'   # Azul para "Outros"
    c_fda = 'tab:orange'   # Laranja para "FDA"
    
    # --- A: Histograma ---
    plt.figure(figsize=(7, 5))
    ax = plt.gca()
    # Substituir 'NetCanDrug_Score' por score_col
    other_scores = df[~df['is_fda_approved']][score_col].dropna()
    fda_scores = df[df['is_fda_approved']][score_col].dropna()
    ax.hist(other_scores, bins=50, alpha=0.7, label='Other drugs', color=c_other)
    ax.hist(fda_scores, bins=20, alpha=0.9, label='FDA-approved', color=c_fda)
    
    perc_95 = df[score_col].quantile(0.95)
    ax.axvline(perc_95, linestyle='--', color='black', label='95th percentile')
        
    ax.set_xlabel('NetCanDrug Score')
    ax.set_ylabel('Frequency')
    ax.set_title('(A) Distribution of NetCanDrug Score', fontweight='bold')
    ax.legend(loc='best', framealpha=0.9) 
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "score_distribution_histogram.png", dpi=300)
    plt.close()
    
    # --- B: Boxplot ---
    plt.figure(figsize=(7, 5))
    ax = plt.gca()
    df_melted = df.melt(id_vars=['is_fda_approved'], value_vars=components, var_name='Component', value_name='Score')
    
    # Paleta explícita para garantir Azul/Laranja
    my_pal = {False: c_other, True: c_fda}
    sns.boxplot(data=df_melted, x='Component', y='Score', hue='is_fda_approved', ax=ax, palette=my_pal)
    
    ax.set_xticklabels(['Topology','Pathway','Expression'])
    ax.set_ylabel('Normalized Score')
    ax.set_title('(B) Score Components: FDA vs Others', fontweight='bold')
    
    # Legenda corrigida
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, ['Other', 'Approved'], title='FDA Status', loc='upper right', framealpha=0.9)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "score_components_boxplot.png", dpi=300)
    plt.close()
    
    # --- C: Mann-Whitney (Barplot) ---
    plt.figure(figsize=(7, 5))
    ax = plt.gca()
    p_values = []
    comp_names = ['Topology','Pathway','Expression']
    
    for comp in components:
        fda_scores = df[df['is_fda_approved']][comp].dropna()
        other_scores = df[~df['is_fda_approved']][comp].dropna()
        if len(fda_scores)>0 and len(other_scores)>0:
            statistic, p_val = stats.mannwhitneyu(fda_scores, other_scores, alternative='greater')
        else:
            p_val = 1.0
        p_values.append(p_val)
    
    bar_lengths = [-np.log10(p) if p>0 else 0 for p in p_values]
    
    # Cores de significância (Verde/Vermelho) mantidas aqui pois é padrão para p-valor
    colors = ['green' if p < 0.05 else 'red' for p in p_values]
    
    bars = ax.barh(comp_names, bar_lengths, color=colors, alpha=0.7)
    
    ax.axvline(-np.log10(0.05), color='black', linestyle='--', label='p=0.05')
    ax.set_xlabel('-log10(p-value)')
    ax.set_title('(C) Mann-Whitney U test (FDA > Others)', fontweight='bold')
    
    # --- FIX CRÍTICO: Espaço para texto ---
    max_len = max(bar_lengths) if bar_lengths else 1
    ax.set_xlim(0, max_len * 1.5) # Muito espaço extra à direita
    
    for i, (length, pval) in enumerate(zip(bar_lengths, p_values)):
        # Se a barra for minúscula (p não significativo), empurra o texto pra direita
        x_pos = length + (max_len * 0.02)
        if length < 0.2: 
             x_pos = max_len * 0.1 # Garante que não encoste no eixo Y
             
        ax.text(x_pos, i, f"p={pval:.2e}", va='center', fontsize=10, fontweight='bold')
    
    ax.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "mann_whitney_test.png", dpi=300)
    plt.close()

    # --- D: CDF ---
    plt.figure(figsize=(7, 5))
    ax = plt.gca()
    # Substituir 'NetCanDrug_Score' por score_col
    fda_scores_sorted = np.sort(df[df['is_fda_approved']][score_col].dropna())
    other_scores_sorted = np.sort(df[~df['is_fda_approved']][score_col].dropna())
    
    if len(fda_scores_sorted)>0:
        ax.plot(fda_scores_sorted, np.linspace(0,1,len(fda_scores_sorted)), label='FDA-approved', color=c_fda, linewidth=2.5)
    if len(other_scores_sorted)>0:
        ax.plot(other_scores_sorted, np.linspace(0,1,len(other_scores_sorted)), label='Other drugs', color=c_other, linewidth=2.5)
        
    ax.set_xlabel('NetCanDrug Score')
    ax.set_ylabel('Cumulative Probability')
    ax.set_title('(D) Cumulative Distribution Function', fontweight='bold')
    ax.legend(loc='lower right')
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "cdf_plot.png", dpi=300)
    plt.close()
    
    return {'p_values': p_values}


# ---------- SECTION 10: HYPERGEOMETRIC & FDA VALIDATION ----------

def hypergeom_fda_validation(final_ranking_path=config.FINAL_RANKING_FILE):
    print("\n[HYPERGEOMETRIC] Validation vs FULL FDA-approved list")
    
    df = pd.read_csv(final_ranking_path)
    # Recalcula se necessário (conforme discutimos antes)
    if 'Validation_Score' not in df.columns:
        df = recalculate_score_without_clinical(df)
        
    df['drug_name'] = normalize_drug_names(df.get('drug_name', pd.Series(dtype=str)))
    df = df.sort_values('Validation_Score', ascending=False)
    
    total_drugs = df.shape[0]
    total_fda = len(FDA_BREAST_DRUGS) # Conjunto COMPLETO
    
    results = []
    for k in [10, 20, 30, 50, 100]:
        top_k = df.head(k)['drug_name'].tolist()
        hits = [d for d in top_k if d in FDA_BREAST_DRUGS]
        hits_count = len(hits)
        precision = hits_count / k
        recall = hits_count / total_fda if total_fda > 0 else 0
        
        # Teste estatístico
        pval = hypergeom.sf(hits_count - 1, total_drugs, total_fda, k)
        
        results.append({
            'Metric_Type': 'Full FDA List',  # <--- NOVA COLUNA
            'K': k, 
            'Precision': round(precision, 4), 
            'Recall': round(recall, 4),
            'Hits': hits_count, 
            'p_value': f"{pval:.2e}"
        })
        
    results_df = pd.DataFrame(results)
    results_df.to_csv(RESULTS_DIR / "validation_full_fda.csv", index=False)
    print(results_df.to_string(index=False))
    return results_df


def check_seed_sensitivity(df_ranking, fda_list, n_seeds=50):
    print(f"\n[ROBUSTNESS] Checking stability across {n_seeds} random splits...")
    precisions = []
    
    # Garante que temos o score final
    df = df_ranking.copy()
    if 'Validation_Score' not in df.columns:
        df = recalculate_score_without_clinical(df)
    df = df.sort_values('Validation_Score', ascending=False)
    
    for i in range(n_seeds):
        # Muda a seed do split a cada iteração (i)
        _, test_drugs_iter = train_test_split(fda_list, test_size=0.4, random_state=i)
        
        # Normaliza nomes
        test_set = set([d.upper().replace(r'[^A-Z0-9]', '') for d in test_drugs_iter])
        
        # Calcula Precision@20 para esse split específico
        top_20 = df.head(20)
        hits = top_20['drug_name'].isin(test_set).sum()
        precisions.append(hits / 20)
        
    mean_p = np.mean(precisions)
    std_p = np.std(precisions)
    
    print(f"Precision@20 Stability: {mean_p:.3f} ± {std_p:.3f}")
    print(f"Min: {np.min(precisions):.3f}, Max: {np.max(precisions):.3f}")
    return mean_p, std_p


# ---------- MAIN ----------

def main():
    print("="*80)
    print("NETCANDRUG CONSOLIDATED VALIDATION SUITE")
    print("="*80)
    
    # Load core data
    df_ranking, df_drug_targets, df_deg = load_data()
    
    # Check seed sensitivity for robustness
    check_seed_sensitivity(df_ranking, FDA_BREAST_DRUGS_LIST)
    
    # Validation (train/test held-out)
    train_drugs, test_drugs, test_results_df = train_test_validation(df_ranking, FDA_BREAST_DRUGS_LIST)
    
    # Ablation (uses df_ranking and is_test_fda column from earlier)
    # ensure is_test_fda exists
    if 'is_test_fda' not in df_ranking.columns:
        df_ranking['is_test_fda'] = df_ranking['drug_name'].isin(test_drugs)
    ablation_df = ablation_study(df_ranking)
    
    # Baseline expression ranking and baseline precision on test set
    baseline_expr_df, baseline_precision_20 = compute_baseline_expression(df_drug_targets, df_deg, test_drugs=[d.upper() for d in test_drugs])
    
    # Use the unified permutation test function
    perm_results_heldout = run_permutation_test(
        df_ranking, 
        target_drug_list=test_drugs,
        mode_label='Held-Out',
        n_iter=10000,
        top_k=20
    )
    
    # Extract values for backwards compatibility
    mean_random = perm_results_heldout.iloc[0]['Random_Mean']
    std_random = perm_results_heldout.iloc[0]['Random_Std']
    netcandrug_p20 = perm_results_heldout.iloc[0]['Model_Precision']
    p_value_empirical = perm_results_heldout.iloc[0]['Empirical_p_value']
    
    print(f"\nRandom mean: {mean_random:.4f} ± {std_random:.4f}")
    print(f"NetCanDrug Precision@20: {netcandrug_p20:.3f}")
    print(f"Fold improvement: {netcandrug_p20/mean_random:.1f}x")
    print(f"Empirical p-value: {p_value_empirical:.4f}")
    
    # Final comparison table
    comparison_df = generate_final_comparison(netcandrug_p20, baseline_precision_20 if baseline_precision_20 is not None else 0.0, mean_random)
    
    # Save validation report text
    report_path = RESULTS_DIR / "validation_report.txt"
    with open(report_path, "w") as f:
        f.write("="*80 + "\n")
        f.write("NETCANDRUG VALIDATION REPORT\n")
        f.write("="*80 + "\n\n")
        f.write("1. TEST SET VALIDATION (Unseen Data)\n")
        f.write(test_results_df.to_string(index=False) + "\n\n")
        f.write("2. ABLATION STUDY (Component Contribution)\n")
        f.write(ablation_df.to_string(index=False) + "\n\n")
        f.write("3. BASELINE COMPARISON\n")
        f.write(comparison_df.to_string(index=False) + "\n\n")
        f.write("4. STATISTICAL SIGNIFICANCE\n")
        f.write(f"   Empirical p-value (permutation test): {p_value_empirical:.4f}\n")
        f.write(f"   Fold improvement over random: {netcandrug_p20/mean_random:.1f}x\n")
        if baseline_precision_20 and baseline_precision_20>0:
            f.write(f"   Fold improvement over baseline: {netcandrug_p20/baseline_precision_20:.1f}x\n\n")
        f.write("="*80 + "\n")
    print("\nReport saved to:", report_path)
    
    # Case study
    try:
        siltuximab_stats = siltuximab_case_study()
    except Exception as e:
        print("Siltuximab case study failed:", e)
        siltuximab_stats = None
    
    # ROC / PR (now includes bootstrap analysis)
    try:
        roc_stats = roc_pr_analysis()
    except Exception as e:
        print("ROC/PR analysis failed:", e)
        roc_stats = None
    
    # Distributional stats
    try:
        stats_summary = distribution_and_statistics()
    except Exception as e:
        print("Distributional analysis failed:", e)
        stats_summary = None
    
    # Hypergeom FDA validation
    try:
        hyper_df = hypergeom_fda_validation()
    except Exception as e:
        print("Hypergeom FDA validation failed:", e)
        hyper_df = None
    
    # Run permutation test for Full FDA list
    try:
        perm_results_full_fda = run_permutation_test(
            df_ranking,
            target_drug_list=FDA_BREAST_DRUGS_LIST,
            mode_label='Full_FDA',
            n_iter=10000,
            top_k=20
        )
        
        # Save all permutation results
        all_perm_results = pd.concat([perm_results_heldout, perm_results_full_fda], ignore_index=True)
        all_perm_results.to_csv(RESULTS_DIR / "permutation_tests.csv", index=False)
        print("Permutation test results saved to permutation_tests.csv")
        
    except Exception as e:
        print("Full FDA permutation test failed:", e)
    
    print("\n✅ All requested analyses attempted. Check the results/ directory for CSVs, PNGs and the report.")
    print("="*80)


if __name__ == "__main__":
    main()

