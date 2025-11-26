import pandas as pd
import numpy as np
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config   

#---------------------------------------------------
# 1. CARREGAR DADOS
#---------------------------------------------------
final_ranking_file = config.FINAL_RANKING_FILE
drug_target_mapping_file = config.DRUG_TARGET_FILE
deg_list_file = config.DEG_FILE
df = pd.read_csv(final_ranking_file)

drug_targets = pd.read_csv(drug_target_mapping_file)
deg = pd.read_csv(deg_list_file)

# Normalizar
drug_targets["drug_name"] = (
    drug_targets["drug_name"]
    .astype(str).str.upper().str.strip()
    .str.replace(r'[^A-Z0-9]', '', regex=True)
)

deg["gene_name"] = deg["gene_name"].astype(str).str.upper()

# Merge
merged = drug_targets.merge(
    deg, 
    left_on="target_gene", 
    right_on="gene_name",
    how="inner"
)

# Remover NaN
merged = merged.dropna(subset=['logFC'])

# Calcular score baseline
baseline_expr = (
    merged.groupby("drug_name")["logFC"]
    .apply(lambda x: np.mean(np.abs(x)))
    .reset_index()
    .rename(columns={"logFC": "baseline_expression_score"})
)

# Ordenar
baseline_expr = baseline_expr.sort_values("baseline_expression_score", ascending=False)

# Salvar
baseline_expr.to_csv("results/baseline_expression_ranking.csv", index=False)

print("\n" + "="*70)
print("BASELINE RANKING BY EXPRESSION (|logFC|)")
print("="*70)
print(baseline_expr.head(20).to_string(index=False))
print("\nFile saved: results/baseline_expression_ranking.csv")
print("="*70)

#---------------------------------------------------
# 2. VALIDAÇÃO COM LISTA FDA (COM SINÔNIMOS)
#---------------------------------------------------
print("\n" + "="*70)
print("VALIDATING BASELINE AGAINST FDA-APPROVED DRUGS")
print("="*70)

# Lista FDA
# Lista FDA correta para câncer de mama (BRCA)
fda_brca_drugs = {
    # Quimioterapia
    'PACLITAXEL', 'DOCETAXEL', 'DOXORUBICIN', 'EPIRUBICIN',
    'CYCLOPHOSPHAMIDE', 'FLUOROURACIL', 'CAPECITABINE',
    'CARBOPLATIN', 'CISPLATIN', 'GEMCITABINE', 'VINORELBINE',

    # Terapias alvo (HER2)
    'TRASTUZUMAB', 'PERTUZUMAB', 'ADOTRASTUZUMAB', 
    'LAPATINIB', 'NERATINIB', 'TUCATINIB', 'MARGETUXIMAB',

    # Inibidores CDK4/6
    'PALBOCICLIB', 'RIBOCICLIB', 'ABEMACICLIB',

    # Terapia hormonal
    'TAMOXIFEN', 'LETROZOLE', 'ANASTROZOLE', 'EXEMESTANE',
    'FULVESTRANT', 'TOREMIFENE',

    # PI3K / mTOR
    'EVEROLIMUS', 'ALPELISIB',

    # Inibidores PARP
    'OLAPARIB', 'TALAZOPARIB',

    # Imunoterapia
    'PEMBROLIZUMAB', 'ATEZOLIZUMAB',

    # Sinônimo válido (apenas este)
    '5FU',   # variante de FLUOROURACIL
}


drug_synonyms = {
    # Sinônimos corretos de CISPLATIN
    'CISDIAMMINEDICHLOROPLATINUM': 'CISPLATIN',
    'CISDIAMMINEDICHLOROPLATINUMII': 'CISPLATIN',
    'CISDICHLORODIAMMINEPLATINUM': 'CISPLATIN',

    # Sinônimos corretos para Fluorouracil / 5FU
    '5FLUOROURACIL': '5FU',
    'FLUOROURACIL': '5FU',
}

# Aplicar sinônimos
baseline_expr["drug_name_mapped"] = baseline_expr["drug_name"].replace(drug_synonyms)

# Calcular Precision@K
K_VALUES = [10, 20, 30, 50, 100]
baseline_results = []

for k in K_VALUES:
    top_k = baseline_expr.head(k)["drug_name_mapped"].tolist()
    hits = [d for d in top_k if d in fda_brca_drugs]
    precision = len(hits) / k
    
    baseline_results.append({
        "K": k,
        "Precision": round(precision, 4),
        "Hits": len(hits),
        "Hit_Drugs": ", ".join(hits[:5])
    })

baseline_results_df = pd.DataFrame(baseline_results)
baseline_results_df.to_csv("results/baseline_expression_validation.csv", index=False)

print("\nBaseline Expression - Validation Results:")
print(baseline_results_df.to_string(index=False))
print("\nFile saved: results/baseline_expression_validation.csv")

#---------------------------------------------------
# 3. COMPARAÇÃO COM NETCANDRUG
#---------------------------------------------------
print("\n" + "="*70)
print("COMPARISON: NetCanDrug vs Baseline Expression")
print("="*70)

netcandrug_precision_20 = 0.15
baseline_precision_20 = baseline_results_df[baseline_results_df['K'] == 20]['Precision'].values[0]

if baseline_precision_20 > 0:
    improvement = netcandrug_precision_20 / baseline_precision_20
    improvement_str = f"{improvement:.2f}x"
else:
    improvement_str = "∞ (infinite)"

comparison = pd.DataFrame({
    "Method": [
        "NetCanDrug (Multi-component)", 
        "Baseline (Expression only)", 
        "Improvement Factor"
    ],
    "Precision@20": [
        f"{netcandrug_precision_20:.4f}", 
        f"{baseline_precision_20:.4f}",
        improvement_str
    ]
})

print("\n")
print(comparison.to_string(index=False))
print("\n" + "="*70)

comparison.to_csv("results/method_comparison.csv", index=False)
print("Comparison saved: results/method_comparison.csv")
print("="*70)