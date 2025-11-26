import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config   


#---------------------------------------------------
# 1. CARREGAR SEU RANKING
#---------------------------------------------------

# ⛔ Ajuste aqui o nome do seu arquivo CSV
final_ranking_file = config.FINAL_RANKING_FILE
df = pd.read_csv(final_ranking_file)
drug_col = "drug_name"

# Normalização ROBUSTA
df[drug_col] = (
    df[drug_col]
    .astype(str)
    .str.upper()
    .str.strip()
    .str.replace(r'[^A-Z0-9]', '', regex=True)  # Remove TUDO exceto letras/números
)

# Lista FDA
fda_brca_drugs = {
    'PACLITAXEL', 'DOCETAXEL', 'DOXORUBICIN', 'EPIRUBICIN',
    'CYCLOPHOSPHAMIDE', 'FLUOROURACIL', 'CAPECITABINE',
    'CARBOPLATIN', 'CISPLATIN', 'GEMCITABINE', 'VINORELBINE',
    'TRASTUZUMAB', 'PERTUZUMAB', 'ADOTRASTUZUMAB', 'LAPATINIB',
    'NERATINIB', 'TUCATINIB', 'MARGETUXIMAB',
    'PALBOCICLIB', 'RIBOCICLIB', 'ABEMACICLIB',
    'TAMOXIFEN', 'LETROZOLE', 'ANASTROZOLE', 'EXEMESTANE',
    'FULVESTRANT', 'TOREMIFENE',
    'EVEROLIMUS', 'ALPELISIB',
    'OLAPARIB', 'TALAZOPARIB',
    'PEMBROLIZUMAB', 'ATEZOLIZUMAB'
}

# Funções
def precision_at_k(df, k):
    top_k = df.head(k)[drug_col].tolist()
    hits = [d for d in top_k if d in fda_brca_drugs]
    return len(hits) / k, hits

def recall_at_k(hits):
    return len(hits) / len(fda_brca_drugs)

def f1_score(precision, recall):
    if precision + recall == 0:
        return 0
    return 2 * (precision * recall) / (precision + recall)

def p_value_hypergeom(hits_count, k, total_drugs, total_fda):
    return hypergeom.sf(hits_count - 1, total_drugs, total_fda, k)

# Calcular métricas
K_VALUES = [10, 20, 30, 50, 100]
results = []
total_drugs = df.shape[0]
total_fda = len(fda_brca_drugs)

for k in K_VALUES:
    precision, hits = precision_at_k(df, k)
    recall = recall_at_k(hits)
    f1 = f1_score(precision, recall)
    pval = p_value_hypergeom(len(hits), k, total_drugs, total_fda)
    
    results.append({
        "K": k,
        "Precision": round(precision, 4),
        "Recall": round(recall, 4),
        "F1": round(f1, 4),
        "Hits": len(hits),
        "p_value": f"{pval:.2e}",
        "Hit_Drugs": ", ".join(hits[:5])  # Mostra só os primeiros 5
    })

# Exportar
results_df = pd.DataFrame(results)
results_df.to_csv("results/validation_fda_results.csv", index=False)

print("\n" + "="*70)
print("VALIDATION RESULTS - NetCanDrug vs FDA-Approved Drugs")
print("="*70)
print(results_df.to_string(index=False))
print("\nFile saved: results/validation_fda_results.csv")
print("="*70)