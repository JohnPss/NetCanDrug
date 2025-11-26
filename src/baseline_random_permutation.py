import pandas as pd
import numpy as np
import sys
import os

# Seed para reprodutibilidade
np.random.seed(42)

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config   


#---------------------------------------------------
# 1. CARREGAR SEU RANKING
#---------------------------------------------------

# ⛔ Ajuste aqui o nome do seu arquivo CSV
final_ranking_file = config.FINAL_RANKING_FILE
df = pd.read_csv(final_ranking_file)

drug_col = "drug_name"

df[drug_col] = (
    df[drug_col].astype(str).str.upper().str.strip()
    .str.replace(r'[^A-Z0-9]', '', regex=True)
)

all_drugs = df[drug_col].unique()
total_drugs = len(all_drugs)

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

# Precision do modelo
def calculate_precision_at_k(k=20):
    top_k = df.head(k)[drug_col].tolist()
    hits = [d for d in top_k if d in fda_brca_drugs]
    return len(hits) / k, len(hits)

precision_model, model_hits = calculate_precision_at_k(20)

# Permutação
N_PERM = 1000
random_precisions = []

for i in range(N_PERM):
    shuffled = np.random.choice(all_drugs, 20, replace=False)
    hits = len(set(shuffled).intersection(fda_brca_drugs))
    random_precisions.append(hits / 20)

random_precisions = np.array(random_precisions)

# Calcular resultados
mean_random = random_precisions.mean()
std_random = random_precisions.std()
fold_improvement = precision_model / mean_random if mean_random > 0 else float('inf')

# P-value empírico
p_value_empirical = (random_precisions >= precision_model).sum() / N_PERM

# Salvar
out = pd.DataFrame({
    "Metric": ["Model_Precision@20", "Random_Mean", "Random_STD", 
               "Fold_Improvement", "P_value_Empirical"],
    "Value": [precision_model, mean_random, std_random, 
              fold_improvement, p_value_empirical]
})

out.to_csv("results/baseline_random_permutation.csv", index=False)

print("\n" + "="*70)
print("BASELINE: Random Permutation Test (1000 iterations)")
print("="*70)
print(out.to_string(index=False))
print("\nFile saved: results/baseline_random_permutation.csv")
print("="*70)