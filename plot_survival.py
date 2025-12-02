import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

# Carregar ranking
df = pd.read_csv("data/processed/final_drug_ranking.csv")

# Lista FDA (seus 13 fármacos de teste)
fda_drugs = ['PACLITAXEL', 'DOCETAXEL', 'DOXORUBICIN', 'TRASTUZUMAB', 
             'TAMOXIFEN', 'LETROZOLE', 'CISPLATIN', 'GEMCITABINE',
             'PALBOCICLIB', 'EVEROLIMUS', 'OLAPARIB', 'PEMBROLIZUMAB',
             'ATEZOLIZUMAB']

# Normalizar
df['drug_name_clean'] = df['drug_name'].str.upper().str.strip()

# Separar FDA vs outros
fda_mask = df['drug_name_clean'].isin(fda_drugs)
fda_scores = df[fda_mask]
other_scores = df[~fda_mask]

# Calcular medianas
print("📊 TOPOLOGY SCORE:")
print(f"   Mediana FDA: {fda_scores['topology_score_norm'].median():.3f}")
print(f"   Mediana Outros: {other_scores['topology_score_norm'].median():.3f}")

# Teste Mann-Whitney
statistic, pvalue = mannwhitneyu(fda_scores['topology_score_norm'], 
                                  other_scores['topology_score_norm'],
                                  alternative='greater')
print(f"   Mann-Whitney p-value: {pvalue:.2e}")

print("\n📊 PATHWAY SCORE:")
print(f"   Mediana FDA: {fda_scores['pathway_score_norm'].median():.3f}")
print(f"   Mediana Outros: {other_scores['pathway_score_norm'].median():.3f}")

statistic, pvalue = mannwhitneyu(fda_scores['pathway_score_norm'], 
                                  other_scores['pathway_score_norm'],
                                  alternative='greater')
print(f"   Mann-Whitney p-value: {pvalue:.2e}")

print("\n📊 EXPRESSION SCORE:")
print(f"   Mediana FDA: {fda_scores['expression_score_norm'].median():.3f}")
print(f"   Mediana Outros: {other_scores['expression_score_norm'].median():.3f}")

statistic, pvalue = mannwhitneyu(fda_scores['expression_score_norm'], 
                                  other_scores['expression_score_norm'],
                                  alternative='greater')
print(f"   Mann-Whitney p-value: {pvalue:.2e}")