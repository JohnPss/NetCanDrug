# ===================================================================
# NETCANDRUG PROJECT - SCRIPT 05: DRUG-TARGET MAPPING
# ===================================================================
#
# OBJECTIVE:
# 1. Load the Drug-Gene interactions database (DGIdb).
# 2. Load the list of proteins from our tumor network.
# 3. Filter interactions to keep only approved drugs
#    and targets present in our network.
# 4. Clean data to remove IDs, supplements and normalize names.
# 5. Save the final Drug-Target mapping.
#
# ===================================================================

import pandas as pd
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
print("SCRIPT 05: DRUG-TARGET MAPPING")
print("=" * 80)

# --- 1. Path Configuration ---
dgidb_file = config.DGIDB_FILE
topology_file = config.TOPOLOGY_FILE
output_file = config.DRUG_TARGET_FILE

# --- 2. Load Base Data ---
print("\n--- STEP 1: Loading base data ---")

topology_df = pd.read_csv(topology_file)
tumor_network_genes = set(topology_df['gene_symbol'])
print(f"  ✅ Loaded {len(tumor_network_genes)} genes from tumor network.")

print("  - Loading DGIdb file (interactions.tsv)...")
dgidb_df = pd.read_csv(dgidb_file, sep='\t')
print(f"  ✅ Loaded {len(dgidb_df)} raw interactions from DGIdb.")

# --- 3. Apply Double Filter ---
print("\n--- STEP 2: Applying double filter on interactions ---")

approval_col = None
if 'drug_approval' in dgidb_df.columns:
    approval_col = 'drug_approval'
elif 'approvals' in dgidb_df.columns:
    approval_col = 'approvals'
elif 'approved' in dgidb_df.columns:
    approval_col = 'approved'

if approval_col:
    approved_drugs_df = dgidb_df[dgidb_df[approval_col].astype(str).str.contains("True|Approved", na=False, case=False)]
    print(f"  - Filter 1 (FDA-Approved): {len(approved_drugs_df)} interactions retained.")
else:
    print("  ⚠️ WARNING: Drug approval column not found. Skipping FDA filter.")
    approved_drugs_df = dgidb_df

filtered_df = approved_drugs_df[approved_drugs_df['gene_name'].isin(tumor_network_genes)]
print(f"  - Filter 2 (Targets in Network): {len(filtered_df)} final interactions retained.")

# --- 4. Clean, Structure and Save ---
print("\n--- STEP 3: Cleaning and structuring final mapping ---")

if filtered_df.empty:
    print("\n❌ WARNING: No drug-target interactions found after filtering.")
    exit()

interaction_col_name = 'interaction_type' if 'interaction_type' in filtered_df.columns else 'interaction_types'
final_mapping = filtered_df[['drug_claim_name', 'gene_name', interaction_col_name]].copy()
final_mapping.rename(columns={
    'drug_claim_name': 'drug_name',
    'gene_name': 'target_gene',
    interaction_col_name: 'interaction_type'
}, inplace=True)

# ##### DATA CLEANING BLOCK #####
print("\n--- STEP 3.5: Applying advanced data cleaning ---")
initial_drug_count = final_mapping['drug_name'].nunique()
print(f"  - Unique drugs before cleaning: {initial_drug_count}")

# 1. Remove database IDs from drug names
id_patterns = config.DATABASE_ID_PATTERNS
final_mapping = final_mapping[~final_mapping['drug_name'].str.contains(id_patterns, case=False, na=False)]
print(f"  - After removing database IDs: {final_mapping['drug_name'].nunique()} drugs remaining.")

# 2. Remove common supplements and natural products
supplements = config.SUPPLEMENTS_TO_EXCLUDE
pattern = '|'.join(supplements)
final_mapping = final_mapping[~final_mapping['drug_name'].str.contains(pattern, case=False, na=False)]
print(f"  - After removing supplements: {final_mapping['drug_name'].nunique()} drugs remaining.")

# 3. Normalize names to UPPERCASE and remove extra spaces
final_mapping['drug_name'] = final_mapping['drug_name'].str.upper().str.strip()
print(f"  - After normalizing names: {final_mapping['drug_name'].nunique()} drugs remaining.")

# 4. Remove final duplicates (e.g., 'Siltuximab' and 'SILTUXIMAB' are now the same)
final_mapping.drop_duplicates(subset=['drug_name', 'target_gene'], inplace=True)
print(f"  - After removing final duplicates: {final_mapping['drug_name'].nunique()} drugs remaining.")
# ##### END OF CLEANING BLOCK #####

# Ensure output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

final_mapping.to_csv(output_file, index=False)
print(f"\n  ✅ Final clean mapping saved at: {output_file}")

# --- 5. Final Report ---
print("\n" + "=" * 80)
print("DRUG-TARGET MAPPING REPORT (POST-CLEANING)")
print("=" * 80)

unique_drugs = final_mapping['drug_name'].nunique()
unique_targets = final_mapping['target_gene'].nunique()

print(f"\n📊 Final Mapping Statistics:")
print(f"  - Found {unique_drugs} unique drugs (CLEAN AND APPROVED)...")
print(f"  - ...that target {unique_targets} different targets within our tumor network.")

print(f"\n🔝 Top 10 most targeted genes in our network:")
print(final_mapping['target_gene'].value_counts().head(10).to_string())

print("\n" + "=" * 80)
print("SCRIPT 05 COMPLETED SUCCESSFULLY!")
print("=" * 80)
print("\n🎯 NEXT STEPS:")
print("  - The bridge between drugs and our network is clean and correct.")
print("  - The next script (06-07) will focus on pathway analysis.")
print("=" * 80)

