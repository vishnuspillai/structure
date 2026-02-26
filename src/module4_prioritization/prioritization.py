import pandas as pd
import numpy as np
import os

df = pd.read_csv('data/processed/chrna7_missense_spatial_annotated.csv')

def compute_score(row):
    score = 0
    # Population
    af = row.get('AF', np.nan)
    if pd.notna(af):
        if af < 1e-5:
            score += 2
        elif 1e-5 <= af < 1e-4:
            score += 1

    # Structural
    struct_score = 0
    if row.get('is_binding_site', False):
        struct_score += 4
    # User mentioned is_pore_core, checking for both possible names from earlier phrasing
    if row.get('is_pore_core', False) or row.get('is_pore_region', False):
        struct_score += 4
    if row.get('is_interface', False):
        struct_score += 2
    
    is_pore = row.get('is_pore_core', False) or row.get('is_pore_region', False)
    if row.get('is_tm_core', False) and not is_pore:
        struct_score += 1
    
    # Cap total structural points
    struct_score = min(struct_score, 4)
    score += struct_score

    # Functional
    cadd = row.get('cadd_phred', np.nan)
    if pd.notna(cadd):
        if cadd >= 30:
            score += 3
        elif 25 <= cadd < 30:
            score += 2
        elif 20 <= cadd < 25:
            score += 1
            
    polyphen = str(row.get('polyphen_pred', '')).lower()
    if 'probably_damaging' in polyphen:
        score += 2
    elif 'possibly_damaging' in polyphen:
        score += 1
        
    sift = str(row.get('sift_pred', '')).lower()
    if 'deleterious' in sift:
        score += 1
        
    return score

df['priority_score'] = df.apply(compute_score, axis=1)

def get_category(s):
    if s >= 10: return 'High'
    if 6 <= s <= 9: return 'Medium'
    if 3 <= s <= 5: return 'Low'
    return 'Very Low'

df['priority_category'] = df['priority_score'].apply(get_category)

df = df.sort_values(by='priority_score', ascending=False)
output_path = 'data/processed/chrna7_ranked_variants.csv'
df.to_csv(output_path, index=False)

output_txt = 'mod4_output.txt'
with open(output_txt, 'w', encoding='utf-8') as f:
    f.write("Score distribution summary:\n")
    f.write(df['priority_score'].describe().to_string() + "\n")
    f.write("\nCount per priority_category:\n")
    f.write(df['priority_category'].value_counts().to_string() + "\n")
    f.write("\nTop 15 variants (rsid, position, domain, score):\n")
    f.write(df[['rsid', 'protein_position', 'domain_region', 'priority_score']].head(15).to_string(index=False) + "\n")

