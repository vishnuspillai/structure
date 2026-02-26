import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
import numpy as np

df = pd.read_csv('data/processed/chrna7_ranked_variants.csv')

features = {
    'is_binding_site': 'is_binding_site',
    'is_pore_core': 'is_pore_region',
    'is_interface': 'is_interface'
}

high_mask = df['priority_category'] == 'High'

output = ["Feature\tOR\tLower_CI\tUpper_CI\tp-value"]
for display_name, col_name in features.items():
    if col_name not in df.columns:
        continue
    
    A = sum((high_mask) & (df[col_name] == True))
    B = sum((high_mask) & (df[col_name] == False))
    C = sum((~high_mask) & (df[col_name] == True))
    D = sum((~high_mask) & (df[col_name] == False))
    
    table = [[A, B],
             [C, D]]
    
    odds, p_value = fisher_exact(table)
    
    # Calculate Exact confidence interval (scipy 1.7+)
    res = odds_ratio(table)
    ci = res.confidence_interval(confidence_level=0.95)
    
    output.append(f"{display_name}\t{odds:.4f}\t{ci.low:.4f}\t{ci.high:.4f}\t{p_value:.4e}")

out_text = "\n".join(output)
print(out_text)
with open('data/processed/chrna7_structural_ci.txt', 'w', encoding='utf-8') as f:
    f.write(out_text)
