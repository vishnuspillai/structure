import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
import os
import yaml
import numpy as np
import json

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
config_path = os.path.join(root_dir, "config", "parameters.yaml")
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)
gene_symbol = config.get("gene_symbol", "CHRNA7").lower()

input_csv = os.path.join(root_dir, f'data/processed/{gene_symbol}_ranked_variants.csv')
df = pd.read_csv(input_csv)

features = {
    'is_pore_core': 'is_pore_region',
    'is_interface': 'is_interface'
}

output = []

# Check if binding site analysis was performed
if 'is_binding_site' in df.columns:
    if 'unknown' in df['is_binding_site'].unique():
        msg = "Binding-site enrichment not computed (no ligand available)"
        print(msg)
        output.append({"feature": "is_binding_site", "status": "skipped", "reason": msg})
    else:
        features['is_binding_site'] = 'is_binding_site'

high_mask = df['priority_category'] == 'High'
not_high_mask = df['priority_category'] != 'High'

for display_name, col_name in features.items():
    if col_name not in df.columns:
        continue
    
    A = len(df[high_mask & (df[col_name] == True)])
    B = len(df[high_mask & (df[col_name] == False)])
    C = len(df[not_high_mask & (df[col_name] == True)])
    D = len(df[not_high_mask & (df[col_name] == False)])
    
    if A < 3 or B < 3 or C < 3 or D < 3:
        print(f"Skipping {display_name}: Insufficient data for enrichment analysis (cell count < 3).")
        output.append({"feature": display_name, "status": "skipped", "reason": "Insufficient data (cell count < 3)"})
        continue

    table = [[A, B],
             [C, D]]
    
    odds, p_value = fisher_exact(table)
    
    # Calculate Exact confidence interval (scipy 1.7+)
    try:
        res = odds_ratio(table)
        ci = res.confidence_interval(confidence_level=0.95)
        if np.isinf(odds) or np.isnan(odds) or np.isinf(ci.low) or np.isinf(ci.high):
            print(f"Skipping {display_name}: OR is infinity or NaN.")
            output.append({"feature": display_name, "status": "skipped", "reason": "OR is infinity or NaN"})
            continue
        output.append({
            "feature": display_name,
            "status": "success",
            "odds_ratio": float(odds),
            "lower_ci": float(ci.low),
            "upper_ci": float(ci.high),
            "p_value": float(p_value)
        })
    except Exception as e:
        print(f"Skipping {display_name}: Stats error - {e}")
        output.append({"feature": display_name, "status": "error", "reason": str(e)})

# Final JSON output
output_path = os.path.join(root_dir, f'data/processed/{gene_symbol}_enrichment_results.json')
with open(output_path, 'w', encoding='utf-8') as f:
    json.dump(output, f, indent=4)

# Print for log
print(json.dumps(output, indent=4))
