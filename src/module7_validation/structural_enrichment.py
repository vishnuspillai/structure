import pandas as pd
from scipy.stats import fisher_exact
import os
import yaml

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
config_path = os.path.join(root_dir, "config", "parameters.yaml")
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)
gene_symbol = config.get("gene_symbol", "CHRNA7").lower()

input_csv = os.path.join(root_dir, f'data/processed/{gene_symbol}_ranked_variants.csv')
df = pd.read_csv(input_csv)

features = {
    'is_binding_site': 'is_binding_site',
    'is_pore_core': 'is_pore_region', # the column in our dataset
    'is_interface': 'is_interface'
}

high_mask = df['priority_category'] == 'High'
not_high_mask = df['priority_category'] != 'High'

output = []
for display_name, col_name in features.items():
    if col_name not in df.columns:
        output.append(f"Column {col_name} not found.")
        continue
    
    #     | Feature=True | Feature=False
    # High|      A       |      B
    # Not |      C       |      D
    
    A = len(df[high_mask & (df[col_name] == True)])
    B = len(df[high_mask & (df[col_name] == False)])
    C = len(df[not_high_mask & (df[col_name] == True)])
    D = len(df[not_high_mask & (df[col_name] == False)])
    
    table = [[A, B],
             [C, D]]
    
    odds, p_value = fisher_exact(table)
    
    output.append(f"--- {display_name} ---")
    output.append(f"A (High, True)      : {A}")
    output.append(f"B (High, False)     : {B}")
    output.append(f"C (Not High, True)  : {C}")
    output.append(f"D (Not High, False) : {D}")
    output.append(f"Odds Ratio          : {odds}")
    output.append(f"p-value             : {p_value}\n")

print("\n".join(output))
