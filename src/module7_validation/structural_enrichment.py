import pandas as pd
from scipy.stats import fisher_exact

df = pd.read_csv('data/processed/chrna7_ranked_variants.csv')

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
