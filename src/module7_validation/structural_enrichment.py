"""
Module 7: Structural Feature Enrichment Analysis (RAREMISS)
Simple Fisher's Exact Test — see structural_enrichment_ci.py for full output.

IMPORTANT METHODOLOGICAL NOTE
------------------------------
This analysis tests whether the high-priority RAREMISS category is enriched
for structural features (binding site, pore, interface) that are *also* used
as scoring inputs.  By construction this creates a circularity: the score
rewards structural features, so high-scoring variants are more likely to have
them.

This analysis should therefore be interpreted as an internal score-consistency
or descriptive check — confirming that the scoring formula behaves as intended —
NOT as independent biological validation of pathogenicity or evidence that
RAREMISS outperforms other predictors.

Contingency-table conventions
-------------------------------
'unknown' values (spatially unresolved or no-ligand) are EXCLUDED from
contingency cells.  Only True (assessed positive) and False (assessed
negative) contribute to the 2×2 table.  The excluded count is reported
separately.
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import os
import sys
import yaml

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from src.raremiss_utils import load_spatial_csv

config_path = os.path.join(root_dir, "config", "parameters.yaml")
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)
gene_symbol = config.get("gene_symbol", "CHRNA7").lower()

input_csv = os.path.join(root_dir, f'data/processed/{gene_symbol}_ranked_variants.csv')
df = load_spatial_csv(input_csv)

features = {
    'is_binding_site': 'is_binding_site',
    'is_pore_core':    'is_pore_region',  # column name in our dataset
    'is_interface':    'is_interface',
}

high_mask     = df['priority_category'] == 'High'
not_high_mask = df['priority_category'] != 'High'

output = [
    "[NOTE] RAREMISS Internal Score-Consistency Analysis (Descriptive check, not independent pathogenicity validation)\n"
]
for display_name, col_name in features.items():
    if col_name not in df.columns:
        output.append(f"Column {col_name} not found.")
        continue

    # Exclude 'unknown' — only True/False participate in the contingency table
    col_true  = df[col_name] == True   # noqa: E712
    col_false = df[col_name] == False  # noqa: E712
    col_unk   = df[col_name] == 'unknown'

    A = int((high_mask & col_true).sum())
    B = int((high_mask & col_false).sum())
    C = int((not_high_mask & col_true).sum())
    D = int((not_high_mask & col_false).sum())
    excluded_unknown = int(col_unk.sum())

    table = [[A, B], [C, D]]

    if A + B + C + D == 0:
        output.append(f"--- {display_name} ---")
        output.append(f"Skipped: all {excluded_unknown} values are 'unknown' (no assessed values)")
        continue

    odds, p_value = fisher_exact(table)

    output.append(f"--- {display_name} ---")
    output.append(f"A (High, True)       : {A}")
    output.append(f"B (High, False)      : {B}")
    output.append(f"C (Not High, True)   : {C}")
    output.append(f"D (Not High, False)  : {D}")
    output.append(f"Excluded (unknown)   : {excluded_unknown}")
    output.append(f"Odds Ratio           : {odds:.4f}")
    output.append(f"p-value              : {p_value:.4e}\n")

print("\n".join(output))
