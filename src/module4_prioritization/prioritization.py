"""
Module 4: Variant Prioritization and Scoring (RAREMISS)

Scoring formula and categories are defined in src/raremiss_utils.compute_score().
This module applies the score, assigns categories, and writes the ranked output.

Score semantics: heuristic prioritization, not calibrated pathogenicity probability.
See src/raremiss_utils.py for the full documented specification.
"""

import pandas as pd
import numpy as np
import os
import sys
import yaml

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from src.raremiss_utils import (
    compute_score, get_category, load_spatial_csv,
    get_af_threshold, DEFAULT_AF_THRESHOLD
)

config_path = os.path.join(root_dir, "config", "parameters.yaml")
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)

gene_symbol  = config.get("gene_symbol", "CHRNA7").lower()
af_threshold = get_af_threshold(config, source="config/parameters.yaml")

input_csv = os.path.join(root_dir, f"data/processed/{gene_symbol}_missense_spatial_annotated.csv")
if not os.path.exists(input_csv):
    raise FileNotFoundError(
        f"Expected input file not found: {input_csv}. Check Module 3 output."
    )

# Use load_spatial_csv to normalise True/False/'unknown' after CSV round-trip.
df = load_spatial_csv(input_csv)

df['priority_score']    = df.apply(compute_score, axis=1)
df['priority_category'] = df['priority_score'].apply(get_category)
df = df.sort_values(by='priority_score', ascending=False)

output_path = os.path.join(root_dir, f"data/processed/{gene_symbol}_ranked_variants.csv")
df.to_csv(output_path, index=False)

# Gene-namespaced text summary (replaces non-namespaced mod4_output.txt)
output_txt = os.path.join(root_dir, f"data/processed/{gene_symbol}_priority_summary.txt")
with open(output_txt, 'w', encoding='utf-8') as f:
    f.write(f"RAREMISS Priority Summary — {gene_symbol.upper()}\n")
    f.write(f"AF threshold used: {af_threshold}\n\n")
    f.write("Score distribution:\n")
    f.write(df['priority_score'].describe().to_string() + "\n")
    f.write("\nCount per category:\n")
    f.write(df['priority_category'].value_counts().to_string() + "\n")
    f.write("\nTop 15 variants:\n")
    f.write(
        df[['rsid', 'protein_position', 'domain_region', 'priority_score']]
        .head(15).to_string(index=False) + "\n"
    )

print(f"Saved: {output_path}")
print(df['priority_category'].value_counts().to_string())
