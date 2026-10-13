"""
Module 7: Structural Feature Enrichment Analysis with CI (RAREMISS)

IMPORTANT METHODOLOGICAL NOTE — INTERNAL CONSISTENCY ANALYSIS
--------------------------------------------------------------
This analysis tests whether the high-priority RAREMISS category is enriched
for structural features that are *also* used as scoring inputs.
By construction this creates a circularity: the score rewards structural
features, so high-scoring variants tend to have them.

Results should therefore be interpreted as an internal score-consistency
check, confirming the scoring formula behaves as designed.  They do NOT
constitute independent biological validation of pathogenicity, nor evidence
that RAREMISS outperforms sequence-based predictors.

Tested hypotheses and multiple-testing correction
--------------------------------------------------
Three features are tested per gene: is_pore_core, is_interface, is_binding_site.
Raw p-values are corrected using the Benjamini–Hochberg (BH) FDR procedure
applied across the set of successfully computed tests.
Skipped tests are NOT included in the correction family; this is noted in output.
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
import os
import sys
import yaml
import json

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from src.raremiss_utils import load_spatial_csv, get_af_threshold, benjamini_hochberg

config_path = os.path.join(root_dir, "config", "parameters.yaml")
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)
gene_symbol  = config.get("gene_symbol", "CHRNA7").lower()
af_threshold = get_af_threshold(config, source="config/parameters.yaml")

input_csv = os.path.join(root_dir, f'data/processed/{gene_symbol}_ranked_variants.csv')
# load_spatial_csv normalises True/False/'unknown' after CSV round-trip
df = load_spatial_csv(input_csv)

features = {
    'is_pore_core': 'is_pore_region',
    'is_interface': 'is_interface',
}

ANALYSIS_TYPE = "internal_score_consistency"
METHODOLOGY_NOTE = "Descriptive internal consistency check; not independent pathogenicity validation"

output = []

# Handle binding site separately: skip only when there are NO assessed values
if 'is_binding_site' in df.columns:
    bs_col = df['is_binding_site']
    bs_has_assessed = bs_col.isin([True, False]).any()
    bs_true_count   = (bs_col == True).sum()   # noqa: E712

    if not bs_has_assessed:
        msg = "Binding-site enrichment skipped (no ligand — all values unresolved)"
        print(msg)
        output.append({
            "feature": "is_binding_site",
            "analysis_type": ANALYSIS_TYPE,
            "methodology_note": METHODOLOGY_NOTE,
            "status": "skipped",
            "reason": msg
        })
    elif bs_true_count == 0:
        msg = "Binding-site enrichment skipped (ligand present but no variants within 5 Å)"
        print(msg)
        output.append({
            "feature": "is_binding_site",
            "analysis_type": ANALYSIS_TYPE,
            "methodology_note": METHODOLOGY_NOTE,
            "status": "skipped",
            "reason": msg
        })
    else:
        features['is_binding_site'] = 'is_binding_site'

high_mask     = df['priority_category'] == 'High'
not_high_mask = df['priority_category'] != 'High'

computed_entries = []  # collect (index_in_output, p_value) for BH correction

for display_name, col_name in features.items():
    if col_name not in df.columns:
        continue

    # Only True and False contribute to the contingency table.
    # 'unknown' values (spatially unresolved, no-ligand) are excluded.
    col_true  = df[col_name] == True   # noqa: E712
    col_false = df[col_name] == False  # noqa: E712
    col_unk   = df[col_name] == 'unknown'

    A = int((high_mask & col_true).sum())
    B = int((high_mask & col_false).sum())
    C = int((not_high_mask & col_true).sum())
    D = int((not_high_mask & col_false).sum())
    excluded_unknown = int(col_unk.sum())

    if A + B + C + D == 0:
        reason = f"Analysis population empty — all {excluded_unknown} values are 'unknown'"
        print(f"Skipping {display_name}: {reason}")
        output.append({
            "feature": display_name,
            "analysis_type": ANALYSIS_TYPE,
            "methodology_note": METHODOLOGY_NOTE,
            "status": "skipped",
            "reason": reason,
            "excluded_unknown": excluded_unknown
        })
        continue

    if A < 3 or B < 3 or C < 3 or D < 3:
        reason = (f"Insufficient data for enrichment analysis "
                  f"(cell count < 3; A={A},B={B},C={C},D={D})")
        print(f"Skipping {display_name}: {reason}")
        output.append({
            "feature": display_name,
            "analysis_type": ANALYSIS_TYPE,
            "methodology_note": METHODOLOGY_NOTE,
            "status": "skipped",
            "reason": reason,
            "excluded_unknown": excluded_unknown,
            "contingency": {"A": A, "B": B, "C": C, "D": D}
        })
        continue

    table = [[A, B], [C, D]]
    odds, p_value = fisher_exact(table)

    try:
        res_ci = odds_ratio(table)
        ci = res_ci.confidence_interval(confidence_level=0.95)
        if (np.isinf(odds) or np.isnan(odds) or
                np.isinf(ci.low) or np.isinf(ci.high)):
            msg = "OR or CI is infinite or NaN (zero/sparse cell)"
            print(f"Skipping {display_name}: {msg}")
            output.append({
                "feature": display_name,
                "analysis_type": ANALYSIS_TYPE,
                "methodology_note": METHODOLOGY_NOTE,
                "status": "skipped",
                "reason": msg,
                "excluded_unknown": excluded_unknown,
                "contingency": {"A": A, "B": B, "C": C, "D": D}
            })
            continue

        entry = {
            "feature":          display_name,
            "analysis_type":    ANALYSIS_TYPE,
            "methodology_note": METHODOLOGY_NOTE,
            "status":           "success",
            "odds_ratio":       float(odds),
            "lower_ci":         float(ci.low),
            "upper_ci":         float(ci.high),
            "p_value":          float(p_value),
            "p_value_bh_adj":   None,    # filled in after BH correction
            "excluded_unknown": excluded_unknown,
            "contingency":      {"A": A, "B": B, "C": C, "D": D},
        }
        computed_entries.append(len(output))
        output.append(entry)

    except Exception as e:
        print(f"Error {display_name}: {e}")
        output.append({
            "feature": display_name,
            "analysis_type": ANALYSIS_TYPE,
            "methodology_note": METHODOLOGY_NOTE,
            "status": "error",
            "reason": str(e),
            "excluded_unknown": excluded_unknown
        })


# --- Benjamini–Hochberg FDR correction over successfully computed tests ---
if computed_entries:
    raw_ps = [output[idx]["p_value"] for idx in computed_entries]
    adj_ps = benjamini_hochberg(raw_ps)
    for idx, adj_val in zip(computed_entries, adj_ps):
        output[idx]["p_value_bh_adj"] = float(adj_val)


output_path = os.path.join(root_dir, f'data/processed/{gene_symbol}_enrichment_results.json')
with open(output_path, 'w', encoding='utf-8') as f:
    json.dump(output, f, indent=4)
print(json.dumps(output, indent=4))
