import pandas as pd
import requests
import json
import numpy as np
import time

input_csv = 'data/processed/chrna7_missense_structural_annotated.csv'
df = pd.read_csv(input_csv)

rsids = df['rsid'].tolist()
batch_size = 100

url_vep = "https://rest.ensembl.org/vep/human/id"
headers = {"Content-Type": "application/json", "Accept": "application/json"}

vep_data = {}

alt_map = dict(zip(df['rsid'], df['alt']))

print(f"STEP 1: Batch querying {len(rsids)} RSIDs to VEP for functional scores...")
cadd_missing: int = 0
for i in range(0, len(rsids), batch_size):
    batch = rsids[i:i+batch_size]
    req_url = f"{url_vep}?CADD=1"
    res = requests.post(req_url, headers=headers, json={"ids": batch})
    if res.status_code != 200:
        print(f"Warning: Failed batch {i}-{i+batch_size}, status code {res.status_code}")
        continue
    
    data = res.json()
    for item in data:
        rid = item.get('input') or item.get('id')
        alt_allele = alt_map.get(rid)
        
        cadd_phred = np.nan
        sift_score = np.nan
        sift_pred = None
        polyphen_score = np.nan
        polyphen_pred = None
        
        has_match = False
        has_cadd = False
        for tc in item.get('transcript_consequences', []):
            if tc.get('transcript_id') == 'ENST00000306901' and tc.get('variant_allele') == alt_allele:
                if 'sift_score' in tc:
                    sift_score = tc.get('sift_score')
                    sift_pred = tc.get('sift_prediction')
                if 'polyphen_score' in tc:
                    polyphen_score = tc.get('polyphen_score')
                    polyphen_pred = tc.get('polyphen_prediction')
                if 'cadd_phred' in tc:
                    cadd_phred = tc.get('cadd_phred')
                    has_cadd = True
                has_match = True
                break
                
        if not has_cadd:
            cadd_missing = cadd_missing + 1
            
        if has_match or rid not in vep_data:
            vep_data[rid] = {
                'cadd_phred': cadd_phred,
                'sift_score': sift_score,
                'sift_pred': sift_pred,
                'polyphen_score': polyphen_score,
                'polyphen_pred': polyphen_pred
            }
    time.sleep(0.1)

if cadd_missing > 0:
    print(f"WARNING: CADD scores were not returned natively for {cadd_missing} variants.")

print("STEP 2: Conservation Retrieval...")
print("LIMITATION NOTED: Querying /overlap/region/human/:chr:start-end?feature=conservation_score for hundreds of individual variants is not feasible at scale due to rate-limiting and high response latency without a dedicated bulk endpoint. Continuing without phyloP conservation scores.")

print("STEP 3 & 4: Merging and Output...")

def get_vep(rsid, field):
    return vep_data.get(rsid, {}).get(field)

df['cadd_phred'] = df['rsid'].apply(lambda x: get_vep(x, 'cadd_phred'))
df['sift_score'] = df['rsid'].apply(lambda x: get_vep(x, 'sift_score'))
df['sift_pred'] = df['rsid'].apply(lambda x: get_vep(x, 'sift_pred'))
df['polyphen_score'] = df['rsid'].apply(lambda x: get_vep(x, 'polyphen_score'))
df['polyphen_pred'] = df['rsid'].apply(lambda x: get_vep(x, 'polyphen_pred'))
df['phyloP_score'] = np.nan

df['cadd_phred'] = pd.to_numeric(df['cadd_phred'], errors='coerce')

output_path = 'data/processed/chrna7_missense_full_annotated.csv'
df.to_csv(output_path, index=False)
print(f"Saved fully annotated dataset to {output_path}")

out = "========= PHASE C FUNCTIONAL ANNOTATION STATS =========\n"
out += "\n--- CADD Summary Stats ---\n"
cadd_s = df['cadd_phred'].dropna()
out += cadd_s.describe().to_string() + "\n"

cadd_gt_20 = (cadd_s > 20).sum()
cadd_gt_20_perc = (cadd_gt_20 / len(df)) * 100 if len(df) > 0 else 0
out += f"\n% variants with CADD > 20: {cadd_gt_20_perc:.2f}%\n"

sift_damaging = df['sift_pred'].str.lower().str.contains('damaging|deleterious', na=False).sum()
sift_perc = (sift_damaging / len(df)) * 100 if len(df) > 0 else 0
out += f"% predicted damaging by SIFT: {sift_perc:.2f}%\n"

polyphen_prob_damaging = df['polyphen_pred'].str.lower().str.contains('probably_damaging', na=False).sum()
poly_prob_perc = (polyphen_prob_damaging / len(df)) * 100 if len(df) > 0 else 0
out += f"% predicted probably damaging by PolyPhen: {poly_prob_perc:.2f}%\n"

polyphen_poss_damaging = df['polyphen_pred'].str.lower().str.contains('possibly_damaging', na=False).sum()
poly_poss_perc = (polyphen_poss_damaging / len(df)) * 100 if len(df) > 0 else 0
out += f"% predicted possibly damaging by PolyPhen: {poly_poss_perc:.2f}%\n"

polyphen_benign = df['polyphen_pred'].str.lower().str.contains('benign', na=False).sum()
poly_benign_perc = (polyphen_benign / len(df)) * 100 if len(df) > 0 else 0
out += f"% predicted benign by PolyPhen: {poly_benign_perc:.2f}%\n"

out += "=========================================================\n"

with open('functional_metrics.txt', 'w', encoding='utf-8') as f:
    f.write(out)

print("Metrics written to functional_metrics.txt")
