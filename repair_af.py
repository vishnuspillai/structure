import requests
import json
import pandas as pd
import numpy as np

print("========== STEP 1: Variation Response Structure ==========")
rsid = 'rs2052300715'
url = 'https://rest.ensembl.org/variation/homo_sapiens?pops=1'
headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
res = requests.post(url, headers=headers, json={"ids": [rsid]})
data = res.json()
pops = data[rsid].get('populations', [])

print(f"response['{rsid}']['populations'] (First 3 entries):")
print(json.dumps(pops[:3], indent=2))
keys = sorted(list(set([p.get('population') for p in pops])))
print("\nAvailable population keys:")
for k in keys:
    print(k)

print("\n========== STEP 2 & 3: Extract True gnomAD AF & Recompute Column ==========")
df = pd.read_csv('data/processed/chrna7_rare_missense_variants.csv')
rsids = df['rsid'].tolist()
alt_map = dict(zip(df['rsid'], df['alt']))
batch_size = 100
af_map = {}

print(f"Batch querying {len(rsids)} variants...")
for i in range(0, len(rsids), batch_size):
    batch = rsids[i:i+batch_size]
    r = requests.post(url, headers=headers, json={"ids": batch})
    if r.status_code != 200:
        continue
    d = r.json()
    for rid, info in d.items():
        if not info:
            continue
        max_af = np.nan
        alt_allele = alt_map.get(rid)
        for pop in info.get('populations', []):
            p_name = pop.get('population', '').lower()
            if ('gnomadg' in p_name or 'gnomade' in p_name) and pop.get('allele') == alt_allele:
                freq = pop.get('frequency', pop.get('allele_frequency'))
                if freq is not None:
                    f_val = float(freq)
                    if pd.isna(max_af) or f_val > max_af:
                        max_af = f_val
        af_map[rid] = max_af

df['gnomAD_AF'] = df['rsid'].map(af_map)

nan_count = df['gnomAD_AF'].isna().sum()
zero_count = (df['gnomAD_AF'] == 0.0).sum()
af_min = df['gnomAD_AF'].min()
af_max = df['gnomAD_AF'].max()
af_mean = df['gnomAD_AF'].mean()

print(f"Count NaN: {nan_count}")
print(f"Count AF == 0.0: {zero_count}")
print(f"Min AF: {af_min}")
print(f"Max AF: {af_max}")
print(f"Mean AF: {af_mean}")

df.to_csv('data/processed/chrna7_missense_master_corrected.csv', index=False)
