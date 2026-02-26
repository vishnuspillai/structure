import pandas as pd
import requests
import yaml
import os

def fetch_uniprot_features(uniprot_id="P36544"):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    res = requests.get(url)
    res.raise_for_status()
    data = res.json()
    
    domains = {}
    transmem_count = 1
    
    for f in data.get('features', []):
        t = f['type']
        if 'location' not in f or 'start' not in f['location'] or 'end' not in f['location']:
            continue
            
        start = f['location']['start']['value']
        end = f['location']['end']['value']
        desc = f.get('description', '').lower()
        
        if t == 'Signal':
            domains['signal_peptide'] = [start, end]
        elif t == 'Topological domain':
            if 'extracellular' in desc:
                if 'extracellular_domain' not in domains:
                    domains['extracellular_domain'] = [start, end]
            elif 'cytoplasmic' in desc:
                domains['intracellular_loop'] = [start, end]
        elif t == 'Transmembrane':
            domains[f'M{transmem_count}'] = [start, end]
            transmem_count += 1
            
    return domains

print("STEP 1: Fetching UniProt features...")
domains = fetch_uniprot_features()

config_path = 'config/parameters.yaml'
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

config['domains'] = domains
with open(config_path, 'w') as file:
    yaml.safe_dump(config, file)
print("Domains saved to config/parameters.yaml")

print("STEP 2 & 3: Assigning structural domains...")
df = pd.read_csv('data/processed/chrna7_missense_rare_corrected.csv')

def assign_domain(pos):
    if pd.isna(pos):
        return "other"
    pos = int(pos)
    for d_name, d_bounds in domains.items():
        if d_bounds[0] <= pos <= d_bounds[1]:
            return d_name
    return "other"

df['domain_region'] = df['protein_position'].apply(assign_domain)
df['is_transmembrane'] = df['domain_region'].isin(['M1', 'M2', 'M3', 'M4'])
df['is_pore_region'] = df['domain_region'] == 'M2'
df['is_extracellular'] = df['domain_region'] == 'extracellular_domain'

annotated_path = 'data/processed/chrna7_missense_structural_annotated.csv'
df.to_csv(annotated_path, index=False)
print("STEP 4: Output")
print(f"Saved annotated dataset to {annotated_path}")

out = "========= PHASE B STRUCTURAL ANNOTATION STATS =========\n"
out += "\n--- Domain Distribution ---\n"
out += df['domain_region'].value_counts(dropna=False).to_string() + "\n"

tm_perc = (df['is_transmembrane'].sum() / len(df)) * 100
pore_perc = (df['is_pore_region'].sum() / len(df)) * 100
ec_perc = (df['is_extracellular'].sum() / len(df)) * 100

out += f"\n% in transmembrane: {tm_perc:.2f}%\n"
out += f"% in M2 (pore): {pore_perc:.2f}%\n"
out += f"% extracellular: {ec_perc:.2f}%\n"
out += "=========================================================\n"

with open('structural_metrics.txt', 'w', encoding='utf-8') as f:
    f.write(out)

print("Metrics written to structural_metrics.txt")
