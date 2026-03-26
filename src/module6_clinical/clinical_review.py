import pandas as pd
import requests
from Bio import Entrez
import time
import urllib.parse
import os
import yaml

Entrez.email = "vishnu@example.com"  # Please configure

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
config_path = os.path.join(PROJECT_ROOT, "config", "parameters.yaml")
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)
gene_symbol = config.get("gene_symbol", "CHRNA7").lower()

input_csv = os.path.join(PROJECT_ROOT, 'data', 'processed', f'{gene_symbol}_top10_structural_metrics.csv')
df = pd.read_csv(input_csv)

def query_clinvar(rsid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={rsid}&retmode=json"
    try:
        r = requests.get(url, timeout=15)
        if r.status_code == 200:
            data = r.json()
        ids = data.get("esearchresult", {}).get("idlist", [])
        if ids:
            fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={','.join(ids)}&retmode=json"
            f_r = requests.get(fetch_url, timeout=15)
            if f_r.status_code == 200:
                f_data = f_r.json()
                result = f_data.get('result', {})
                for uid in ids:
                    record = result.get(uid, {})
                    cg = record.get('clinical_significance', {})
                    sig = cg.get('description', 'Not Provided')
                    
                    traits = record.get('trait_set', [])
                    condition = "Not Provided"
                    if traits:
                        condition = traits[0].get('trait_name', 'Not Provided')
                        
                    return sig, condition
    except Exception as e:
        print(f"ClinVar query failed for {rsid}: {e}")
    return "Not found", "Not found"

def query_pubmed(aa_change):
    query = f"{config.get('gene_symbol', 'CHRNA7')} AND {aa_change}"
    encoded = urllib.parse.quote(query)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={encoded}&retmode=json"
    try:
        r = requests.get(url, timeout=15)
        if r.status_code == 200:
            return int(r.json().get('esearchresult', {}).get('count', 0))
    except Exception as e:
        print(f"PubMed query failed for {aa_change}: {e}")
    return 0

results = []
for index, row in df.iterrows():
    rsid = row['rsid']
    aa_change = row['amino_acid_change']
    
    cv_sig, cv_cond = query_clinvar(rsid)
    pm_hits = query_pubmed(aa_change)
    
    results.append({
        'rsid': rsid,
        'aa_change': aa_change,
        'ClinVar_significance': cv_sig,
        'ClinVar_condition': cv_cond,
        'PubMed_hits': pm_hits
    })
    time.sleep(0.4) # Entrez rate limit

res_df = pd.DataFrame(results)

output_path = os.path.join(PROJECT_ROOT, 'data', 'processed', f'{gene_symbol}_top10_clinical_lit.csv')
res_df.to_csv(output_path, index=False)

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', None)
print(res_df.to_string(index=False))
