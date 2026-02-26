import pandas as pd
import requests
import xml.etree.ElementTree as ET
from Bio import Entrez
import time
import urllib.parse
import os

Entrez.email = "vishnu@example.com"  # Please configure

df = pd.read_csv('data/processed/chrna7_top10_structural_metrics.csv')

def query_clinvar(rsid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={rsid}&retmode=json"
    r = requests.get(url)
    if r.status_code == 200:
        data = r.json()
        ids = data.get("esearchresult", {}).get("idlist", [])
        if ids:
            fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={','.join(ids)}&retmode=json"
            f_r = requests.get(fetch_url)
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
    return "Not found", "Not found"

def query_pubmed(aa_change):
    query = f"CHRNA7 AND {aa_change}"
    encoded = urllib.parse.quote(query)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={encoded}&retmode=json"
    r = requests.get(url)
    if r.status_code == 200:
        return int(r.json().get('esearchresult', {}).get('count', 0))
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

output_path = 'data/processed/chrna7_top10_clinical_lit.csv'
res_df.to_csv(output_path, index=False)

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', None)
print(res_df.to_string(index=False))
