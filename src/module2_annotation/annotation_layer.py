import os
import sys
import logging
import yaml
import requests
import json
import pandas as pd
import numpy as np
from typing import Dict, List, Any

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def setup_logging(log_file: str) -> logging.Logger:
    logger = logging.getLogger('module2_annotation')
    logger.setLevel(logging.INFO)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

def phase_a_coordinate_correction(df: pd.DataFrame, logger: logging.Logger, output_path: str) -> pd.DataFrame:
    """Phase A: Fetch GRCh38 coordinates via /variation API, update chrom, pos and audit AF."""
    url = "https://rest.ensembl.org/variation/homo_sapiens?pops=1"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    rsids = df['rsid'].tolist()
    batch_size = 100
    coord_map = {}
    af_map = {}
    alt_map = dict(zip(df['rsid'], df['alt']))
    
    logger.info(f"Phase A: Batch querying {len(rsids)} RSIDs for coordinate correction and AF repair.")
    for i in range(0, len(rsids), batch_size):
        batch = rsids[i:i+batch_size]
        payload = {"ids": batch}
        res = requests.post(url, headers=headers, json=payload)
        if res.status_code != 200:
            logger.warning(f"Failed to fetch batch {i}-{i+batch_size}: {res.status_code}")
            continue
            
        data = res.json()
        for rsid, info in data.items():
            if not info:
                continue
            mappings = info.get('mappings', [])
            for m in mappings:
                if m.get('assembly_name') == 'GRCh38':
                    coord_map[rsid] = {
                        'chrom': m.get('seq_region_name'),
                        'pos': m.get('start'),
                        'ref': m.get('allele_string', '').split('/')[0] if '/' in m.get('allele_string', '') else '',
                        'alt': m.get('allele_string', '').split('/')[-1] if '/' in m.get('allele_string', '') else ''
                    }
                    break
                    
            max_af = np.nan
            alt_allele = alt_map.get(rsid)
            for pop in info.get('populations', []):
                p_name = pop.get('population', '').lower()
                if ('gnomadg' in p_name or 'gnomade' in p_name) and pop.get('allele') == alt_allele:
                    freq = pop.get('frequency', pop.get('allele_frequency'))
                    if freq is not None:
                        f_val = float(freq)
                        if pd.isna(max_af) or f_val > max_af:
                            max_af = f_val
            af_map[rsid] = max_af
                    
    logger.info(f"Phase A: Successfully mapped coordinates and AF for {len(coord_map)} RSIDs.")
    
    # Map back to DataFrame
    df = df.copy()
    def get_coord(rsid, field):
        return coord_map.get(rsid, {}).get(field)
    
    df['chrom'] = df['rsid'].apply(lambda x: get_coord(x, 'chrom') or df.loc[df['rsid'] == x, 'chrom'].values[0])
    df['pos'] = df['rsid'].apply(lambda x: get_coord(x, 'pos') or df.loc[df['rsid'] == x, 'pos'].values[0])
    df['gnomAD_AF'] = df['rsid'].map(af_map)
    
    # Ensure pos is numeric
    df['pos'] = pd.to_numeric(df['pos'], errors='coerce')
    
    # Audit AF
    logger.info("Phase A: Auditing gnomAD_AF distribution...")
    zero_count = (df['gnomAD_AF'] == 0.0).sum()
    null_count = df['gnomAD_AF'].isnull().sum()
    af_min = df['gnomAD_AF'].min()
    af_max = df['gnomAD_AF'].max()
    af_mean = df['gnomAD_AF'].mean()
    
    logger.info(f"AF == 0.0: {zero_count} ({zero_count/len(df)*100:.2f}%)")
    logger.info(f"AF is null: {null_count}")
    logger.info(f"AF min: {af_min}, max: {af_max}, mean: {af_mean}")
    
    if (zero_count / len(df)) > 0.80:
        logger.warning("WARNING: AF unreliable (Extremely High percentage of 0.0 AF values reported). Stop sequence engaged.")
        print("WARNING: AF unreliable. Stopping execution.")
        sys.exit(1)
        
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    logger.info(f"Phase A: Corrected CSV saved to {output_path}")
    return df

def fetch_uniprot_features(uniprot_id: str, logger: logging.Logger) -> Dict[str, Any]:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    logger.info(f"Phase B: Fetching canonical UniProt record for {uniprot_id}")
    res = requests.get(url)
    res.raise_for_status()
    data = res.json()
    
    domains = {}
    transmem_count = 1
    
    for f in data.get('features', []):
        t = f['type']
        start = f['location']['start']['value']
        end = f['location']['end']['value']
        desc = f.get('description', '').lower()
        
        if t == 'Signal':
            domains['signal_peptide'] = [start, end]
        elif t == 'Topological domain':
            if 'extracellular' in desc:
                # Capture the largest ECD, usually before M1
                if 'extracellular_domain' not in domains:
                    domains['extracellular_domain'] = [start, end]
            elif 'cytoplasmic' in desc:
                domains['intracellular_loop'] = [start, end]
        elif t == 'Transmembrane':
            domains[f'm{transmem_count}'] = [start, end]
            transmem_count += 1
            
    logger.info(f"Phase B: Parsed domain boundaries: {domains}")
    return domains

def phase_b_structural_domains(df: pd.DataFrame, logger: logging.Logger, config_path: str) -> pd.DataFrame:
    """Phase B: Parse UniProt limits, alter params, assign domains."""
    uniprot_id = "P36544" # Canonical CHRNA7
    domains = fetch_uniprot_features(uniprot_id, logger)
    
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
        
    config['domains'] = domains
    with open(config_path, 'w') as file:
        yaml.safe_dump(config, file)
    logger.info(f"Phase B: Updated {config_path} with Domain bounds.")
    
    def assign_domain(pos):
        if pd.isna(pos):
            return "Unknown"
        pos = int(pos)
        for d_name, d_bounds in domains.items():
            if d_bounds[0] <= pos <= d_bounds[1]:
                return d_name
        return "Unassigned"
        
    df = df.copy()
    df['domain_region'] = df['protein_position'].apply(assign_domain)
    df['is_transmembrane'] = df['domain_region'].isin(['m1', 'm2', 'm3', 'm4'])
    df['is_pore_region'] = df['domain_region'] == 'm2'
    df['is_extracellular'] = df['domain_region'] == 'extracellular_domain'
    
    return df

def phase_c_functional_scores(df: pd.DataFrame, logger: logging.Logger, output_path: str) -> pd.DataFrame:
    """Phase C: VEP querying for CADD and Conservation."""
    url = "https://rest.ensembl.org/vep/human/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    params = {"CADD": 1}
    
    rsids = df['rsid'].tolist()
    batch_size = 100
    vep_map = {}
    
    logger.info(f"Phase C: Batch querying {len(rsids)} RSIDs to VEP for CADD.")
    for i in range(0, len(rsids), batch_size):
        batch = rsids[i:i+batch_size]
        payload = {"ids": batch}
        # post query string via url payload? No, VEP API accepts URL params in the URL itself
        req_url = f"{url}?CADD=1"
        res = requests.post(req_url, headers=headers, json=payload)
        
        if res.status_code != 200:
            logger.warning(f"Failed to fetch batch {i}-{i+batch_size}: {res.status_code}")
            continue
            
        data = res.json()
        for item in data:
            rsid = item.get('id')
            cadd_phred = None
            sift = None
            polyphen = None
            
            # VEP returns multiple transcript consequences. Try to grab the max or one matching canonical.
            for tc in item.get('transcript_consequences', []):
                # SIFT/PolyPhen
                if 'sift_score' in tc and sift is None:
                    sift = tc.get('sift_score')
                if 'polyphen_score' in tc and polyphen is None:
                    polyphen = tc.get('polyphen_score')
                    
            if 'transcript_consequences' in item and len(item['transcript_consequences']) > 0:
                tc = item['transcript_consequences'][0]
                if 'cadd_phred' in tc:
                    cadd_phred = tc.get('cadd_phred')
            
            # Note: Conservation not cleanly returned by VEP ID without specific plugins
            vep_map[rsid] = {
                'cadd_phred': cadd_phred,
                'sift': sift,
                'polyphen': polyphen
            }
            
    print("Conservation not returned via VEP; requires region endpoint.")
    
    def get_vep(rsid, field):
        return vep_map.get(rsid, {}).get(field)
        
    df = df.copy()
    df['cadd_phred'] = df['rsid'].apply(lambda x: get_vep(x, 'cadd_phred'))
    df['sift_score'] = df['rsid'].apply(lambda x: get_vep(x, 'sift'))
    df['polyphen_score'] = df['rsid'].apply(lambda x: get_vep(x, 'polyphen'))
    
    # Save output
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    logger.info(f"Phase C: Annotated CSV saved to {output_path}")
    
    # Print stats
    print("\n========= MODULE 2 SUMMARY STATASTICS =========")
    print("\n--- Domain Distribution ---")
    print(df['domain_region'].value_counts(dropna=False).to_string())
    
    tm_perc = (df['is_transmembrane'].sum() / len(df)) * 100
    pore_perc = (df['is_pore_region'].sum() / len(df)) * 100
    print(f"\n% in Transmembrane region: {tm_perc:.2f}%")
    print(f"% in Pore region (M2): {pore_perc:.2f}%")
    
    print("\n--- CADD Summary Stats ---")
    df_cadd = pd.to_numeric(df['cadd_phred'], errors='coerce')
    print(df_cadd.describe().to_string())
    print("===============================================\n")

    return df

def main():
    root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    config_path = os.path.join(root_dir, "config", "parameters.yaml")
    log_file = os.path.join(root_dir, "logs", "module2_annotation.log")
    
    logger = setup_logging(log_file)
    logger.info("Initializing Module 2: Functional Annotation Layer")
    
    try:
        input_csv = os.path.join(root_dir, "data/processed/chrna7_rare_missense_variants.csv")
        df = pd.read_csv(input_csv)
        
        # Phase A
        corrected_csv_path = os.path.join(root_dir, "data/processed/chrna7_missense_master_corrected.csv")
        df_corrected = phase_a_coordinate_correction(df, logger, corrected_csv_path)
        
        # Phase B
        df_domain = phase_b_structural_domains(df_corrected, logger, config_path)
        
        # Phase C
        annotated_csv_path = os.path.join(root_dir, "data/processed/chrna7_missense_annotated.csv")
        df_final = phase_c_functional_scores(df_domain, logger, annotated_csv_path)
        
        logger.info("Module 2 Execution Complete.")
        
    except Exception as e:
        logger.error(f"Module 2 execution failed: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
