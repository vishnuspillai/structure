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

def get_session():
    session = requests.Session()
    from urllib3.util.retry import Retry
    from requests.adapters import HTTPAdapter
    retries = Retry(total=5, backoff_factor=2, status_forcelist=[ 429, 500, 502, 503, 504 ])
    session.mount('https://', HTTPAdapter(max_retries=retries))
    return session

def phase_a_coordinate_correction(df: pd.DataFrame, logger: logging.Logger, output_path: str) -> pd.DataFrame:
    """Phase A: Fetch GRCh38 coordinates via /variation API, update chrom, pos and audit AF."""
    url = "https://rest.ensembl.org/variation/homo_sapiens?pops=1"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    rsids = df['rsid'].tolist()
    batch_size = 200
    coord_map = {}
    af_map = {}
    alt_map = dict(zip(df['rsid'], df['alt']))
    
    logger.info(f"Phase A: Batch querying {len(rsids)} RSIDs for coordinate correction and AF repair.")
    session = get_session()
    total_batches = (len(rsids) + batch_size - 1) // batch_size
    
    for i in range(0, len(rsids), batch_size):
        batch_num = (i // batch_size) + 1
        print(f"Phase A: Processing batch {batch_num}/{total_batches}...")
        batch = rsids[i:i+batch_size]
        payload = {"ids": batch}
        try:
            res = session.post(url, headers=headers, json=payload, timeout=60)
            if res.status_code != 200:
                logger.warning(f"Failed to fetch batch {batch_num}: {res.status_code}. Skipping.")
                continue
            res.raise_for_status()
            data = res.json()
        except Exception as e:
            logger.warning(f"Failed to fetch batch {batch_num} due to exception: {e}. Skipping.")
            continue
            
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
    
    # Rare AF filtering: drop variants where updated AF is > 0.001
    original_count = len(df)
    df = df[(df['gnomAD_AF'].isna()) | (df['gnomAD_AF'] < 0.001)].copy()
    logger.info(f"Phase A: AF re-filtering removed {original_count - len(df)} variants > 0.001 AF.")
    
    df.to_csv(output_path, index=False)
    logger.info(f"Phase A: Corrected CSV saved to {output_path}")
    return df

def fetch_uniprot_features(uniprot_id: str, logger: logging.Logger) -> Dict[str, Any]:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    logger.info(f"Phase B: Fetching canonical UniProt record for {uniprot_id}")
    session = get_session()
    domains = {}
    try:
        res = session.get(url, timeout=60)
        res.raise_for_status()
        data = res.json()
        
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
    except Exception as e:
        logger.warning(f"Phase B: Failed to parse UniProt domains for {uniprot_id} ({e}). Defaulting to unassigned.")
        
    return domains

def phase_b_structural_domains(df: pd.DataFrame, logger: logging.Logger, config_path: str, config: dict) -> pd.DataFrame:
    """Phase B: Parse UniProt limits, alter params, assign domains."""
    uniprot_id = config.get('uniprot_id')
    if not uniprot_id:
        raise ValueError("Missing uniprot_id in config")
        
    domains = fetch_uniprot_features(uniprot_id, logger)
    
    with open(config_path, 'r') as file:
        file_config = yaml.safe_load(file)

        
    file_config['domains'] = domains
    with open(config_path, 'w') as file:
        yaml.safe_dump(file_config, file)
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

def phase_c_functional_scores(df: pd.DataFrame, logger: logging.Logger, output_path: str, config: dict) -> pd.DataFrame:
    """Phase C: VEP querying for CADD and Conservation."""
    url = "https://rest.ensembl.org/vep/human/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    params = {"CADD": 1}
    
    rsids = df['rsid'].tolist()
    batch_size = 200
    vep_map = {}
    
    alt_map = dict(zip(df['rsid'], df['alt']))
    
    logger.info(f"Phase C: Batch querying {len(rsids)} RSIDs to VEP for CADD.")
    session = get_session()
    total_batches = (len(rsids) + batch_size - 1) // batch_size
    
    for i in range(0, len(rsids), batch_size):
        batch_num = (i // batch_size) + 1
        print(f"Phase C: Processing batch {batch_num}/{total_batches}...")
        batch = rsids[i:i+batch_size]
        payload = {"ids": batch}
        # post query string via url payload? No, VEP API accepts URL params in the URL itself
        req_url = f"{url}?CADD=1"
        try:
            res = session.post(req_url, headers=headers, json=payload, timeout=120)
            if res.status_code != 200:
                logger.warning(f"Failed to fetch batch {batch_num}: {res.status_code}. Skipping.")
                continue
            res.raise_for_status()
            data = res.json()
        except Exception as e:
            logger.warning(f"Failed to fetch batch {batch_num} due to exception: {e}. Skipping.")
            continue
            
        for item in data:
            rid = item.get('input') or item.get('id')
            alt_allele = alt_map.get(rid)
            
            cadd_phred = None
            sift = None
            sift_pred = None
            polyphen = None
            polyphen_pred = None
            
            has_match = False
            transcript_id = config.get('transcript_id')
            for tc in item.get('transcript_consequences', []):
                if tc.get('transcript_id') == transcript_id and tc.get('variant_allele') == alt_allele:
                    if 'sift_score' in tc:
                        sift = tc.get('sift_score')
                        sift_pred = tc.get('sift_prediction')
                    if 'polyphen_score' in tc:
                        polyphen = tc.get('polyphen_score')
                        polyphen_pred = tc.get('polyphen_prediction')
                    if 'cadd_phred' in tc:
                        cadd_phred = tc.get('cadd_phred')
                    has_match = True
                    break
                    
            # Fallback if no matching transcript was found
            if not has_match:
                for tc in item.get('transcript_consequences', []):
                    if tc.get('variant_allele') == alt_allele:
                        if 'sift_score' in tc and sift is None:
                            sift = tc.get('sift_score')
                            sift_pred = tc.get('sift_prediction')
                        if 'polyphen_score' in tc and polyphen is None:
                            polyphen = tc.get('polyphen_score')
                            polyphen_pred = tc.get('polyphen_prediction')
                        if 'cadd_phred' in tc:
                            cadd_phred = tc.get('cadd_phred')
                            
            if rid not in vep_map or has_match:
                vep_map[rid] = {
                    'cadd_phred': cadd_phred,
                    'sift': sift,
                    'sift_pred': sift_pred,
                    'polyphen': polyphen,
                    'polyphen_pred': polyphen_pred
                }
            
    print("Conservation not returned via VEP; requires region endpoint.")
    
    def get_vep(rsid, field):
        return vep_map.get(rsid, {}).get(field)
        
    df = df.copy()
    df['cadd_phred'] = df['rsid'].apply(lambda x: get_vep(x, 'cadd_phred'))
    df['sift_score'] = df['rsid'].apply(lambda x: get_vep(x, 'sift'))
    df['sift_pred'] = df['rsid'].apply(lambda x: get_vep(x, 'sift_pred'))
    df['polyphen_score'] = df['rsid'].apply(lambda x: get_vep(x, 'polyphen'))
    df['polyphen_pred'] = df['rsid'].apply(lambda x: get_vep(x, 'polyphen_pred'))
    
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
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
            
        gene_symbol = config.get('gene_symbol', 'CHRNA7').lower()
        
        input_csv = os.path.join(root_dir, f"data/processed/{gene_symbol}_rare_missense_variants.csv")
        df = pd.read_csv(input_csv)
        
        # Phase A
        corrected_csv_path = os.path.join(root_dir, f"data/processed/{gene_symbol}_missense_master_corrected.csv")
        df_corrected = phase_a_coordinate_correction(df, logger, corrected_csv_path)
        
        # Phase B
        df_domain = phase_b_structural_domains(df_corrected, logger, config_path, config)
        
        # Phase C
        annotated_csv_path = os.path.join(root_dir, f"data/processed/{gene_symbol}_missense_annotated.csv")
        df_final = phase_c_functional_scores(df_domain, logger, annotated_csv_path, config)
        
        logger.info("Module 2 Execution Complete.")
        
    except Exception as e:
        logger.error(f"Module 2 execution failed: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
