import os
import sys
import logging
import yaml
import requests
import pandas as pd
from typing import Dict, List, Any

# Ensure project structure allows relative imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def setup_logging(log_file: str) -> logging.Logger:
    """Configures project-level logging."""
    logger = logging.getLogger('module1_variant_mining')
    logger.setLevel(logging.INFO)

    # File handler
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger

def load_config(config_path: str) -> Dict[str, Any]:
    """Loads YAML configuration file."""
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def get_ensembl_gene_info(symbol: str, logger: logging.Logger) -> tuple[str, str, str]:
    """Fetches Gene ID, Canonical Transcript ID, and Translation ID from Ensembl."""
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?expand=1"
    headers = {"Content-Type": "application/json"}
    
    logger.info(f"API Call: Fetching target gene {symbol} details from {url}")
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    data = response.json()
    
    gene_id = data.get('id')
    canonical_transcript = None
    translation_id = None
    
    for transcript in data.get('Transcript', []):
        if transcript.get('is_canonical') == 1:
            canonical_transcript = transcript.get('id')
            translation_id = transcript.get('Translation', {}).get('id')
            break
            
    if not gene_id or not canonical_transcript or not translation_id:
        raise ValueError(f"Could not resolve IDs for {symbol}")
        
    logger.info(f"Resolved {symbol} to Gene: {gene_id}, Transcript: {canonical_transcript}, Translation: {translation_id}")
    return gene_id, canonical_transcript, translation_id

def get_transcript_variants(translation_id: str, logger: logging.Logger) -> List[Dict[str, Any]]:
    """Fetches all variants overlapping the canonical translation region."""
    url = f"https://rest.ensembl.org/overlap/translation/{translation_id}?feature=transcript_variation"
    headers = {"Content-Type": "application/json"}
    
    logger.info(f"API Call: Fetching transcript variants from {url}")
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    variants = response.json()
    logger.info(f"Fetched {len(variants)} overlapping sequences for {translation_id}")
    return variants

def get_variant_gnomad_af(rsids: List[str], logger: logging.Logger) -> Dict[str, float]:
    """Batch queries the gnomAD API or Ensembl variation API for allele frequencies."""
    url = "https://rest.ensembl.org/variation/homo_sapiens"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    logger.info(f"API Call: Batch fetching AF for {len(rsids)} RSIDs.")
    af_map = {}
    
    # Ensembl Batch API limits POST payloads
    batch_size = 200
    for i in range(0, len(rsids), batch_size):
        batch = rsids[i:i+batch_size]
        payload = {"ids": batch}
        
        response = requests.post(url, headers=headers, json=payload)
        
        if response.status_code != 200:
            logger.warning(f"Failed to fetch batch {i}-{i+batch_size}: {response.status_code}")
            continue
            
        data = response.json()
        
        for rsid, info in data.items():
            if not info:
                 continue
            
            freq = 0.0
            for col in info.get("colocated_variants", []):
               # we check population frequencies if available
               if "frequencies" in col:
                   for allele, pops in col["frequencies"].items():
                       # Attempt to get gnomAD wide frequency
                       af = pops.get("gnomad", pops.get("gnomad", 0.0))
                       if af:
                           freq = max(freq, float(af))
                           
            af_map[rsid] = freq
            
    logger.info(f"Resolved AF for {len(af_map)} variants.")
    return af_map

def process_variants(variants: List[Dict[str, Any]], gene_symbol: str, transcript_id: str, logger: logging.Logger, af_threshold: float) -> pd.DataFrame:
    """Filters for missense, extracts required columns, and queries AF."""
    extracted = []
    
    # Filter missense and extract columns
    for var in variants:
        consequence = var.get('type')
        if not consequence or 'missense_variant' not in consequence:
            continue
            
        rsid = var.get('id')
        if not rsid or not rsid.startswith('rs'):
             continue
             
        # Extract features
        row = {
            "gene": gene_symbol,
            "transcript_id": transcript_id,
            "rsid": rsid,
            "chrom": var.get('seq_region_name'),
            "pos": var.get('start'),
            "ref": var.get('allele', '').split('/')[0] if '/' in var.get('allele', '') else '',
            "alt": var.get('allele', '').split('/')[-1] if '/' in var.get('allele', '') else '',
            "protein_position": var.get('start'), # in translation, start=end=prot pos
            "amino_acid_change": f"{var.get('residues', '').split('/')[0] if '/' in var.get('residues', '') else ''}{var.get('start')}{var.get('residues', '').split('/')[-1] if '/' in var.get('residues', '') else ''}",
            "consequence": "missense_variant",
        }
        extracted.append(row)
        
    df = pd.DataFrame(extracted)
    logger.info(f"Filtered to {len(df)} missense variants with RSIDs.")
    
    if df.empty:
         return df
    
    # Fetch AFs
    rsids = df['rsid'].tolist()
    af_map = get_variant_gnomad_af(rsids, logger)
    
    df['gnomAD_AF'] = df['rsid'].map(af_map).fillna(0.0)
    
    # Filter by Rare frequency
    df_rare = df[df['gnomAD_AF'] < af_threshold].copy()
    logger.info(f"Filtered down to {len(df_rare)} rare variants (AF < {af_threshold}).")
    
    return df_rare

def main():
    root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    config_path = os.path.join(root_dir, "config", "parameters.yaml")
    log_file = os.path.join(root_dir, "logs", "module1_api.log")
    
    # Setup
    logger = setup_logging(log_file)
    logger.info("Initializing Module 1: Variant Mining")
    
    try:
        config = load_config(config_path)
        gene_symbol = config.get("gene_symbol", "CHRNA7")
        af_threshold = config.get("af_threshold", 0.001)
        output_dir = os.path.join(root_dir, config.get("output_dir", "data/processed"))
        
        # execution
        gene_id, transcript_id, translation_id = get_ensembl_gene_info(gene_symbol, logger)
        variants = get_transcript_variants(translation_id, logger)
        df_rare = process_variants(variants, gene_symbol, transcript_id, logger, af_threshold)
        
        # Save output
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f"{gene_symbol.lower()}_rare_missense_variants.csv")
        df_rare.to_csv(output_file, index=False)
        logger.info(f"Successfully saved results to {output_file}")
        
    except Exception as e:
        logger.error(f"Module 1 execution failed: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
