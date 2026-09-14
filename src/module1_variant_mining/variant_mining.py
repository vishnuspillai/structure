import os
import sys
import logging
import yaml
import requests
import pandas as pd
from typing import Dict, List, Any

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def setup_logging(log_file: str) -> logging.Logger:
    logger = logging.getLogger('module1_variant_mining')
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

def load_config(config_path: str) -> Dict[str, Any]:
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def get_session():
    session = requests.Session()
    from urllib3.util.retry import Retry
    from requests.adapters import HTTPAdapter
    retries = Retry(total=5, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504])
    session.mount('https://', HTTPAdapter(max_retries=retries))
    return session

def get_ensembl_gene_info(symbol: str, species: str, logger: logging.Logger) -> tuple[str, str, str, str]:
    url = f"https://rest.ensembl.org/lookup/symbol/{species}/{symbol}?expand=1"
    headers = {"Content-Type": "application/json"}
    logger.info(f"Fetching gene details for {symbol}")
    session = get_session()
    response = session.get(url, headers=headers, timeout=120)
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

    uniprot_id = None
    xrefs_url = f"https://rest.ensembl.org/xrefs/id/{translation_id}?external_db=Uniprot/SWISSPROT;content-type=application/json"
    x_res = session.get(xrefs_url, timeout=120)
    if x_res.status_code == 200 and len(x_res.json()) > 0:
        uniprot_id = x_res.json()[0]['primary_id']
    else:
        xrefs_url_tr = f"https://rest.ensembl.org/xrefs/id/{translation_id}?external_db=Uniprot/SPTREMBL;content-type=application/json"
        x_tr_res = session.get(xrefs_url_tr, timeout=120)
        if x_tr_res.status_code == 200 and len(x_tr_res.json()) > 0:
            uniprot_id = x_tr_res.json()[0]['primary_id']

    if not uniprot_id:
        raise ValueError(f"Could not resolve UniProt ID for {translation_id}")

    logger.info(f"Resolved: Gene={gene_id}, Transcript={canonical_transcript}, Translation={translation_id}, UniProt={uniprot_id}")
    return gene_id, canonical_transcript, translation_id, uniprot_id

def get_transcript_variants(translation_id: str, logger: logging.Logger) -> List[Dict[str, Any]]:
    url = f"https://rest.ensembl.org/overlap/translation/{translation_id}?feature=transcript_variation"
    headers = {"Content-Type": "application/json"}
    logger.info(f"Fetching transcript variants for {translation_id}")
    session = get_session()
    response = session.get(url, headers=headers, timeout=120)
    response.raise_for_status()
    variants = response.json()
    logger.info(f"Fetched {len(variants)} overlapping variants")
    return variants

def get_variant_gnomad_af(rsids: List[str], logger: logging.Logger) -> Dict[str, float]:
    url = "https://rest.ensembl.org/variation/homo_sapiens"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    logger.info(f"Fetching AF for {len(rsids)} RSIDs")
    af_map = {}
    session = get_session()
    batch_size = 200
    total_batches = (len(rsids) + batch_size - 1) // batch_size
    for i in range(0, len(rsids), batch_size):
        batch_num = (i // batch_size) + 1
        print(f"Processing batch {batch_num}/{total_batches}...")
        batch = rsids[i:i + batch_size]
        try:
            response = session.post(url, headers=headers, json={"ids": batch}, timeout=120)
            if response.status_code != 200:
                logger.warning(f"Failed batch {batch_num}: HTTP {response.status_code}")
                continue
            data = response.json()
        except Exception as e:
            logger.warning(f"Failed batch {batch_num}: {e}")
            continue
        for rsid, info in data.items():
            if not info:
                continue
            freq = 0.0
            for col in info.get("colocated_variants", []):
                if "frequencies" in col:
                    for allele, pops in col["frequencies"].items():
                        af = pops.get("gnomad", 0.0)
                        if af:
                            freq = max(freq, float(af))
            af_map[rsid] = freq
    logger.info(f"Resolved AF for {len(af_map)} variants")
    return af_map

def process_variants(variants: List[Dict[str, Any]], gene_symbol: str, transcript_id: str,
                     logger: logging.Logger, af_threshold: float, max_variants: int = None) -> pd.DataFrame:
    extracted = []
    for var in variants:
        consequence = var.get('type')
        if not consequence or 'missense_variant' not in consequence:
            continue
        rsid = var.get('id')
        if not rsid or not rsid.startswith('rs'):
            continue
        row = {
            "gene": gene_symbol,
            "transcript_id": transcript_id,
            "rsid": rsid,
            "chrom": var.get('seq_region_name'),
            "pos": var.get('start'),
            "ref": var.get('allele', '').split('/')[0] if '/' in var.get('allele', '') else '',
            "alt": var.get('allele', '').split('/')[-1] if '/' in var.get('allele', '') else '',
            "protein_position": var.get('start'),
            "amino_acid_change": (
                f"{var.get('residues', '').split('/')[0] if '/' in var.get('residues', '') else ''}"
                f"{var.get('start')}"
                f"{var.get('residues', '').split('/')[-1] if '/' in var.get('residues', '') else ''}"
            ),
            "consequence": "missense_variant",
        }
        extracted.append(row)

    df = pd.DataFrame(extracted)
    logger.info(f"Filtered to {len(df)} missense variants with RSIDs")
    print(f"Total variants: {len(df)}")

    if df.empty:
        return df

    if len(df) > 2000:
        logger.warning("High variant load detected — applying chunked processing")

    if max_variants and len(df) > max_variants:
        logger.info(f"Downsampling from {len(df)} to {max_variants} variants")
        df = df.sample(n=max_variants, random_state=42).copy()

    rsids = df['rsid'].tolist()
    af_map = get_variant_gnomad_af(rsids, logger)
    df['gnomAD_AF'] = df['rsid'].map(af_map).fillna(0.0)

    df_rare = df[df['gnomAD_AF'] < af_threshold].copy()
    logger.info(f"Filtered to {len(df_rare)} rare variants (AF < {af_threshold})")
    print(f"Rare variants: {len(df_rare)}")
    return df_rare

def main():
    root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    config_path = os.path.join(root_dir, "config", "parameters.yaml")
    log_file = os.path.join(root_dir, "logs", "module1_api.log")
    logger = setup_logging(log_file)
    logger.info("Initializing Module 1: Variant Mining")
    try:
        config = load_config(config_path)
        gene_symbol = config.get("gene_symbol", "CHRNA7")
        species = config.get("species", "homo_sapiens")
        af_threshold = config.get("af_threshold", 0.001)
        max_variants = config.get("max_variants", None)
        output_dir = os.path.join(root_dir, config.get("output_dir", "data/processed"))

        gene_id, transcript_id, translation_id, uniprot_id = get_ensembl_gene_info(gene_symbol, species, logger)

        config['transcript_id'] = transcript_id
        config['translation_id'] = translation_id
        config['uniprot_id'] = uniprot_id
        with open(config_path, 'w') as f:
            yaml.safe_dump(config, f)

        variants = get_transcript_variants(translation_id, logger)
        df_rare = process_variants(variants, gene_symbol, transcript_id, logger, af_threshold, max_variants)

        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f"{gene_symbol.lower()}_rare_missense_variants.csv")
        df_rare.to_csv(output_file, index=False)
        logger.info(f"Saved results to {output_file}")

    except Exception as e:
        logger.error(f"Module 1 failed: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
