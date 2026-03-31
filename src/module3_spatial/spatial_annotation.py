import os
import requests
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
import ssl
import urllib.request
import warnings
from Bio import BiopythonWarning
import yaml

warnings.simplefilter('ignore', BiopythonWarning)
ssl._create_default_https_context = ssl._create_unverified_context

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
config_path = os.path.join(root_dir, "config", "parameters.yaml")
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)

pdb_id = config.get("structure_id", "7kox")
uniprot_id = config.get("uniprot_id")
if not uniprot_id:
    raise ValueError("Missing uniprot_id in config")
gene_symbol = config.get("gene_symbol", "CHRNA7").lower()
ligand_mode = config.get("ligand_mode", "auto")

raw_dir = os.path.join(root_dir, "data/raw")
os.makedirs(raw_dir, exist_ok=True)
pdb_file = os.path.join(raw_dir, f"{pdb_id.upper()}.pdb")

print("=====================================================")
print("STEP 1 - Structure Retrieval")
print("=====================================================")

rcsb_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
res = requests.get(rcsb_url)
metadata = res.json()
resolution = metadata.get("rcsb_entry_info", {}).get("resolution_combined", ["N/A"])[0]
exp_method = metadata.get("exptl", [{}])[0].get("method", "N/A")

pdb_download_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
if not os.path.exists(pdb_file):
    with open(pdb_file, "w") as f:
        f.write(requests.get(pdb_download_url).text)

parser = PDBParser(QUIET=True)
structure = parser.get_structure(pdb_id.upper(), pdb_file)

chains = list(set([chain.id for model in structure for chain in model]))
chains.sort()

# Ligands
ligands = set()
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[0] != ' ' and residue.id[0] != 'W' and residue.resname != 'HOH': 
                ligands.add(residue.resname)

print(f"Resolution: {resolution} Angstroms")
print(f"Experimental method: {exp_method}")
print(f"Chains: {', '.join(chains)}")
print(f"Ligand identifiers: {', '.join(ligands)}")

print("\n=====================================================")
print("STEP 2 - SIFTS Mapping")
print("=====================================================")

sifts_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
sifts_res = requests.get(sifts_url)
if sifts_res.status_code != 200:
    print("Failed to fetch SIFTS mapping.")
    exit(1)

uniprot_to_pdb = {}
unp_mappings = sifts_res.json().get(pdb_id.lower(), {}).get("UniProt", {}).get(uniprot_id, {}).get("mappings", [])

for m in unp_mappings:
    if m["chain_id"] == "A":
        unp_start = m["unp_start"]
        unp_end = m["unp_end"]
        pdb_start = m["start"]["residue_number"]
        pdb_end = m["end"]["residue_number"]
        
        if (unp_end - unp_start) == (pdb_end - pdb_start):
            for i in range(unp_end - unp_start + 1):
                uniprot_to_pdb[unp_start + i] = pdb_start + i
        else:
            print(f"Warning: Segment mismatch - UNP {unp_start}-{unp_end} vs PDB {pdb_start}-{pdb_end}")

input_csv = os.path.join(root_dir, f"data/processed/{gene_symbol}_missense_annotated.csv")
df = pd.read_csv(input_csv)
variant_positions = df['protein_position'].dropna().astype(int).unique()

mappable = [p for p in variant_positions if p in uniprot_to_pdb]
unmapped = [p for p in variant_positions if p not in uniprot_to_pdb]

coverage = (len(mappable) / len(variant_positions)) * 100 if len(variant_positions) > 0 else 0

print(f"Mapping coverage: {coverage:.2f}% ({len(mappable)}/{len(variant_positions)} variants mapped)")
print(f"Unmapped positions: {len(unmapped)}")

if coverage < 85:
    print(f"Warning: Coverage {coverage:.2f}% is below 85%, proceeding according to policy.")

import json
mapping_report = {
    "total_variants": int(len(variant_positions)),
    "mapped_variants": int(len(mappable)),
    "unmapped_variants": int(len(unmapped)),
    "mapping_coverage_percentage": float(coverage)
}
with open(os.path.join(root_dir, f"data/processed/{gene_symbol}_mapping_report.json"), "w") as jf:
    json.dump(mapping_report, jf, indent=4)

print("\n=====================================================")
print("STEP 3 - Spatial Flags")
print("=====================================================")

model = structure[0]
chain_A = model["A"]

chain_A_atoms = list(chain_A.get_atoms())
other_chain_atoms = [a for chain in model for a in chain.get_atoms() if chain.id != "A"]
ligand_atoms = [a for chain in model for res in chain for a in res.get_atoms() if res.id[0] != ' ' and res.id[0] != 'W' and res.resname != 'HOH']

ns_other_chains = NeighborSearch(other_chain_atoms) if other_chain_atoms else None
ns_ligands = NeighborSearch(ligand_atoms) if ligand_atoms else None

def is_mappable(unp_pos):
    return pd.notna(unp_pos) and int(unp_pos) in uniprot_to_pdb

def check_is_binding_site(unp_pos):
    if ligand_mode == "force_off" or len(ligands) == 0:
        return 'unknown'
    if not is_mappable(unp_pos): return False
    pdb_pos = uniprot_to_pdb.get(int(unp_pos))
    if pdb_pos is None or ns_ligands is None: return False
    try:
        res = chain_A[(' ', pdb_pos, ' ')]
        for atom in res.get_atoms():
            if ns_ligands.search(atom.coord, 5.0):
                return True
    except KeyError:
        return False
    return False

def check_is_interface(unp_pos):
    if not is_mappable(unp_pos): return False
    pdb_pos = uniprot_to_pdb.get(int(unp_pos))
    if pdb_pos is None or ns_other_chains is None: return False
    try:
        res = chain_A[(' ', pdb_pos, ' ')]
        for atom in res.get_atoms():
            if ns_other_chains.search(atom.coord, 5.0):
                return True
    except KeyError:
        return False
    return False

df['spatially_unresolved'] = ~df['protein_position'].apply(is_mappable)
df['is_binding_site'] = df['protein_position'].apply(lambda x: check_is_binding_site(x) if pd.notna(x) else 'unknown')
df['is_interface'] = df['protein_position'].apply(lambda x: check_is_interface(x) if pd.notna(x) else False)
df['is_tm_core'] = df['is_transmembrane'].fillna(False).astype(bool) & ~df['spatially_unresolved']

# Enforce strict defaults for structurally unmapped variants to gracefully pass analysis safely
df['is_binding_site'] = df['is_binding_site'].astype(object)
df.loc[df['spatially_unresolved'], 'is_binding_site'] = 'unknown'
df.loc[df['spatially_unresolved'], 'is_interface'] = False
df.loc[df['spatially_unresolved'], 'is_tm_core'] = False
if 'is_pore_region' in df.columns:
    df.loc[df['spatially_unresolved'], 'is_pore_region'] = False
if 'is_transmembrane' in df.columns:
    df.loc[df['spatially_unresolved'], 'is_transmembrane'] = False

if ligand_mode == "force_off" or len(ligands) == 0:
    print("Warning: No ligand detected (or forced off) — binding site annotation skipped")

output_file = os.path.join(root_dir, f"data/processed/{gene_symbol}_missense_spatial_annotated.csv")
df.to_csv(output_file, index=False)
print(f"Saved spatial annotations to: {output_file}")
print("Finished calculating spatial flags.")

print("\n=====================================================")
print("STEP 4 - Output")
print("=====================================================")

bs_perc = (df['is_binding_site'].eq(True).sum() / len(df)) * 100
int_perc = (df['is_interface'].eq(True).sum() / len(df)) * 100
tm_perc = (df['is_tm_core'].eq(True).sum() / len(df)) * 100
unres_perc = (df['spatially_unresolved'].eq(True).sum() / len(df)) * 100

print(f"% binding_site: {bs_perc:.2f}%")
print(f"% interface: {int_perc:.2f}%")
print(f"% tm_core: {tm_perc:.2f}%")
print(f"Spatially unresolved: {unres_perc:.2f}%")
