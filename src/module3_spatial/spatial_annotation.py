"""
Module 3: Spatial Structural Annotation (RAREMISS)

Three-state annotation semantics
----------------------------------
Spatial feature columns use an explicit three-state representation:

    True      = feature assessed POSITIVE  (variant is at the feature)
    False     = feature assessed NEGATIVE  (variant is definitely not at the feature)
    'unknown' = feature NOT assessable     (residue spatially unresolved, or
                                            no ligand available for binding-site)

IMPORTANT: 'unknown' must NEVER be treated as False.
  - Unresolved variants are excluded from enrichment contingency cells.
  - Unresolved variants receive 0 structural scoring points.

CSV round-trip: use src.raremiss_utils.parse_feature_cols() after any pd.read_csv()
call on files produced by this module.
"""

import os
import sys
import requests
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
import ssl
import urllib.request
import warnings
from Bio import BiopythonWarning
import yaml
import json

warnings.simplefilter('ignore', BiopythonWarning)
ssl._create_default_https_context = ssl._create_unverified_context

# Add project root to path for raremiss_utils import
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from src.raremiss_utils import validate_af_threshold, get_af_threshold, DEFAULT_AF_THRESHOLD

config_path = os.path.join(root_dir, "config", "parameters.yaml")
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)

pdb_id = config.get("structure_id", "7kox")
uniprot_id = config.get("uniprot_id")
if not uniprot_id:
    raise ValueError("Missing uniprot_id in config")
gene_symbol = config.get("gene_symbol", "CHRNA7").lower()
af_threshold = get_af_threshold(config, source="config/parameters.yaml")
ligand_mode = config.get("ligand_mode", "auto")

raw_dir = os.path.join(root_dir, "data/raw")
os.makedirs(raw_dir, exist_ok=True)
pdb_file = os.path.join(raw_dir, f"{pdb_id.upper()}.pdb")

print("STEP 1 - Structure Retrieval")
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
chains = sorted(set([chain.id for model in structure for chain in model]))

ligands = set()
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[0] != ' ' and residue.id[0] != 'W' and residue.resname != 'HOH':
                ligands.add(residue.resname)

print(f"Resolution: {resolution} Å  |  Method: {exp_method}")
print(f"Chains: {', '.join(chains)}")
print(f"Ligands: {', '.join(ligands) if ligands else 'None'}")
print(f"Chain used for spatial annotation: A (default for homomers; all chains: {', '.join(chains)})")

print("\nSTEP 2 - SIFTS Mapping")
sifts_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
sifts_res = requests.get(sifts_url)
if sifts_res.status_code != 200:
    print("Failed to fetch SIFTS mapping.")
    sys.exit(1)

target_chain_id = config.get("chain_id", "A")
if target_chain_id not in chains and len(chains) > 0:
    target_chain_id = chains[0]

uniprot_to_pdb = {}
unp_mappings = sifts_res.json().get(pdb_id.lower(), {}).get("UniProt", {}).get(uniprot_id, {}).get("mappings", [])
for m in unp_mappings:
    if m["chain_id"] == target_chain_id:
        unp_start, unp_end = m["unp_start"], m["unp_end"]
        pdb_start, pdb_end = m["start"]["residue_number"], m["end"]["residue_number"]
        if (unp_end - unp_start) == (pdb_end - pdb_start):
            for i in range(unp_end - unp_start + 1):
                uniprot_to_pdb[unp_start + i] = pdb_start + i
        else:
            print(f"Warning: Segment length mismatch — UNP {unp_start}-{unp_end} "
                  f"({unp_end-unp_start+1} aa) vs PDB {pdb_start}-{pdb_end} "
                  f"({pdb_end-pdb_start+1} aa). Segment skipped; residues in this "
                  f"range marked unresolved.")

input_csv = os.path.join(root_dir, f"data/processed/{gene_symbol}_missense_annotated.csv")
df = pd.read_csv(input_csv)
variant_positions = df['protein_position'].dropna().astype(int).unique()
mappable = [p for p in variant_positions if p in uniprot_to_pdb]
unmapped = [p for p in variant_positions if p not in uniprot_to_pdb]
coverage = (len(mappable) / len(variant_positions)) * 100 if len(variant_positions) > 0 else 0

print(f"Coverage: {coverage:.2f}% ({len(mappable)}/{len(variant_positions)} positions mapped)")
print(f"Unmapped positions: {len(unmapped)}")
if coverage < 85:
    print(f"Warning: Coverage {coverage:.2f}% is below 85%")

mapping_report = {
    "gene_symbol": gene_symbol.upper(),
    "af_threshold_used": float(af_threshold) if 'af_threshold' in dir() and af_threshold is not None else float(config.get('af_threshold', 0.01)),
    "total_variants": int(len(variant_positions)),
    "mapped_variants": int(len(mappable)),
    "unmapped_variants": int(len(unmapped)),
    "mapping_coverage_percentage": float(coverage),
    "ligand_available": len(ligands) > 0 and ligand_mode != "force_off",
    "ligand_identifiers": list(ligands),
    "ligand_mode": ligand_mode,
    "chain_used": target_chain_id,
    "all_chains": chains,
    "pdb_id": pdb_id.upper(),
    "resolution_angstrom": resolution,
    "experimental_method": exp_method,
}
with open(os.path.join(root_dir, f"data/processed/{gene_symbol}_mapping_report.json"), "w") as jf:
    json.dump(mapping_report, jf, indent=4)

print("\nSTEP 3 - Spatial Flags")
model = structure[0]
target_chain = model[target_chain_id]
chain_atoms = list(target_chain.get_atoms())
other_chain_atoms = [a for chain in model for a in chain.get_atoms() if chain.id != target_chain_id]
ligand_atoms = [
    a for chain in model for res in chain for a in res.get_atoms()
    if res.id[0] != ' ' and res.id[0] != 'W' and res.resname != 'HOH'
]

ns_other_chains = NeighborSearch(other_chain_atoms) if other_chain_atoms else None
ns_ligands = NeighborSearch(ligand_atoms) if ligand_atoms else None


def is_mappable(unp_pos):
    return pd.notna(unp_pos) and int(unp_pos) in uniprot_to_pdb


def check_is_binding_site(unp_pos):
    """
    Returns:
        True      — residue is within 5 Å of a hetero-residue atom
        False     — residue is mapped and NOT within 5 Å
        'unknown' — residue is not mappable, OR ligand_mode is force_off,
                    OR no hetero-residues present
    Note: ligand_mode='auto' includes crystallographic additives (SO4, glycerol,
    ions).  The caller must apply domain knowledge to interpret biological
    relevance.
    """
    if ligand_mode == "force_off" or len(ligands) == 0:
        return 'unknown'
    if not is_mappable(unp_pos):
        return 'unknown'
    pdb_pos = uniprot_to_pdb.get(int(unp_pos))
    if pdb_pos is None or ns_ligands is None:
        return 'unknown'
    try:
        res = target_chain[(' ', pdb_pos, ' ')]
        for atom in res.get_atoms():
            if ns_ligands.search(atom.coord, 5.0):
                return True
    except KeyError:
        return 'unknown'
    return False


def check_is_interface(unp_pos):
    """
    Returns:
        True      — residue is within 5 Å of an atom from another chain
        False     — residue is mapped and NOT within 5 Å
        'unknown' — residue is not mappable
    """
    if not is_mappable(unp_pos):
        return 'unknown'
    pdb_pos = uniprot_to_pdb.get(int(unp_pos))
    if pdb_pos is None or ns_other_chains is None:
        return False  # no other chains → interface impossible → False (not unknown)
    try:
        res = target_chain[(' ', pdb_pos, ' ')]
        for atom in res.get_atoms():
            if ns_other_chains.search(atom.coord, 5.0):
                return True
    except KeyError:
        return 'unknown'
    return False


# Compute flags
df['spatially_unresolved'] = ~df['protein_position'].apply(is_mappable)

df['is_binding_site'] = df['protein_position'].apply(
    lambda x: check_is_binding_site(x) if pd.notna(x) else 'unknown'
)
df['is_interface'] = df['protein_position'].apply(
    lambda x: check_is_interface(x) if pd.notna(x) else 'unknown'
)

# is_tm_core: only meaningful for mapped residues that are transmembrane
df['is_tm_core'] = df['is_transmembrane'].fillna(False).astype(bool) & ~df['spatially_unresolved']

# Cast to object before assigning mixed-type values
df['is_binding_site'] = df['is_binding_site'].astype(object)
df['is_interface']    = df['is_interface'].astype(object)
df['is_tm_core']      = df['is_tm_core'].astype(object)
if 'is_pore_region' in df.columns:
    df['is_pore_region'] = df['is_pore_region'].astype(object)
if 'is_transmembrane' in df.columns:
    df['is_transmembrane'] = df['is_transmembrane'].astype(object)

# Overwrite all spatial feature values for unresolved residues with 'unknown'
# True = feature assessed POSITIVE
# False = feature assessed NEGATIVE
# 'unknown' = cannot assess (unresolved)
unresolved_mask = df['spatially_unresolved']
df.loc[unresolved_mask, 'is_binding_site'] = 'unknown'
df.loc[unresolved_mask, 'is_interface']    = 'unknown'
df.loc[unresolved_mask, 'is_tm_core']      = 'unknown'
if 'is_pore_region' in df.columns:
    df.loc[unresolved_mask, 'is_pore_region'] = 'unknown'
if 'is_transmembrane' in df.columns:
    df.loc[unresolved_mask, 'is_transmembrane'] = 'unknown'

if ligand_mode == "force_off" or len(ligands) == 0:
    print("Warning: No ligand detected — binding site annotation skipped")

output_file = os.path.join(root_dir, f"data/processed/{gene_symbol}_missense_spatial_annotated.csv")
df.to_csv(output_file, index=False)
print(f"Saved: {output_file}")

print("\nSTEP 4 - Output")
bs_count  = df['is_binding_site'].eq(True).sum()
int_count = df['is_interface'].eq(True).sum()
tm_count  = df['is_tm_core'].eq(True).sum()
unres_count = df['spatially_unresolved'].eq(True).sum()
unk_bs    = df['is_binding_site'].eq('unknown').sum()
unk_int   = df['is_interface'].eq('unknown').sum()

print(f"Total variants      : {len(df)}")
print(f"Spatially unresolved: {unres_count} ({unres_count/len(df)*100:.1f}%)")
print(f"binding_site=True   : {bs_count}  |  unknown: {unk_bs}")
print(f"interface=True      : {int_count}  |  unknown: {unk_int}")
print(f"tm_core=True        : {tm_count}")
