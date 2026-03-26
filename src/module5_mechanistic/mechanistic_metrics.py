import pandas as pd
import numpy as np
import requests
import os
import yaml
from Bio.PDB import PDBParser
import warnings
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
config_path = os.path.join(root_dir, "config", "parameters.yaml")
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)

gene_symbol = config.get("gene_symbol", "CHRNA7").lower()
structure_id = config.get("structure_id", "7kox").lower()
uniprot_id = config.get("uniprot_id")
ligand_mode = config.get("ligand_mode", "auto")

print("=====================================================")
print("STEP 1 - Extract Top 10")
print("=====================================================")

input_csv = os.path.join(root_dir, f'data/processed/{gene_symbol}_ranked_variants.csv')
df = pd.read_csv(input_csv)

top10 = df.head(10).copy()

# Keep requested columns
cols_to_keep = [
    'rsid', 'protein_position', 'amino_acid_change', 'domain_region',
    'gnomAD_AF', 'cadd_phred', 'sift_pred', 'polyphen_pred',
    'is_binding_site', 'is_interface', 'is_pore_region', 'is_tm_core',
    'spatially_unresolved', 'priority_score'
]

# Rename according to instructions
rename_dict = {
    'sift_pred': 'sift_prediction',
    'polyphen_pred': 'polyphen_prediction',
    'is_pore_region': 'is_pore_core'
}

top10 = top10[cols_to_keep].rename(columns=rename_dict)
top10_output = os.path.join(root_dir, f'data/processed/{gene_symbol}_top10_variants.csv')
top10.to_csv(top10_output, index=False)

print("=====================================================")
print(f"COMPUTING STRUCTURAL METRICS ({structure_id.upper()})")
print("=====================================================")

# Load PDB
pdb_file = os.path.join(root_dir, f'data/raw/{structure_id.upper()}.pdb')
parser = PDBParser(QUIET=True)
structure = parser.get_structure(structure_id.upper(), pdb_file)
model = structure[0]

# Load SIFTS mappings
sifts_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{structure_id}"
sifts_res = requests.get(sifts_url)
uniprot_to_pdb = {}
if sifts_res.status_code == 200:
    unp_mappings = sifts_res.json().get(structure_id, {}).get("UniProt", {}).get(uniprot_id, {}).get("mappings", [])
    for m in unp_mappings:
        if m["chain_id"] == "A":
            unp_start = m["unp_start"]
            unp_end = m["unp_end"]
            pdb_start = m["start"]["residue_number"]
            pdb_end = m["end"]["residue_number"]
            if (unp_end - unp_start) == (pdb_end - pdb_start):
                for i in range(unp_end - unp_start + 1):
                    uniprot_to_pdb[unp_start + i] = pdb_start + i

# Compute COM
all_heavy_coords = []
for chain in model:
    for res in chain:
        if res.id[0] == ' ':
            for atom in res:
                if atom.element != 'H':
                    all_heavy_coords.append(atom.coord)
pentamer_com = np.mean(all_heavy_coords, axis=0)

chain_A = model['A']
ligand_atoms = [a for chain in model for res in chain for a in res if res.id[0] != ' ' and res.id[0] != 'W' and res.resname != 'HOH' and a.element != 'H']
other_chain_atoms = [a for chain in model for res in chain for a in res if chain.id != 'A' and res.id[0] == ' ' and a.element != 'H']

def get_metrics(row):
    unp_pos = row['protein_position']
    domain = row['domain_region']
    
    if pd.isna(unp_pos) or row['spatially_unresolved']:
        return pd.Series([np.nan, np.nan, np.nan])
        
    pdb_pos = uniprot_to_pdb.get(int(unp_pos))
    if pdb_pos is None:
        return pd.Series([np.nan, np.nan, np.nan])
        
    try:
        target_res = chain_A[(' ', pdb_pos, ' ')]
    except KeyError:
        return pd.Series([np.nan, np.nan, np.nan])
        
    target_atoms = [a for a in target_res if a.element != 'H']
    if not target_atoms:
        return pd.Series([np.nan, np.nan, np.nan])
        
    if ligand_mode == "force_off" or len(ligand_atoms) == 0:
        min_ligand = np.nan
    else:
        min_ligand = min((np.linalg.norm(ta.coord - la.coord) for ta in target_atoms for la in ligand_atoms), default=np.nan)
        
    min_interface = min((np.linalg.norm(ta.coord - oa.coord) for ta in target_atoms for oa in other_chain_atoms), default=np.nan)
    
    radial_dist = np.nan
    if domain == 'M2':
        res_com = np.mean([a.coord for a in target_atoms], axis=0)
        radial_dist = np.linalg.norm(res_com - pentamer_com)
        
    return pd.Series([min_ligand, min_interface, radial_dist])

metrics = top10.apply(get_metrics, axis=1)
metrics.columns = ['min_dist_ligand', 'min_dist_interface', 'radial_dist_center']

final_df = pd.concat([top10, metrics], axis=1)

output_path = os.path.join(root_dir, f'data/processed/{gene_symbol}_top10_structural_metrics.csv')
final_df.to_csv(output_path, index=False)

print("\n=====================================================")
print("FINAL TABLE OUTPUT")
print("=====================================================")
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

output_txt = os.path.join(root_dir, 'output_mod5.txt')
with open(output_txt, 'w', encoding='utf-8') as f:
    f.write(final_df.to_string(index=False))

print("Metrics written to output_mod5.txt")
