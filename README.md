# CHRNA7 Variant Structural Prioritization Framework

A fully reproducible, modular bioinformatics pipeline to prioritize rare missense variants in the human CHRNA7 gene using public datasets (Ensembl, UniProt, RCSB PDB, ClinVar, PubMed).

## Project Structure
- `data/`: Raw downloaded data and pipeline processed files.
- `src/`: Modular codebase encompassing 7 core pipeline stages.
  - `module1_variant_mining`: Retrieve variants from Ensembl REST API and filter by gnomAD AF < 0.001.
  - `module2_annotation`: Validate mapping, reconstruct domain bounds, and batch query CADD/SIFT/PolyPhen via VEP.
  - `module3_spatial`: Align subset of variables to 7KOX structure via PDBe SIFTS to compute mapping bounds.
  - `module4_prioritization`: Assign scoring matrices scaling structurally critical elements and pathogenic functional predictions.
  - `module5_mechanistic`: Compute high-precision 3D structural boundaries to ligands, chains, and overall center.
  - `module6_clinical`: Review clinical traits across ClinVar and map to PubMed literature presence.
  - `module7_validation`: Run rigorous Fisher Exact Tests evaluating pipeline feature enrichment accuracy against priority ranking.
- `logs/`: Application execution logs.
- `config/`: Pipeline YAML configuration parameters.

## Setup
1. Create a virtual environment and install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
2. Execute modules sequentially from the `src/` directory.

## Pipeline Version
**v1.0.0** - First validated structural-genomic prioritization framework release.
