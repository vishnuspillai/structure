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

## Setup and Usage

1. **Environment Setup**: Create a virtual environment and install dependencies:
   ```bash
   python -m venv venv
   source venv/bin/activate  # Or venv\Scripts\activate on Windows
   pip install -r requirements.txt
   ```

2. **Pipeline Execution Workflow**: The framework is deeply modular. Execute the scripts sequentially from the root directory:

   ```bash
   # Module 1: Extract and filter rare baseline variants
   python src/module1_variant_mining/variant_mining.py
   
   # Module 2: Phase A/B/C Annotations (Fixes coordinates, assigns domains, fetches scores)
   python src/module2_annotation/annotation_layer.py
   
   # Module 3: Structural Mapping against PDB 7KOX
   python src/module3_spatial/spatial_annotation.py
   
   # Module 4: Calculate Priority Scores and Categorize 
   python src/module4_prioritization/prioritization.py
   
   # Module 5: Isolate Mechanistic Distances for the Top 10 variants
   python src/module5_mechanistic/mechanistic_metrics.py
   
   # Module 6: Query NCBI ClinVar and PubMed for Top 10 Evidence
   python src/module6_clinical/clinical_review.py
   
   # Module 7: Validate Feature Enrichment via Fisher Exact Test
   python src/module7_validation/structural_enrichment_ci.py
   ```

3. **Outputs**: Track logs in `logs/` and review generated tabular payload datasets in `data/processed/`.

## Pipeline Version
**v1.0.0** - First validated structural-genomic prioritization framework release.
