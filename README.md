# RAREMISS: Structure-Aware Prioritization of Rare Missense Variants

##  Overview

RAREMISS (Rare Missense Variant Prioritization Framework) is a generalized computational framework for prioritizing rare missense variants by integrating:

- Population allele frequency (gnomAD)
- Functional prediction scores (CADD, SIFT, PolyPhen)
- Three-dimensional structural context (PDB)

Unlike conventional pipelines that rely solely on sequence-based predictions, RAREMISS incorporates protein structure to identify variants enriched in functionally critical regions.

---

##  Key Insight

We demonstrate that:

> Structural context enables protein-specific prioritization of rare variants, revealing that enrichment in functional regions is not universal but depends on protein architecture.



##  Features

- ✅ Gene-agnostic variant mining (Ensembl API)
- ✅ Allele-specific VEP annotation (CADD, SIFT, PolyPhen)
- ✅ Dynamic UniProt domain parsing
- ✅ Structure-aware mapping via SIFTS
- ✅ Spatial annotation:
  - binding site (ligand-aware / ligand-agnostic)
  - interface regions
  - pore/core regions
- ✅ Robust prioritization scoring framework
- ✅ Statistical validation (Fisher’s Exact Test with safeguards)
- ✅ Scalable batching for large genes
- ✅ Interactive web interface (Mol* visualization)

---

##  Validation Summary

RAREMISS was validated across structurally distinct proteins:

| Protein | Result |
|--------|--------|
| CHRNA7 | Strong enrichment in binding/pore regions |
| KCNQ2 | Moderate structural enrichment |
| SCN1A | No significant enrichment (biologically consistent) |

This demonstrates that structural enrichment is **protein-dependent**, not algorithmically biased.

---

##  Installation

```bash
git clone https://github.com/vishnuspillai/structure.git
cd RAREMISS
pip install -r requirements.txt
```

---

## ⌨️ Running Locally

### 1. Start the Backend API
The backend handles pipeline orchestration and data serving.
```bash
python src/api/main.py
```
The API will be available at `http://localhost:8000`.

### 2. Start the Frontend Dashboard
The frontend provides the interactive visualization and control panel.
```bash
cd ui
npm install
npm run dev
```
The dashboard will be available at `http://localhost:5173`.

### 3. CLI Execution (Alternative)
You can also run the full pipeline directly from the command line:
```bash
python run_pipeline.py
```
Configure your parameters in `config/parameters.yaml` before running.
