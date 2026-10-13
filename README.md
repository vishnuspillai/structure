# RAREMISS: Structure-Aware Prioritization of Rare Missense Variants

> [!IMPORTANT]
> **EXPERIMENTAL RESEARCH USE ONLY — NOT FOR CLINICAL DIAGNOSIS**  
> RAREMISS is an exploratory research framework designed for structural-genomic variant hypothesis generation. It is **NOT** a clinically validated diagnostic tool and must **NOT** be used for patient diagnosis, prognosis, treatment selection, or any clinical decision-making.

---

## Overview

RAREMISS (Rare Missense Variant Prioritization Framework) is an experimental computational pipeline for exploratory prioritization of rare missense variants by integrating:

- **Population Allele Frequency:** Filtered and scored via gnomAD (default threshold $\le 0.01$; ultra-rare scoring tiers $< 10^{-5}$ and $< 10^{-4}$)
- **In Silico Pathogenicity Predictors:** CADD (Phred), SIFT, and PolyPhen-2
- **3D Structural Context:** Solved coordinate structures (PDB) via SIFTS mapping, computing binding site proximity, multimer interface contacts, transmembrane core topology, and pore domains

Unlike pipelines that rely solely on 1D sequence conservation, RAREMISS incorporates 3D spatial coordinate models to examine whether variants cluster within key functional and architectural domains.

---

## Intended Use & Methodological Limitations

- **Software Engineering Status:** The software architecture, data retrieval pipelines, 3-state spatial coordinate logic, and Benjamini–Hochberg FDR multiple-testing corrections have been verified via a deterministic automated test suite (`python tests/test_integrity.py`).
- **Internal Consistency vs. Clinical Validation:** The structural enrichment calculations in Module 7 assess **internal consistency** between prioritized scores and structural features. Because structural features are directly included in the priority score, structural enrichment demonstrates that the prioritization reflects structural inputs; **it does NOT constitute independent clinical or biological validation of pathogenicity**.
- **Clinical Label Scarcity:** In many target proteins (including *CHRNA7*), public clinical databases (ClinVar) lack sufficient high-confidence benign missense variants ($\ge 1$-star review, no conflicts) to support rigorous prospective ROC AUC / PR AUC clinical benchmarking. In such cases, RAREMISS reports an `insufficient_data` evaluation status.
- **Predictor Training Overlap:** Component predictors (CADD v1.6, PolyPhen-2, SIFT) utilized historical variant catalogs during parameter training. Because those closed training sets cannot be audited locally, potential training overlap is an unresolved possibility. RAREMISS makes no claim of leakage-free clinical evaluation.

---

## Scoring Architecture (Heuristic Integer Formulation)

The RAREMISS Priority Score is an integer heuristic metric with a theoretical range of $[0, 12]$:

$$\text{Priority Score} = S_{\text{AF}} + S_{\text{Structure}} + S_{\text{CADD}} + S_{\text{PolyPhen}} + S_{\text{SIFT}}$$

| Component | Condition / Criterion | Points | Sub-Cap |
| :--- | :--- | :---: | :---: |
| **Allele Frequency ($S_{\text{AF}}$)** | $\text{AF} < 10^{-5}$<br>$10^{-5} \le \text{AF} < 10^{-4}$<br>$\text{AF} \ge 10^{-4}$ or Missing (`NaN`) | $+2$<br>$+1$<br>$0$ | Max 2 |
| **Structural Context ($S_{\text{Structure}}$)** | `is_binding_site == True`<br>`is_pore_region == True`<br>`is_interface == True`<br>`is_tm_core == True` (if not pore)<br>`False` or `'unknown'` | $+4$<br>$+4$<br>$+2$<br>$+1$<br>$0$ | **Max 4 (Hard Cap)** |
| **CADD Phred ($S_{\text{CADD}}$)** | $\text{CADD} \ge 30$<br>$25 \le \text{CADD} < 30$<br>$20 \le \text{CADD} < 25$<br>$\text{CADD} < 20$ or Missing | $+3$<br>$+2$<br>$+1$<br>$0$ | Max 3 |
| **PolyPhen-2 ($S_{\text{PolyPhen}}$)** | `probably_damaging`<br>`possibly_damaging`<br>`benign` or Missing | $+2$<br>$+1$<br>$0$ | Max 2 |
| **SIFT ($S_{\text{SIFT}}$)** | `deleterious`<br>`tolerated` or Missing | $+1$<br>$0$ | Max 1 |
| **Total Possible** | | | **Max 12** |

### Priority Risk Categories
- **High:** $\text{Score} \ge 10$
- **Medium:** $6 \le \text{Score} \le 9$
- **Low:** $3 \le \text{Score} \le 5$
- **Very Low:** $\text{Score} < 3$

---

## Installation & Setup

### Requirements
- Python 3.10+ (tested on Python 3.11)
- Node.js 18+ (for optional React/Mol* web frontend)

### Clone & Dependencies
```bash
git clone https://github.com/vishnuspillai/structure.git RAREMISS
cd RAREMISS
pip install -r requirements.txt
```

---

## Deterministic Verification & Testing

Before running the pipeline or web services, verify system integrity using the deterministic test suite:

```bash
python tests/test_integrity.py
```

Expected output: `145/145 passed (0 failed, 8 skipped)` (skipped tests correspond to absent external datasets for comparison genes).

---

## Running Locally

### 1. Command-Line Execution
Run the full automated pipeline directly from the command line:
```bash
python run_pipeline.py --gene [gene] --structure [pdb] --af 0.01
```
Parameters are loaded from `config/parameters.yaml` unless overridden via CLI flags.

### 2. Backend REST API
The FastAPI backend orchestrates pipeline execution, serves variant datasets, and hosts the AI assistant:
```bash
python src/api/main.py
```
The interactive API documentation is available at `http://localhost:8000/docs`.

### 3. Interactive Web Dashboard (Optional)
The React + Vite dashboard provides 3D Mol* molecular visualization and variant exploration:
```bash
cd ui
npm install
npm run dev
```
Dashboard available at `http://localhost:5173`.

---

## License & Attribution

This software is released for academic and research evaluation purposes only. See `LICENSE` for details.
