# RAREMISS Release Manifest

## Release Overview

| Attribute | Specification |
| :--- | :--- |
| **Product Name** | RAREMISS (Rare Missense Variant Prioritization Framework) |
| **Release Version** | **`v1.4.0-rc1`** *(Experimental Testing Pre-Release)* |
| **Previous Tagged Releases** | `v1.0.0`, `v1.1.0`, `v1.2.0`, `v1.3.0` |
| **Release Channel** | Pre-Release / Experimental Testing |
| **Target Architecture** | Platform-independent Python (tested on Windows 11 & Linux x86_64) |
| **Base Git Commit** | `b07e4f4` |
| **Python Version** | `>= 3.10` (verified on Python 3.11.9) |
| **Dependencies** | Listed in [requirements.txt](file:///c:/Users/Vishnu%20S%20Pillai/.gemini/antigravity/scratch/CHRNA7_project/requirements.txt) |
| **Regulatory & Clinical Status** | **EXPERIMENTAL RESEARCH USE ONLY — NOT VALIDATED FOR CLINICAL DIAGNOSIS** |

---

## Intended Purpose & Operational Scope

RAREMISS is designed exclusively as an exploratory research framework for structural-genomic variant hypothesis generation. It computes multi-parametric heuristic priority rankings by integrating population allele frequency, computational pathogenicity scores, and 3D macromolecular coordinate models.

### Clinical Disclaimer
RAREMISS is **not** a diagnostic test, has **not** been cleared by any regulatory authority, and must **not** be used to guide clinical care, medical diagnosis, patient screening, or therapy selection.

---

## Canonical Configuration Defaults

The frozen default configuration is defined in [config/parameters.yaml](file:///c:/Users/Vishnu%20S%20Pillai/.gemini/antigravity/scratch/CHRNA7_project/config/parameters.yaml):

```yaml
af_threshold: 0.01
consequence_filter: missense_variant
gene_symbol: CHRNA7
genome_build: GRCh38
log_dir: logs
output_dir: data/processed
species: homo_sapiens
structure_id: 7KOX
transcript_id: ENST00000305886
translation_id: ENSP00000305018
uniprot_id: P36544
```

### Parameter Validation Rules
- `af_threshold`: Must satisfy $0.0 < \text{threshold} \le 1.0$. Thresholds $\le 0.0$ or $> 1.0$ are strictly rejected by [validate_af_threshold()](file:///c:/Users/Vishnu%20S%20Pillai/.gemini/antigravity/scratch/CHRNA7_project/src/raremiss_utils.py).
- `consequence_filter`: Defaults to `missense_variant`.
- `species`: Defaults to `homo_sapiens`.

---

## Scoring System Specification (Frozen Formulation)

The RAREMISS Priority Score is an integer heuristic metric on a bounded theoretical scale of $[0, 12]$:

$$\text{Priority Score} = S_{\text{AF}} + S_{\text{Structure}} + S_{\text{CADD}} + S_{\text{PolyPhen}} + S_{\text{SIFT}}$$

| Component | Condition / Criterion | Assigned Points | Component Sub-Cap |
| :--- | :--- | :---: | :---: |
| **Allele Frequency ($S_{\text{AF}}$)** | $\text{AF} < 10^{-5}$<br>$10^{-5} \le \text{AF} < 10^{-4}$<br>$\text{AF} \ge 10^{-4}$ or Missing (`NaN`) | $+2$<br>$+1$<br>$0$ | Max 2 |
| **Structural Context ($S_{\text{Structure}}$)** | `is_binding_site == True`<br>`is_pore_region == True`<br>`is_interface == True`<br>`is_tm_core == True` (if not pore)<br>`False` or `'unknown'` | $+4$<br>$+4$<br>$+2$<br>$+1$<br>$0$ | **Max 4 (Hard Cap)** |
| **CADD Phred ($S_{\text{CADD}}$)** | $\text{CADD} \ge 30$<br>$25 \le \text{CADD} < 30$<br>$20 \le \text{CADD} < 25$<br>$\text{CADD} < 20$ or Missing | $+3$<br>$+2$<br>$+1$<br>$0$ | Max 3 |
| **PolyPhen-2 ($S_{\text{PolyPhen}}$)** | `probably_damaging`<br>`possibly_damaging`<br>`benign` or Missing | $+2$<br>$+1$<br>$0$ | Max 2 |
| **SIFT ($S_{\text{SIFT}}$)** | `deleterious`<br>`tolerated` or Missing | $+1$<br>$0$ | Max 1 |
| **Theoretical Maximum** | | | **Max 12** |

### Priority Risk Categories
- **High:** $\text{Score} \ge 10$
- **Medium:** $6 \le \text{Score} \le 9$
- **Low:** $3 \le \text{Score} \le 5$
- **Very Low:** $\text{Score} < 3$

> **Reconciliation Note:** Earlier narrative summaries in draft QA logs referenced an obsolete fractional point formulation (e.g. SIFT 1.5, PolyPhen 1.5/1.0, Pore 2.5, categories $\ge 10, 7-9, 4-6, <4$). The production codebase, utility functions ([src/raremiss_utils.py](file:///c:/Users/Vishnu%20S%20Pillai/.gemini/antigravity/scratch/CHRNA7_project/src/raremiss_utils.py)), test suite ([tests/test_integrity.py](file:///c:/Users/Vishnu%20S%20Pillai/.gemini/antigravity/scratch/CHRNA7_project/tests/test_integrity.py)), generated output datasets, and formal specification ([docs/SCIENTIFIC_METHODS.md](file:///c:/Users/Vishnu%20S%20Pillai/.gemini/antigravity/scratch/CHRNA7_project/docs/SCIENTIFIC_METHODS.md)) are fully reconciled and frozen on the clean integer $[0, 12]$ formulation above.

---

## Verification & Test Evidence

### Deterministic Test Suite
- **Command:** `python tests/test_integrity.py`
- **Result:** **145 Passed | 0 Failed | 8 Skipped**
- **Exit Status:** `0`
- **Coverage:** Exercises production AF threshold validation, 3-state spatial representation, scoring boundary edge cases, category mapping, column preservation, Benjamini–Hochberg FDR correction, clinical benchmark selection bias guards, and run provenance metadata generation.

### Smoke Test Battery
1. **CLI Help:** `python run_pipeline.py --help` $\rightarrow$ Exit 0.
2. **Prioritization Execution:** `python src/module4_prioritization/prioritization.py` $\rightarrow$ Exit 0; generates `chrna7_ranked_variants.csv` (438 variants).
3. **Mechanistic Metrics:** `python src/module5_mechanistic/mechanistic_metrics.py` $\rightarrow$ Exit 0; generates `chrna7_top10_structural_metrics.csv`.
4. **Structural Enrichment:** `python src/module7_validation/structural_enrichment_ci.py` $\rightarrow$ Exit 0; evaluates Fisher $2\times 2$ contingency tables with cell count safeguards and emits `chrna7_enrichment_results.json`.
5. **Synthetic Clinical Benchmark:** Evaluates Wilcoxon ROC AUC rank-sum calculation with `validation_scope: "test_of_metric_calculation_only"` $\rightarrow$ Exit 0.
6. **Backend REST API:** `TestClient(app).get('/config')` $\rightarrow$ 200 OK; `TestClient(app).get('/data/chrna7_ranked_variants.csv')` $\rightarrow$ 200 OK (438 rows); `TestClient(app).get('/steps')` $\rightarrow$ 200 OK.

---

## Known Scientific & Technical Limitations

1. **Absence of Independent Clinical Validation:**
   - In *CHRNA7*, ClinVar contains zero high-confidence benign/likely benign missense variants meeting $\ge 1$-star review criteria without conflicts.
   - Independent ROC AUC and PR AUC clinical discrimination cannot be empirically calculated without both true positive and true negative cohorts. The pipeline reports `insufficient_data` for *CHRNA7*.
2. **Unresolved Predictor Training Overlap (Data Leakage):**
   - Embedded sequence predictors (CADD v1.6, PolyPhen-2, SIFT) utilized historical clinical variant databases during model training.
   - Because those training sets are closed and not bundled locally, overlap cannot be verified or excluded from available local evidence. RAREMISS explicitly records this limitation as `unresolved_possibility`.
3. **Internal Score Consistency vs. Independent Evidence:**
   - Module 7 structural enrichment tests confirm that priority scores correlate with structural annotations. Because structural annotations are directly included in the score, this is an **internal consistency check**, not an independent proof of clinical pathogenicity.
4. **Coordinate Coverage & Electron Density:**
   - Disordered loops and N/C-terminal residues lacking PDB atom coordinates are represented as `'unknown'` (spatially unresolved) and excluded from spatial feature enrichment denominators.
5. **Historical Git Repository Tracking of `venv/`:**
   - Historical commit `b07e4f4` previously tracked 9,238 files under `venv/`. The local working directory has cleaned these files from disk, and `.gitignore` prevents future tracking. A clean Git commit/push must not re-stage virtual environment directories.

---

## File Inventory (Release-Critical Components)

```
.
├── README.md                                      # Primary documentation & clinical disclaimer
├── RELEASE_MANIFEST.md                            # Release manifest & frozen specifications
├── requirements.txt                               # Pinned Python package dependencies
├── run_pipeline.py                                # End-to-end pipeline execution entrypoint
├── config/
│   └── parameters.yaml                            # Canonical run configuration
├── src/
│   ├── raremiss_utils.py                          # Core utility library (scoring, AF, FDR, metadata)
│   ├── module1_variant_mining/variant_mining.py   # Module 1: Ensembl API variant extraction
│   ├── module2_annotation/annotation_layer.py     # Module 2: VEP & UniProt annotation
│   ├── module3_spatial/spatial_annotation.py      # Module 3: SIFTS & 3D spatial mapping
│   ├── module4_prioritization/prioritization.py   # Module 4: Priority scoring & ranking
│   ├── module5_mechanistic/mechanistic_metrics.py # Module 5: Mechanistic structural metrics
│   ├── module6_clinical/clinical_review.py        # Module 6: Top-10 clinical literature review
│   ├── module7_validation/structural_enrichment.py# Module 7: Structural feature enrichment
│   ├── module7_validation/structural_enrichment_ci.py # Module 7 CI: Safeguarded enrichment
│   ├── module7_validation/clinical_benchmark.py   # Module 7b: Leakage-aware clinical benchmark
│   └── api/
│       ├── main.py                                # FastAPI backend REST endpoints
│       └── orchestrator.py                        # Pipeline process orchestrator
├── tests/
│   └── test_integrity.py                          # Automated deterministic verification suite
└── docs/
    ├── REPAIR_AND_AUDIT_REPORT.md                 # Complete audit log & defect resolution history
    ├── SCIENTIFIC_METHODS.md                      # Mathematical & biological specifications
    └── VALIDATION_PLAN.md                         # Validation protocol & leakage controls
```
