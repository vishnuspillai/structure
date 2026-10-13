# RAREMISS: Validation Plan & Benchmark Strategy

**Version:** 2.0  
**Status:** Approved for Codebase Verification & Pre-Publication Audit  

---

## 1. Objectives & Verification Matrix

The RAREMISS Validation Plan outlines the protocols for verifying the computational integrity, reproducibility, and clinical translation of the pipeline.

| Verification Level | Scope | Success Criteria | Execution Method |
| :--- | :--- | :--- | :--- |
| **Level 1: Deterministic Code Integrity** | Regression tests (BUG-01 to BUG-09), formula boundaries, three-state CSV round-tripping | 100% test pass rate in `test_integrity.py` | Automated test suite (`python tests/test_integrity.py`) |
| **Level 2: Pipeline Integration & Runner** | End-to-end multi-gene execution (CHRNA7, BRCA1) | Zero unhandled exceptions, non-stale output generation, accurate logs | Automated runner (`python test_runner.py`) |
| **Level 3: Statistical Rigor** | Contingency matrices, FDR monotonicity, denominator tracking | Zero cell-count violations, adjusted $p \ge$ raw $p$ | Automated reference validation in Module 7 |
| **Level 4: Independent Clinical Benchmark** | Performance against ClinVar/HGMD diagnostic ground truths | Leakage-aware ROC/PR reporting, transparent disclosure of data deficits | Clinical benchmark module (`src/module7_validation/clinical_benchmark.py`) |

---

## 2. Benchmark Architecture: Leakage-Aware Evaluation

### 2.1 The Circularity Challenge in Computational Genomics
Evaluating computational pathogenicity algorithms on public variant repositories is frequently confounded by:
1. **Model circularity (Type I leakage):** Component predictors (e.g., CADD, PolyPhen-2, SIFT) have seen past ClinVar/HGMD variants during training.
2. **Internal consistency circularity (Type II leakage):** High-priority categories defined by structural features will trivially show structural feature enrichment.
3. **Ascertainment bias:** In ion channel genes like *CHRNA7*, reported missense variants in ClinVar are heavily skewed toward VUS (Variants of Uncertain Significance) or unreviewed submissions, with virtually zero verified Benign missense controls.

### 2.2 Benchmarking Protocol
For each target gene:
1. **Extract Curated Labels:**
   - Filter ClinVar variants for unambiguous annotations:
     - Positive class: `Pathogenic`, `Likely Pathogenic`
     - Negative class: `Benign`, `Likely Benign`
   - Explicitly discard: `Uncertain significance`, `Conflicting interpretations`, `Not provided`.
2. **Minimum Data Threshold:**
   - Evaluation proceeds **only** if $N_{\text{Pathogenic}} \ge 3$ and $N_{\text{Benign}} \ge 3$.
   - When sample counts are below 3, RAREMISS outputs status `insufficient_data` and suppresses ROC AUC / PR AUC calculations.
3. **Metrics Computed (When Sample Size Permitted):**
   - Area Under Receiver Operating Characteristic Curve (ROC AUC)
   - Area Under Precision-Recall Curve (PR AUC)
   - Diagnostic Sensitivity, Specificity, and Positive Predictive Value at threshold cutoffs ($S \ge 10$, $S \ge 6$).

---

## 3. Test Runner & Pipeline Isolation Guarantees

1. **Non-destructive Directory Clearing:**
   - The test runner must never purge the entire `data/processed/` directory.
   - Deletion is strictly scoped to the test gene prefix (`f"{gene.lower()}_*"`), ensuring reference datasets and other genes remain untouched.
2. **Stale Data Prevention:**
   - Pre-existing files for the active gene are cleared prior to running.
   - Output validation verifies that newly generated files bear the current run timestamp and valid file sizes.
3. **Runner Failure Reporting:**
   - Any non-zero exit code from submodules immediately halts execution and flags the run as `FAILED`. No false `PASS` is reported if intermediate stages fail.
