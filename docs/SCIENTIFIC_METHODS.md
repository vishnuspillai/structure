# RAREMISS: Scientific Methods & Computational Specification

**Version:** 2.0 (Post-Audit Release)  
**Authors/Affiliation:** Lead Computational Genetics & Bioinformatics Methods Group  
**Project Repository:** [RAREMISS](https://github.com/vishnuspillai/structure)

---

## 1. Abstract & Scope

RAREMISS is an automated pipeline designed to prioritize rare missense genomic variants by synthesizing population allele frequency (gnomAD), sequence-based pathogenicity conservation algorithms (CADD, PolyPhen-2, SIFT), and 3D macromolecular protein structure spatial coordinates (PDB/SIFTS).

This document details the exact mathematical formulas, discrete algorithmic definitions, missing data imputation policies, spatial coordinate geometry conventions, and multiple-testing statistical frameworks implemented in RAREMISS.

---

## 2. Allele Frequency (AF) Filtering & Threshold Policies

### 2.1 Default & Ultra-Rare Thresholds
- **Canonical Default Threshold:** `af_threshold: 0.01` (1% population frequency in gnomAD).
- **Ultra-Rare Threshold:** `af_threshold: 0.001` (0.1% population frequency), supported as an explicit configuration option for severe, highly penetrant early-onset Mendelian phenotypes.
- **Valid Boundary Range:** The threshold must strictly satisfy $0 < \text{AF} \le 1.0$.
- **Zero AF Policy:** An input value of `0.0` is explicitly rejected with a `ValueError`. Because gnomAD reports non-zero frequencies for observed variants, filtering on `AF == 0.0` would erroneously drop all verified variants or select only unsequenced sites.
- **Missing AF Policy:** Variants with missing (`NaN` or `None`) AF in gnomAD are **retained** in upstream mining stages to prevent losing novel, unobserved mutations. However, in scoring, missing AF is **not** rewarded with rarity bonus points (points = 0).

---

## 3. Heuristic Priority Scoring Architecture

The RAREMISS Priority Score is an integer heuristic ranking metric with a theoretical range of $[0, 12]$.

$$\text{Priority Score} = S_{\text{AF}} + S_{\text{Structure}} + S_{\text{CADD}} + S_{\text{PolyPhen}} + S_{\text{SIFT}}$$

### 3.1 Component Definitions & Weight Allocations

| Component | Condition / Criterion | Points | Sub-Cap |
| :--- | :--- | :---: | :---: |
| **Allele Frequency ($S_{\text{AF}}$)** | $\text{AF} < 10^{-5}$<br>$10^{-5} \le \text{AF} < 10^{-4}$<br>$\text{AF} \ge 10^{-4}$ or Missing (`NaN`) | $+2$<br>$+1$<br>$0$ | Max 2 |
| **Structural Context ($S_{\text{Structure}}$)** | `is_binding_site == True`<br>`is_pore_region == True`<br>`is_interface == True`<br>`is_tm_core == True` (if not pore)<br>`False` or `'unknown'` | $+4$<br>$+4$<br>$+2$<br>$+1$<br>$0$ | **Max 4 (Hard Cap)** |
| **CADD Phred ($S_{\text{CADD}}$)** | $\text{CADD} \ge 30$<br>$25 \le \text{CADD} < 30$<br>$20 \le \text{CADD} < 25$<br>$\text{CADD} < 20$ or Missing | $+3$<br>$+2$<br>$+1$<br>$0$ | Max 3 |
| **PolyPhen-2 ($S_{\text{PolyPhen}}$)** | `probably_damaging`<br>`possibly_damaging`<br>`benign` or Missing | $+2$<br>$+1$<br>$0$ | Max 2 |
| **SIFT ($S_{\text{SIFT}}$)** | `deleterious`<br>`tolerated` or Missing | $+1$<br>$0$ | Max 1 |
| **Total Possible** | | | **Max 12** |

### 3.2 Score Categorization
Variants are mapped into 4 qualitative risk categories:
- **High:** $\text{Score} \ge 10$
- **Medium:** $6 \le \text{Score} \le 9$
- **Low:** $3 \le \text{Score} \le 5$
- **Very Low:** $\text{Score} < 3$

> **Methodological Caveat:** Because in silico algorithms (CADD, PolyPhen-2, SIFT) correlate with each other, the score is a heuristic prioritization aid, **not** a calibrated posterior probability of pathogenicity.

---

## 4. Three-State Spatial Structural Semantics & Proximity Rules

### 4.1 Explicit Three-State Representation
To eliminate the classic bioinformatics bug where unmapped residues are erroneously classified as `False` (absent), RAREMISS enforces an explicit three-state object representation:
1. `True`: Feature assessed and verified present within 5.0 Å distance threshold.
2. `False`: Residue structurally resolved and mapped, but located $> 5.0$ Å away from feature.
3. `'unknown'`: Feature unassessable (residue is unresolved in the electron density / crystallographic model, or required ligand is absent).

### 4.2 Proximity & Proximity Calculation
- **Binding Site (`is_binding_site`):** Heavy atoms of the target residue are within $5.0$ Å of any hetero-atom (excluding water molecules and crystallographic ions unless biologically relevant).
- **Subunit Interface (`is_interface`):** Heavy atoms of the target residue on the evaluated chain are within $5.0$ Å of any heavy atom belonging to adjacent protein chains in the multimer.
- **Transmembrane Pore & Core:** Defined by UniProt topological domain annotations mapped onto structurally resolved regions.

### 4.3 Missing Structural Coverage
Residues that fall outside the PDB coordinate model (e.g., disordered loops, N/C-terminal tails) are flagged as `spatially_unresolved = True`. All their spatial feature flags are set to `'unknown'`. They receive $0$ structural points, and are strictly excluded from contingency analysis denominators.

---

## 5. Statistical Framework for Structural Feature Enrichment

### 5.1 Internal Score-Consistency Analysis
Enrichment tests (Module 7) assess whether variants categorized as **High Priority** ($\text{Score} \ge 10$) are statistically enriched for structural features compared to non-High variants.

Because structural features contribute directly to the priority score, **this analysis is inherently an internal score-consistency check** demonstrating that the scoring engine concentrates structural features in high tiers as designed. It must **not** be presented as independent clinical validation.

### 5.2 Contingency Matrix Formulation
For each feature, variants with `'unknown'` status are completely excluded:

| Group | Feature Present (`True`) | Feature Absent (`False`) |
| :--- | :---: | :---: |
| **High Priority** ($\text{Score} \ge 10$) | $A$ | $B$ |
| **Not High Priority** ($\text{Score} < 10$) | $C$ | $D$ |

- Odds Ratio: $\text{OR} = \frac{A \cdot D}{B \cdot C}$
- $p$-value: Calculated via two-sided Fisher's Exact Test (`scipy.stats.fisher_exact`).
- 95% Confidence Interval: Calculated via exact unconditional/conditional odds ratio bounds.
- Minimum Cell Rule: If $\min(A, B, C, D) < 3$, asymptotic assumptions are violated and the test status is marked `'skipped'`.

### 5.3 Multiple-Testing Correction (Benjamini–Hochberg FDR)
Across all $m$ successfully computed hypothesis tests within an analysis:
1. Sort raw $p$-values: $p_{(1)} \le p_{(2)} \le \dots \le p_{(m)}$.
2. Compute adjusted values:
   $$p^{\text{adj}}_{(i)} = \min\left(1.0, \frac{p_{(i)} \cdot m}{i}\right)$$
3. Enforce step-up monotonicity:
   $$p^{\text{adj}}_{(i)} = \min\left(p^{\text{adj}}_{(i)}, p^{\text{adj}}_{(i+1)}\right)$$
4. Strictly enforce $p^{\text{adj}}_{(i)} \ge p_{(i)}$ and preserve full floating-point precision (avoid rounding collapse to 0.0).

---

## 6. Independent Clinical Validation Principles

To achieve scientific validity suitable for peer-reviewed release, any external diagnostic benchmarking must adhere to:
1. **Curated Clinical Labels:** Golden standards from ClinVar (Pathogenic/Likely Pathogenic vs Benign/Likely Benign with $\ge 2$ submitter stars).
2. **Leakage Awareness:** Acknowledgment that CADD, SIFT, and PolyPhen-2 were trained on historical clinical variants.
3. **No Invented Metrics:** If a gene possesses insufficient benign or pathogenic controls, the analysis must report data scarcity rather than imputing artificial performance metrics.
