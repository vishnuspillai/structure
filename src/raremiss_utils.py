"""
raremiss_utils.py — Shared utilities for RAREMISS pipeline.

Three-state spatial annotation semantics
-----------------------------------------
Spatial feature columns (is_binding_site, is_interface, is_tm_core,
is_pore_region, is_transmembrane) use an explicit three-state representation:

    True   = feature assessed and POSITIVE  (annotated as present)
    False  = feature assessed and NEGATIVE  (annotated as absent)
    'unknown' = feature NOT assessable      (residue spatially unresolved,
                                             or no ligand available)

CSV round-trip hazard
---------------------
pandas writes True → "True", False → "False", 'unknown' → "unknown".
On reload, pandas infers the column as object dtype with strings "True",
"False", and "unknown".  Comparisons like col == True then fail silently
because "True" (str) != True (bool).

Use  parse_feature_col(series)  to normalise any spatial feature column
after loading a CSV.  The output series has dtype object with values:
    True (bool), False (bool), or 'unknown' (str).

Allele-frequency threshold semantics
-------------------------------------
Default threshold: 0.01
Ultra-rare setting: 0.001 (must be passed explicitly).
Missing AF (NaN): treated as UNKNOWN — variants are NOT automatically
assumed rare.  The handling policy must be documented per-analysis step.
"""

from __future__ import annotations

import os
import sys
import json
import warnings
import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Three-state feature column parsing
# ---------------------------------------------------------------------------

_TRUTHY_STRINGS  = {"true", "1", "1.0", "yes"}
_FALSY_STRINGS   = {"false", "0", "0.0", "no"}
_UNKNOWN_STRINGS = {"unknown", "nan", "none", ""}


def parse_feature_value(value) -> "bool | str":
    """
    Normalise a single spatial-feature value to True, False, or 'unknown'.

    Accepts:
        bool True/False, int 0/1, float 0.0/1.0, str variants of the above,
        the string 'unknown', NaN.

    Raises ValueError for any other value.
    """
    if value is True or value is False:
        return value
    if isinstance(value, float):
        if np.isnan(value):
            return 'unknown'
        if value == 1.0:
            return True
        if value == 0.0:
            return False
        raise ValueError(f"Unexpected float feature value: {value!r}")
    if isinstance(value, int):
        if value == 1:
            return True
        if value == 0:
            return False
        raise ValueError(f"Unexpected int feature value: {value!r}")
    if isinstance(value, str):
        s = value.strip().lower()
        if s in _TRUTHY_STRINGS:
            return True
        if s in _FALSY_STRINGS:
            return False
        if s in _UNKNOWN_STRINGS:
            return 'unknown'
        raise ValueError(f"Unexpected string feature value: {value!r}")
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return 'unknown'
    raise ValueError(f"Unexpected feature value type {type(value)!r}: {value!r}")


def parse_feature_col(series: pd.Series, col_name: str = "") -> pd.Series:
    """
    Normalise a spatial feature column after CSV load.

    Returns an object-dtype Series with values True (bool), False (bool),
    or 'unknown' (str).  Any unrecognised value raises ValueError with
    the column name and offending value.
    """
    def _apply(v):
        try:
            return parse_feature_value(v)
        except ValueError as exc:
            label = f"column '{col_name}': " if col_name else ""
            raise ValueError(f"Malformed annotation in {label}{exc}") from exc

    return series.apply(_apply)


def parse_feature_cols(df: pd.DataFrame, cols: list[str] | None = None) -> pd.DataFrame:
    """
    In-place normalise spatial annotation columns in a DataFrame.

    Parameters
    ----------
    df   : DataFrame to modify (a copy is NOT made — caller is responsible).
    cols : list of column names to normalise.  Defaults to the standard
           spatial feature set.

    Columns absent from the DataFrame are silently skipped.
    """
    if cols is None:
        cols = [
            "is_binding_site", "is_interface", "is_tm_core",
            "is_pore_region", "is_transmembrane",
        ]
    df = df.copy()
    for col in cols:
        if col in df.columns:
            df[col] = parse_feature_col(df[col], col_name=col)
    return df


def load_spatial_csv(path: str, cols: list[str] | None = None) -> pd.DataFrame:
    """
    Load a spatial-annotated CSV and normalise all feature columns.

    This is the preferred way to load any file written by Module 3 or
    Module 4 so that downstream comparisons like col == True are safe.
    """
    df = pd.read_csv(path)
    return parse_feature_cols(df, cols)


# ---------------------------------------------------------------------------
# Feature-column summary helpers
# ---------------------------------------------------------------------------

def feature_counts(df: pd.DataFrame, col: str) -> dict:
    """
    Return {'positive': int, 'negative': int, 'unknown': int, 'total': int}
    for a spatial feature column.  Raises KeyError if col absent.
    """
    s = df[col]
    pos     = int((s == True).sum())   # noqa: E712 — explicit identity check
    neg     = int((s == False).sum())  # noqa: E712
    unk     = int((s == 'unknown').sum())
    total   = len(s)
    if pos + neg + unk != total:
        warnings.warn(
            f"Feature column '{col}': {total - pos - neg - unk} values are "
            "neither True, False, nor 'unknown'. Run parse_feature_cols first."
        )
    return {"positive": pos, "negative": neg, "unknown": unk, "total": total}


# ---------------------------------------------------------------------------
# Allele-frequency threshold validation
# ---------------------------------------------------------------------------

DEFAULT_AF_THRESHOLD = 0.01
_ULTRA_RARE_AF      = 0.001


def validate_af_threshold(value, source: str = "") -> float:
    """
    Validate and return an allele-frequency threshold.

    Accepts float or int in (0, 1].  Raises ValueError for out-of-range or
    non-numeric values.  A value of 0.0 is explicitly rejected with an
    explanatory message (it would filter all known variants).

    Parameters
    ----------
    value  : candidate threshold value (from config, CLI arg, API call, …)
    source : description of where the value came from (for error messages)
    """
    try:
        v = float(value)
    except (TypeError, ValueError):
        tag = f" (from {source})" if source else ""
        raise ValueError(
            f"AF threshold{tag} must be a number, got {value!r}"
        )
    if v == 0.0:
        tag = f" (from {source})" if source else ""
        raise ValueError(
            f"AF threshold{tag} is exactly 0.0.  This would filter out every "
            "variant with a known frequency.  Use a small positive value "
            f"(e.g. {DEFAULT_AF_THRESHOLD}) or {_ULTRA_RARE_AF} for ultra-rare "
            "analysis."
        )
    if not (0.0 < v <= 1.0):
        tag = f" (from {source})" if source else ""
        raise ValueError(
            f"AF threshold{tag} must be in (0, 1], got {v}"
        )
    return v


def get_af_threshold(config: dict, source: str = "config") -> float:
    """
    Extract and validate the AF threshold from a config dict.

    Returns DEFAULT_AF_THRESHOLD (0.01) when key is absent.
    Raises ValueError for invalid values.
    """
    raw = config.get("af_threshold", DEFAULT_AF_THRESHOLD)
    return validate_af_threshold(raw, source=source)


# ---------------------------------------------------------------------------
# Scoring helpers (shared with Module 4)
# ---------------------------------------------------------------------------

def compute_score(row: "pd.Series | dict") -> int:
    """
    Compute the RAREMISS heuristic priority score for a single variant.

    This is the single authoritative implementation of the scoring formula.
    Module 4 (prioritization.py) imports and calls this function.

    Score components
    ----------------
    Allele frequency (max 2 pts):
        AF < 1e-5              → +2
        1e-5 ≤ AF < 1e-4      → +1
        AF ≥ 1e-4 or missing  → +0
        (Missing AF is NOT treated as evidence of rarity.)

    Structural context (max 4 pts, capped):
        is_binding_site == True  → +4
        is_pore_region  == True  → +4  (also counts as pore for is_tm_core)
        is_pore_core    == True  → +4  (alias; either column accepted)
        is_interface    == True  → +2
        is_tm_core      == True  (and not pore) → +1
        'unknown' or False values → +0

    CADD phred (max 3 pts):
        CADD ≥ 30              → +3
        25 ≤ CADD < 30        → +2
        20 ≤ CADD < 25        → +1

    PolyPhen-2 prediction (max 2 pts):
        probably_damaging      → +2
        possibly_damaging      → +1

    SIFT prediction (max 1 pt):
        deleterious            → +1

    Maximum possible score: 12
    Priority categories: High ≥ 10, Medium 6–9, Low 3–5, Very Low < 3.

    Note: CADD, SIFT and PolyPhen derive from correlated evidence sources.
    The combined score should be interpreted as a heuristic ranking aid,
    not a calibrated probability of pathogenicity.
    """
    score = 0

    # --- Allele frequency ---
    af = row.get('gnomAD_AF', np.nan)
    try:
        af = float(af) if af is not None else np.nan
    except (ValueError, TypeError):
        af = np.nan
    if pd.notna(af):
        if af < 1e-5:
            score += 2
        elif af < 1e-4:
            score += 1

    # --- Structural context ---
    struct_score = 0
    if row.get('is_binding_site') is True:
        struct_score += 4
    is_pore = (row.get('is_pore_core') is True) or (row.get('is_pore_region') is True)
    if is_pore:
        struct_score += 4
    if row.get('is_interface') is True:
        struct_score += 2
    if (row.get('is_tm_core') is True) and not is_pore:
        struct_score += 1
    struct_score = min(struct_score, 4)
    score += struct_score

    # --- CADD ---
    cadd = row.get('cadd_phred', np.nan)
    try:
        cadd = float(cadd) if cadd is not None else np.nan
    except (ValueError, TypeError):
        cadd = np.nan
    if pd.notna(cadd):
        if cadd >= 30:
            score += 3
        elif cadd >= 25:
            score += 2
        elif cadd >= 20:
            score += 1

    # --- PolyPhen-2 ---
    polyphen = str(row.get('polyphen_pred', '') or '').lower()
    if 'probably_damaging' in polyphen:
        score += 2
    elif 'possibly_damaging' in polyphen:
        score += 1

    # --- SIFT ---
    sift = str(row.get('sift_pred', '') or '').lower()
    if 'deleterious' in sift:
        score += 1

    return score


def get_category(score: int) -> str:
    """Map a numeric priority score to a category label."""
    if score >= 10:
        return 'High'
    if score >= 6:
        return 'Medium'
    if score >= 3:
        return 'Low'
    return 'Very Low'


# ---------------------------------------------------------------------------
# Multiple-testing correction: Benjamini–Hochberg FDR
# ---------------------------------------------------------------------------

def benjamini_hochberg(p_values: list[float] | np.ndarray) -> list[float]:
    """
    Benjamini–Hochberg False Discovery Rate (FDR) step-up multiple testing correction.

    Parameters
    ----------
    p_values : list or 1D array of unadjusted p-values (each in [0, 1]).

    Returns
    -------
    list of adjusted p-values matching input order, preserving floating precision.
    Guarantees monotonicity: p_adj >= p_raw and step-up cumulative minimum.
    """
    if len(p_values) == 0:
        return []
    m = len(p_values)
    raw_sorted = sorted(enumerate(p_values), key=lambda x: float(x[1]))
    adj = [None] * m
    for rank, (orig_idx, p) in enumerate(raw_sorted):
        p_flt = float(p)
        adj[rank] = min(p_flt * m / (rank + 1), 1.0)
    for i in range(m - 2, -1, -1):
        adj[i] = min(adj[i], adj[i + 1])
    result = [None] * m
    for rank, (orig_idx, p) in enumerate(raw_sorted):
        result[orig_idx] = float(max(float(p), adj[rank]))
    return result


# ---------------------------------------------------------------------------
# End-to-end run provenance & metadata
# ---------------------------------------------------------------------------

def record_run_metadata(
    output_dir: str,
    config: dict,
    gene_symbol: str,
    extra_info: dict | None = None
) -> dict:
    """
    Record execution provenance metadata for reproducibility.
    Captures:
      - Timestamp (UTC ISO 8601)
      - Git commit hash & branch
      - Effective configuration & validated AF threshold
      - Gene, transcript, translation, UniProt identity
      - Genome build / assembly
      - Structure & chain ID
      - Software environment / library versions
      - External data sources & endpoint URLs
    """
    import subprocess
    import datetime

    commit = "unknown"
    branch = "unknown"
    try:
        r_c = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, timeout=5)
        if r_c.returncode == 0:
            commit = r_c.stdout.strip()
        r_b = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"], capture_output=True, text=True, timeout=5)
        if r_b.returncode == 0:
            branch = r_b.stdout.strip()
    except Exception:
        pass

    af_val = get_af_threshold(config, source="run_metadata")
    metadata = {
        "gene_symbol": gene_symbol.upper(),
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "git_commit": commit,
        "git_branch": branch,
        "genome_assembly": config.get("genome_build", "GRCh38"),
        "species": config.get("species", "homo_sapiens"),
        "transcript_id": config.get("transcript_id", "unspecified"),
        "translation_id": config.get("translation_id", "unspecified"),
        "uniprot_id": config.get("uniprot_id", "unspecified"),
        "structure_id": str(config.get("structure_id", "unspecified")).upper(),
        "chain_id": config.get("chain_id", "A"),
        "af_threshold_effective": af_val,
        "af_threshold_default": DEFAULT_AF_THRESHOLD,
        "consequence_filter": config.get("consequence_filter", "missense_variant"),
        "software_versions": {
            "python": sys.version.split()[0],
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
        "external_data_sources": {
            "ensembl_rest": "https://rest.ensembl.org",
            "uniprot_rest": "https://rest.uniprot.org",
            "pdbe_sifts": "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot",
            "rcsb_pdb": "https://files.rcsb.org/download",
            "ncbi_entrez": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        }
    }
    if extra_info:
        metadata.update(extra_info)

    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, f"{gene_symbol.lower()}_run_metadata.json")
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=4)
    return metadata

