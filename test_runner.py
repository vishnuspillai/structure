"""
Terminal-based validation runner for generalized RAREMISS pipeline.
Runs specified test cases sequentially, capturing logs and validating outputs.
"""
import subprocess
import sys
import os
import yaml
import shutil
import pandas as pd
import json
import time

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.join(ROOT_DIR, "config", "parameters.yaml")
DATA_DIR = os.path.join(ROOT_DIR, "data", "processed")
RAW_DIR = os.path.join(ROOT_DIR, "data", "raw")

def set_config(gene_symbol, structure_id, af_threshold=0.001, species="homo_sapiens"):
    """Write test parameters to config."""
    config = {
        "af_threshold": af_threshold,
        "consequence_filter": "missense_variant",
        "gene_symbol": gene_symbol,
        "genome_build": "GRCh38",
        "log_dir": "logs",
        "output_dir": "data/processed",
        "species": species,
        "structure_id": structure_id,
        "transcript_id": None,
        "translation_id": None,
        "uniprot_id": None,
    }
    with open(CONFIG_PATH, "w") as f:
        yaml.safe_dump(config, f)
    print(f"  Config set: gene={gene_symbol}, structure={structure_id}")

def clear_processed():
    """Remove all files in data/processed/."""
    if os.path.exists(DATA_DIR):
        for f in os.listdir(DATA_DIR):
            fp = os.path.join(DATA_DIR, f)
            if os.path.isfile(fp):
                os.remove(fp)
    print("  Cleared data/processed/")

def run_pipeline():
    """Execute the full pipeline and capture output. Returns (success, log_text)."""
    result = subprocess.run(
        [sys.executable, "run_pipeline.py"],
        capture_output=True,
        text=True,
        cwd=ROOT_DIR,
        encoding="utf-8",
        errors="replace",
    )
    full_log = result.stdout + "\n" + result.stderr
    return result.returncode == 0, full_log, result.returncode

def validate_outputs(gene_symbol):
    """Validate the output files for a given gene. Returns (pass, report_dict)."""
    gs = gene_symbol.lower()
    report = {}
    warnings = []

    # 1. Check ranked variants exist
    ranked_csv = os.path.join(DATA_DIR, f"{gs}_ranked_variants.csv")
    if not os.path.exists(ranked_csv):
        return False, {"error": f"Missing {gs}_ranked_variants.csv"}, ["FATAL: no ranked CSV"]

    df = pd.read_csv(ranked_csv)
    report["variant_count"] = len(df)
    if len(df) == 0:
        warnings.append("variant_count is 0")

    # 2. Check rare variant filtering
    if "gnomAD_AF" in df.columns:
        above_threshold = (df["gnomAD_AF"] > 0.001).sum()
        report["variants_above_af_threshold"] = int(above_threshold)
        if above_threshold > 0:
            warnings.append(f"{above_threshold} variants above AF threshold in ranked output")

    # 3. VEP annotations
    for col in ["cadd_phred", "sift_pred", "polyphen_pred"]:
        if col in df.columns:
            non_null = df[col].notna().sum()
            report[f"{col}_non_null"] = int(non_null)
        else:
            report[f"{col}_non_null"] = 0
            warnings.append(f"Missing column: {col}")

    # 4. Structural mapping coverage
    spatial_csv = os.path.join(DATA_DIR, f"{gs}_missense_spatial_annotated.csv")
    if os.path.exists(spatial_csv):
        sdf = pd.read_csv(spatial_csv)
        if "spatially_unresolved" in sdf.columns:
            resolved = (~sdf["spatially_unresolved"]).sum()
            total = len(sdf)
            coverage = (resolved / total * 100) if total > 0 else 0
            report["mapping_coverage_pct"] = round(coverage, 2)
            if coverage < 60:
                warnings.append(f"Mapping coverage {coverage:.1f}% < 60%")
        else:
            warnings.append("No spatially_unresolved column in spatial CSV")
    else:
        report["mapping_coverage_pct"] = 0
        warnings.append("Missing spatial annotated CSV")

    # 5. Ligand detection (binding site column)
    if "is_binding_site" in df.columns:
        report["binding_site_count"] = int(df["is_binding_site"].sum())
    else:
        report["binding_site_count"] = 0
        warnings.append("No is_binding_site column")

    # 6. Scoring distribution
    if "priority_score" in df.columns:
        report["score_mean"] = round(df["priority_score"].mean(), 2)
        report["score_max"] = int(df["priority_score"].max())
        if df["priority_score"].sum() == 0:
            warnings.append("All priority scores are 0")
    else:
        warnings.append("No priority_score column")

    # 7. High-priority count
    if "priority_category" in df.columns:
        report["high_priority_count"] = int((df["priority_category"] == "High").sum())
    else:
        report["high_priority_count"] = 0

    # 8. Enrichment ran
    ci_file = os.path.join(DATA_DIR, f"{gs}_structural_ci.txt")
    report["enrichment_ci_exists"] = os.path.exists(ci_file)
    if not os.path.exists(ci_file):
        warnings.append("Enrichment CI file missing")

    report["warnings"] = warnings
    passed = len([w for w in warnings if w.startswith("FATAL")]) == 0
    return passed, report, warnings


def run_test(gene_symbol, structure_id, test_num):
    """Full test cycle for one gene/structure pair."""
    print(f"\n{'='*60}")
    print(f"TEST {test_num}: {gene_symbol} / {structure_id}")
    print(f"{'='*60}")

    # 1. Clear & configure
    clear_processed()
    set_config(gene_symbol, structure_id)

    # 2. Run pipeline
    print("  Running pipeline...")
    start = time.time()
    success, log_text, exit_code = run_pipeline()
    elapsed = time.time() - start
    print(f"  Pipeline finished in {elapsed:.1f}s with exit code {exit_code}")

    # Save log
    log_path = os.path.join(ROOT_DIR, "logs", f"test_{gene_symbol.lower()}.log")
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "w", encoding="utf-8", errors="replace") as f:
        f.write(log_text)
    print(f"  Full log saved to {log_path}")

    if not success:
        # Print last 40 lines of log for debugging
        lines = log_text.strip().split("\n")
        print(f"\n  PIPELINE FAILED! Last 40 lines of output:")
        for line in lines[-40:]:
            print(f"    {line}")
        return {"gene": gene_symbol, "structure": structure_id, "status": "CRASHED",
                "exit_code": exit_code, "log_tail": "\n".join(lines[-40:])}

    # 3. Validate
    print("  Validating outputs...")
    passed, report, warnings = validate_outputs(gene_symbol)
    report["gene"] = gene_symbol
    report["structure"] = structure_id
    report["status"] = "PASS" if passed else "FAIL"
    report["elapsed_seconds"] = round(elapsed, 1)

    print(f"  Status: {report['status']}")
    print(f"  Variant count: {report.get('variant_count', 'N/A')}")
    print(f"  Mapping coverage: {report.get('mapping_coverage_pct', 'N/A')}%")
    print(f"  High priority: {report.get('high_priority_count', 'N/A')}")
    if warnings:
        print(f"  Warnings:")
        for w in warnings:
            print(f"    - {w}")

    return report


if __name__ == "__main__":
    tests = [
        # ("TP53", "6M0J"),
        # ("KCNQ2", "7CR0"),
        ("BRCA1", "1T15"),
    ]

    results = []
    for i, (gene, struct) in enumerate(tests, 1):
        result = run_test(gene, struct, i)
        results.append(result)

    # Final summary
    print(f"\n{'='*60}")
    print("FINAL SUMMARY")
    print(f"{'='*60}")
    for r in results:
        status = r.get("status", "UNKNOWN")
        gene = r.get("gene")
        struct = r.get("structure")
        vc = r.get("variant_count", "N/A")
        mc = r.get("mapping_coverage_pct", "N/A")
        hp = r.get("high_priority_count", "N/A")
        warns = len(r.get("warnings", []))
        print(f"  {gene}/{struct}: {status} | variants={vc} | coverage={mc}% | high_priority={hp} | warnings={warns}")

    # Save summary JSON
    summary_path = os.path.join(ROOT_DIR, "logs", "validation_summary.json")
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSummary saved to {summary_path}")
