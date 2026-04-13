import subprocess
import sys
import argparse
import yaml
import os

def run_pipeline():
    parser = argparse.ArgumentParser(description="Run the RAREMISS pipeline.")
    parser.add_argument('--gene', type=str, help="Gene symbol to process")
    parser.add_argument('--structure', type=str, help="PDB Structure ID")
    parser.add_argument('--af', type=float, help="Allele frequency threshold")
    args, unknown = parser.parse_known_args()

    # Update config if args provided
    if args.gene or args.structure or args.af:
        config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config", "parameters.yaml")
        if os.path.exists(config_path):
            with open(config_path, "r") as f:
                config = yaml.safe_load(f) or {}
        else:
            config = {}
            
        if args.gene:
            config['gene_symbol'] = args.gene
        if args.structure:
            config['structure_id'] = args.structure
        if args.af:
            config['af_threshold'] = args.af
            
        with open(config_path, "w") as f:
            yaml.safe_dump(config, f)
        print(f"Updated config/parameters.yaml with: {args}")
    steps = [
        ("Module 1: Extract and filter rare baseline variants", "src/module1_variant_mining/variant_mining.py"),
        ("Module 2: Phase A/B/C Annotations", "src/module2_annotation/annotation_layer.py"),
        ("Module 3: Structural Mapping", "src/module3_spatial/spatial_annotation.py"),
        ("Module 4: Calculate Priority Scores and Categorize", "src/module4_prioritization/prioritization.py"),
        ("Module 5: Isolate Mechanistic Distances", "src/module5_mechanistic/mechanistic_metrics.py"),
        ("Module 6: Query NCBI ClinVar and PubMed", "src/module6_clinical/clinical_review.py"),
        ("Module 7: Validate Feature Enrichment", "src/module7_validation/structural_enrichment.py"),
        ("Module 7: Validate Feature Enrichment (CI)", "src/module7_validation/structural_enrichment_ci.py")
    ]

    print("=== Starting RAREMISS Pipeline Execution =========================================")
    
    for step_num, (description, script_path) in enumerate(steps, 1):
        print(f"\n=> [{step_num}/{len(steps)}] Running {description}...")
        print(f"   Executing script: {script_path}")
        
        try:
            # Run the script and stream output to the console
            subprocess.run([sys.executable, script_path], check=True)
            print(f"   [ SUCCESS ]")
        except subprocess.CalledProcessError as e:
            print(f"\n   [ ERROR ] Step {step_num} failed with exit code {e.returncode}.")
            print("================================================================================\n")
            print("Pipeline execution aborted. Fix the underlying error and retry.")
            sys.exit(1)
            
    print("\n================================================================================")
    print("=== Pipeline Execution Completed Successfully! ===")

if __name__ == "__main__":
    run_pipeline()
