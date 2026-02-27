import subprocess
import sys

def run_pipeline():
    steps = [
        ("Module 1: Extract and filter rare baseline variants", "src/module1_variant_mining/variant_mining.py"),
        ("Module 2: Phase A/B/C Annotations", "src/module2_annotation/annotation_layer.py"),
        ("Module 3: Structural Mapping against PDB 7KOX", "src/module3_spatial/spatial_annotation.py"),
        ("Module 4: Calculate Priority Scores and Categorize", "src/module4_prioritization/prioritization.py"),
        ("Module 5: Isolate Mechanistic Distances", "src/module5_mechanistic/mechanistic_metrics.py"),
        ("Module 6: Query NCBI ClinVar and PubMed", "src/module6_clinical/clinical_review.py"),
        ("Module 7: Validate Feature Enrichment", "src/module7_validation/structural_enrichment.py"),
        ("Module 7: Validate Feature Enrichment (CI)", "src/module7_validation/structural_enrichment_ci.py")
    ]

    print("=== Starting CHRNA7 Pipeline Execution =========================================")
    
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
