# CHRNA7 Variant Structural Prioritization Framework

A fully reproducible, modular bioinformatics pipeline to prioritize rare missense variants in the human CHRNA7 gene using public datasets.

## Project Structure
- `data/`: Raw and processed datasets
- `results/`: Tables and figures
- `src/`: Modular codebase
  - `module1_variant_mining`
  - `module2_annotation`
  - `module3_structural_mapping`
  - `module4_prioritization`
- `logs/`: Application logs
- `config/`: Configuration files

## Setup
1. Create a virtual environment and install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
2. Execute modules sequentially from the `src/` directory.

## Current Module
**Module 1: Variant Mining**
Retrieves missense variants mapped to the CHRNA7 canonical transcript, filtering by gnomAD Allele Frequency < 0.001.
