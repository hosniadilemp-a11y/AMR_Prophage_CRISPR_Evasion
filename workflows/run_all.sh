#!/bin/bash
# =============================================================================
#          Paper 2 Reproducibility Pipeline Execution Orchestrator
# =============================================================================
set -e

echo "=== [Paper 2 Pipeline Started] ==="

# Create standard output directories if they don't exist
mkdir -p results logs figures data

# Step 0: Download Reference Bacteriophage DLP12
echo "[Step 0/6] Downloading reference bacteriophage DLP12..."
python scripts/download_reference_phage.py > logs/00_download_ref_phage.log 2>&1 || {
    echo "[WARNING] Reference download failed (usually due to network). Skipping download..."
}

# Step 1: Scan for Prophage Insertion Junctions (attL/attR)
echo "[Step 1/6] Running sliding-window attL/attR scanner..."
python scripts/01_find_att_sites.py > logs/01_att_sites_scan.log 2>&1

# Step 2: Extract & Annotate Prophage Cargo Genes
echo "[Step 2/6] Extracting prophage cargo coordinates..."
python scripts/02_extract_prophage_cargo.py > logs/02_extract_prophage_cargo.log 2>&1

# Step 3: Run AMR & QAC Co-Selection Jaccard Networks
echo "[Step 3/6] Running Jaccard resistome co-occurrence network..."
python scripts/03_run_amr_cooccurrence.py > logs/03_amr_cooccurrence.log 2>&1

# Step 4: Profile CRISPR Evasion & CARD Phenotypes
echo "[Step 4/6] Profiling CRISPR arrays and CARD phenotype metrics..."
python scripts/04_run_experiments_2_5_6.py > logs/04_experiments_2_5_6.log 2>&1

# Step 5: Perform MD Trajectory Analysis (Ehly_61 & Antirepressor)
echo "[Step 5/6] Calculating RMSD/Rg/RMSF for Molecular Dynamics simulations..."
python scripts/05_analyze_paper2_mds.py > logs/05_trajectory_analysis.log 2>&1

# Step 6: Generate Publication-Quality Figures & Plots
echo "[Step 6/6] Generating paper figures..."
python scripts/06_generate_paper2_plots.py > logs/06_generate_plots.log 2>&1

echo "=== [Paper 2 Pipeline Completed Successfully] ==="
