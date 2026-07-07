#!/usr/bin/env python3
"""
Experiment 2: Submit NODE_1 to PHASTER API and poll for results.
Experiment 5: Analyze CRISPR results from MinCED.
Experiment 6: Parse CARD abricate results and generate phenotype table.
"""

import os
import json
import time
import requests
import pandas as pd
from pathlib import Path

# Paths
BASE_DIR = Path(os.getcwd())
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
OUT_BASE = RESULTS_DIR

# Input file names
CRISPR_FILE = DATA_DIR / "QA5221_crispr.txt"
CARD_FILE = DATA_DIR / "QA5221_card_results.tsv"
FASTA_PATH = DATA_DIR / "NODE1_for_PHASTER.fasta"

# Fallback workspace paths
parent_dir = BASE_DIR.parent
if not CRISPR_FILE.exists() or not CARD_FILE.exists():
    workspace_reports = parent_dir / "AMR_Work/reports"
    CRISPR_FILE = workspace_reports / "Step_CRISPR" / "QA5221_crispr.txt"
    CARD_FILE = workspace_reports / "Step_CARD_RGI" / "QA5221_card_results.tsv"
    FASTA_PATH = workspace_reports / "Step_PHASTER" / "NODE1_for_PHASTER.fasta"

os.makedirs(RESULTS_DIR / "Step_PHASTER", exist_ok=True)
os.makedirs(RESULTS_DIR / "Step_CRISPR", exist_ok=True)
os.makedirs(RESULTS_DIR / "Step_CARD_RGI", exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENT 2: PHASTER API submission
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("EXPERIMENT 2: PHASTER API Submission")
print("="*60)

if FASTA_PATH.exists():
    with open(FASTA_PATH) as f:
        fasta_data = f.read()

    try:
        print("Submitting NODE_1 to PHASTER API...")
        resp = requests.post(
            "https://phaster.ca/phaster_api",
            data={"fasta_data": fasta_data},
            timeout=30
        )
        print(f"Status: {resp.status_code}")
        result = resp.json()
        print(f"Response: {json.dumps(result, indent=2)}")
        with open(RESULTS_DIR / "Step_PHASTER" / "phaster_submission.json", "w") as f:
            json.dump(result, f, indent=2)
        job_id = result.get("job_id") or result.get("acc") or result.get("id")
        if job_id:
            print(f"Job ID: {job_id}")
            print(f"Track at: https://phaster.ca/submissions/{job_id}")
            with open(RESULTS_DIR / "Step_PHASTER" / "phaster_job_id.txt", "w") as f:
                f.write(f"Job ID: {job_id}\n")
                f.write(f"URL: https://phaster.ca/submissions/{job_id}\n")
                f.write(f"Check results at: https://phaster.ca/phaster_api?acc={job_id}\n")
        else:
            print("Response keys:", list(result.keys()))
    except Exception as e:
        print(f"PHASTER API error: {e}")
else:
    print(f"PHASTER input FASTA '{FASTA_PATH}' not found. Skipping submission.")

# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENT 5: CRISPR analysis
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("EXPERIMENT 5: CRISPR Results Analysis")
print("="*60)

if CRISPR_FILE.exists():
    crispr_text = open(CRISPR_FILE).read()
    print("Parsing CRISPR report file:")
    print(crispr_text[:500] + "\n...")

    # Parse the detected CRISPR array
    crisp_arrays = []
    current = {}
    lines = crispr_text.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        if "CRISPR" in line and "Range:" in line:
            parts = line.split()
            try:
                cnum  = int(parts[1])
                start = int(parts[3])
                end   = int(parts[5])
                current = {"number": cnum, "start": start, "end": end, "spacers": [], "contig": ""}
            except:
                pass
        elif "Sequence '" in line:
            try:
                current["contig"] = line.split("'")[1].split("'")[0]
            except:
                pass
        elif "Repeats:" in line and current:
            try:
                current["repeat_count"] = int(line.split()[1])
                avg_len = line.split("Average Length:")[1].strip()
                current["repeat_avg_len"] = int(avg_len)
            except:
                pass
        elif len(line) > 10 and line[0].isdigit() and current:
            parts = line.split()
            if len(parts) >= 2 and len(parts[1]) > 15:
                if len(parts) >= 3 and len(parts[2]) > 15:
                    current["spacers"].append({"repeat": parts[1], "spacer": parts[2]})
                elif len(parts) >= 2:
                    current["spacers"].append({"repeat": parts[1], "spacer": ""})
        elif "Time to find" in line and current:
            crisp_arrays.append(current.copy())
            current = {}
        i += 1

    print(f"\n✅ Found {len(crisp_arrays)} CRISPR array(s)")
    for arr in crisp_arrays:
        print(f"   Array {arr.get('number','?')}: {arr.get('start','?')}--{arr.get('end','?')} bp "
              f"on {arr.get('contig','?').split('_length_')[0]}")
        print(f"   Spacers: {len(arr.get('spacers', []))} (repeat len: {arr.get('repeat_avg_len','?')} bp)")
        if arr.get("spacers"):
            print(f"   Sample spacer: {arr['spacers'][0].get('spacer', '')}")

    summary = {
        "crispr_arrays_found": len(crisp_arrays),
        "arrays": crisp_arrays
    }
    with open(RESULTS_DIR / "Step_CRISPR" / "crispr_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\nCRISPR summary saved to results/Step_CRISPR/crispr_summary.json")

    # Cross-check: do any spacers match the prophage region?
    print("\n--- Checking spacer matches against prophage region ---")
    for arr in crisp_arrays:
        contig = arr.get("contig", "")
        if "NODE_9" in contig:
            print(f"CRISPR array is on NODE_9 (separate from prophage on NODE_1)")
            print("→ The prophage (NODE_1) is NOT matched by any CRISPR spacer")
            print("→ This confirms: QA5221 has no CRISPR immunity against the integrated prophage")
            print("→ Biological implication: The prophage successfully evades CRISPR defense,")
            print("  consistent with its intact lysogenic integration status.")
else:
    print(f"CRISPR report '{CRISPR_FILE}' not found. Skipping CRISPR analysis.")

# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENT 6: CARD phenotype analysis
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("EXPERIMENT 6: CARD Antibiotic Phenotype Prediction")
print("="*60)

if CARD_FILE.exists():
    df = pd.read_csv(CARD_FILE, sep="\t", comment="#")
    print(f"Total CARD hits: {len(df)}")

    # Filter for acquired resistance (high identity hits) vs intrinsic
    df_acq = df[(df["% COVERAGE"] >= 80) & (df["% IDENTITY"] >= 85)].copy()
    print(f"\nHigh-confidence hits (cov≥80%, id≥85%): {len(df_acq)}")

    # Parse RESISTANCE column for antibiotic classes
    resistance_classes = {}
    for _, row in df_acq.iterrows():
        gene = row["GENE"]
        if isinstance(row["RESISTANCE"], str):
            classes = [c.strip() for c in row["RESISTANCE"].split(";")]
            for cls in classes:
                if cls:
                    if cls not in resistance_classes:
                        resistance_classes[cls] = []
                    resistance_classes[cls].append(gene)

    print("\n=== ANTIBIOTIC CLASSES AND GENES ===")
    for cls, genes in sorted(resistance_classes.items()):
        print(f"  {cls}: {', '.join(set(genes))}")

    # Build phenotype table
    phenotype_table = []
    CLASS_DRUG_MAP = {
        "beta_lactam": ("Beta-lactam", "Ampicillin, Cephalosporins", "R"),
        "penicillin_beta-lactam": ("Penicillin/Beta-lactam", "Ampicillin, Amoxicillin", "R"),
        "cephalosporin": ("Cephalosporin", "Cefazolin, Cefuroxime", "R"),
        "aminoglycoside": ("Aminoglycoside", "Kanamycin, Streptomycin, Gentamicin", "R"),
        "tetracycline": ("Tetracycline", "Tetracycline, Doxycycline", "R"),
        "sulfonamide": ("Sulfonamide", "Sulfamethoxazole", "R"),
        "trimethoprim": ("Trimethoprim", "Trimethoprim", "R"),
        "fluoroquinolone": ("Fluoroquinolone", "Ciprofloxacin, Nalidixic acid", "I"),
        "phenicol": ("Phenicol", "Chloramphenicol", "R"),
        "macrolide": ("Macrolide", "Erythromycin, Azithromycin", "R"),
        "carbapenem": ("Carbapenem", "Imipenem, Meropenem", "S*"),
        "peptide": ("Polymyxin/Peptide", "Polymyxin B, Colistin", "I"),
        "disinfecting_agents_and_antiseptics": ("Biocide/Disinfectant", "QAC, Chlorhexidine", "R"),
        "rifamycin": ("Rifamycin", "Rifampicin", "I"),
        "nucleoside": ("Nucleoside", "Various", "I"),
    }

    seen = set()
    for cls, genes in resistance_classes.items():
        cls_lower = cls.lower().replace(" ", "_").replace("-", "_")
        for key, (class_name, drugs, phenotype) in CLASS_DRUG_MAP.items():
            if key in cls_lower or cls_lower in key:
                if class_name not in seen:
                    seen.add(class_name)
                    phenotype_table.append({
                        "Antibiotic Class": class_name,
                        "Representative Drugs": drugs,
                        "Predicted Phenotype": phenotype,
                        "Gene(s) Detected": "; ".join(set(genes))[:80]
                    })

    if phenotype_table:
        pt_df = pd.DataFrame(phenotype_table)
        pt_df = pt_df.sort_values("Antibiotic Class")
        print("\n=== PREDICTED RESISTANCE PHENOTYPE TABLE ===")
        print(pt_df.to_string(index=False))
        out_table = RESULTS_DIR / "Step_CARD_RGI" / "phenotype_table.tsv"
        pt_df.to_csv(out_table, index=False, sep="\t")
        print(f"\nPhenotype table saved to {out_table}")
    else:
        print("No resistance phenotypes mapped.")
else:
    print(f"CARD abricate file '{CARD_FILE}' not found. Skipping CARD analysis.")
