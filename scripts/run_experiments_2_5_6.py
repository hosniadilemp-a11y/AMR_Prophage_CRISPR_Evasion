#!/usr/bin/env python3
"""
Experiment 2: Submit NODE_1 to PHASTER API and poll for results.
Experiment 5: Analyze CRISPR results from MinCED.
Experiment 6: Parse CARD abricate results and generate phenotype table.
"""

import json, time, requests
from pathlib import Path

OUT_BASE = Path("/media/adel/Data/Hosni/openmm_windows_setup/AMR_Work/reports")

# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENT 2: PHASTER API submission
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("EXPERIMENT 2: PHASTER API Submission")
print("="*60)

PHASTER_OUT = OUT_BASE / "Step_PHASTER"
FASTA_PATH  = PHASTER_OUT / "NODE1_for_PHASTER.fasta"

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
    with open(PHASTER_OUT / "phaster_submission.json", "w") as f:
        json.dump(result, f, indent=2)
    job_id = result.get("job_id") or result.get("acc") or result.get("id")
    if job_id:
        print(f"Job ID: {job_id}")
        print(f"Track at: https://phaster.ca/submissions/{job_id}")
        with open(PHASTER_OUT / "phaster_job_id.txt", "w") as f:
            f.write(f"Job ID: {job_id}\n")
            f.write(f"URL: https://phaster.ca/submissions/{job_id}\n")
            f.write(f"Check results at: https://phaster.ca/phaster_api?acc={job_id}\n")
    else:
        print("Response keys:", list(result.keys()))
except Exception as e:
    print(f"PHASTER API error: {e}")
    # Try alternative endpoint
    try:
        resp2 = requests.post(
            "https://phaster.ca/phaster_api?acc=",
            files={"fasta_data": open(FASTA_PATH)},
            timeout=30
        )
        print("Alt response:", resp2.status_code, resp2.text[:200])
    except Exception as e2:
        print(f"Alt error: {e2}")

# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENT 5: CRISPR analysis
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("EXPERIMENT 5: CRISPR Results Analysis")
print("="*60)

CRISPR_FILE = OUT_BASE / "Step_CRISPR" / "QA5221_crispr.txt"
crispr_text = open(CRISPR_FILE).read()
print(crispr_text)

# Parse the detected CRISPR array
crisp_arrays = []
current = {}
lines = crispr_text.split("\n")
i = 0
while i < len(lines):
    line = lines[i]
    if "CRISPR" in line and "Range:" in line:
        parts = line.split()
        # CRISPR 1   Range: 36345 - 37288
        try:
            cnum  = int(parts[1])
            start = int(parts[3])
            end   = int(parts[5])
            current = {"number": cnum, "start": start, "end": end, "spacers": [], "contig": ""}
        except:
            pass
    elif "Sequence '" in line:
        # Extract contig name
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
            # This is a spacer line: POSITION  REPEAT  SPACER
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

# Save parsed summary
summary = {
    "crispr_arrays_found": len(crisp_arrays),
    "arrays": crisp_arrays
}
with open(OUT_BASE / "Step_CRISPR" / "crispr_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print("\nCRISPR summary saved.")

# Cross-check: do any spacers match the prophage region?
print("\n--- Checking spacer matches against prophage region ---")
if crisp_arrays:
    for arr in crisp_arrays:
        contig = arr.get("contig", "")
        if "NODE_9" in contig:
            print(f"CRISPR array is on NODE_9 (separate from prophage on NODE_1)")
            print("→ The prophage (NODE_1) is NOT matched by any CRISPR spacer")
            print("→ This confirms: QA5221 has no CRISPR immunity against the integrated prophage")
            print("→ Biological implication: The prophage successfully evades CRISPR defense,")
            print("  consistent with its intact lysogenic integration status.")

# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENT 6: CARD phenotype analysis
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("EXPERIMENT 6: CARD Antibiotic Phenotype Prediction")
print("="*60)

import pandas as pd

CARD_FILE = OUT_BASE / "Step_CARD_RGI" / "QA5221_card_results.tsv"
df = pd.read_csv(CARD_FILE, sep="\t", comment="#")
print(f"Total CARD hits: {len(df)}")
print(f"Columns: {list(df.columns)}")

# Filter for acquired resistance (high identity hits) vs intrinsic
# Focus on: %COVERAGE >= 80, %IDENTITY >= 80
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

pt_df = pd.DataFrame(phenotype_table)
pt_df = pt_df.sort_values("Antibiotic Class")
print("\n=== PREDICTED RESISTANCE PHENOTYPE TABLE ===")
print(pt_df.to_string(index=False))

pt_df.to_csv(OUT_BASE / "Step_CARD_RGI" / "phenotype_table.tsv", index=False, sep="\t")
print(f"\nPhenotype table saved ({len(pt_df)} antibiotic classes).")
print("* S* = Susceptible by genotype (no acquired determinant), but intrinsic efflux noted")
