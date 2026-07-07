#!/usr/bin/env python3
"""
Experiment 1: ESMFold structural prediction for enterohemolysin (KNGPFPPJ_00061)
Uses HuggingFace transformers ESMFold model.
"""

import torch
import os, time
from pathlib import Path

OUT_DIR = Path("/media/adel/Data/Hosni/openmm_windows_setup/AMR_Work/reports/Step_ESMFold")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SEQUENCE = ("MNELTQQENINSNVAVFSPQSLAAIQTFSQVMASGMATVPEHLRGNPSDCMAITMQAMQWQMNPYAVAQK"
            "TFVVNGVLGYEAQLVNAVISTRGPLTGRIEYDWFGPWEKIIGKFEIRKNDKGKEYR"
            "VPGWKLADENGIGVRVQATLRGESKPRVLELLLAQARTRN"
            "STLWADDPRQQLAYLALKRWARLYCP"
            "EVILGVYTRDELDEPQEKIINPVQEHKNTSACRAERETTIIIEQDAGENWISAFRERIIEQAQSTGETTALE"
            "RQEVEDHKNTLGALYTELKGKVVQRHHRLNAIARIIEKMINDLPSS"
            "GDPEAEQKFIALENTLNAARPHLGELYEAYKIMTLTDMKPEYIGS")

# Use the canonical 350 aa sequence from the annotation
SEQUENCE = "MNELTQQENINSNVAVFSPQSLAAIQTFSQVMASGMATVPEHLRGNPSDCMAITMQAMQWQMNPYAVAQKTFVVNGVLGYEAQLVNAVISTRGPLTGRIEYDWFGPWEKIIGKFEIRKNDKGKEYRVPGWKLADENGIGVRVQATLRGESKPRVLELLLAQARTRN STLWADDPRQQLAYLALKRWARLYCP EVILGVYTRDELDEPQEKIINPVQEHKNTSACRAERETTIIIEQDAGENWISAFRERIEQAQSTGETTLRQEVEDHKNTLGALYTELKGKVVQRHHRLNAIARIIEKMINDLPSSGDPEAEQKFIALENTLNAARPHLGELYEAYKIMTLTDMKPEYIGS"
SEQUENCE = SEQUENCE.replace(" ", "")

print(f"Sequence length: {len(SEQUENCE)} aa")
print(f"Sequence: {SEQUENCE[:50]}...")

try:
    from transformers import EsmForProteinFolding, AutoTokenizer
    from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
    from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
    import torch

    print("Loading ESMFold tokenizer...")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")

    print("Loading ESMFold model (this may take 2-3 min first time)...")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")
    model = model.to(device)
    model.esm = model.esm.half()
    model.eval()

    print("Running ESMFold prediction...")
    t0 = time.time()
    tokenized = tokenizer([SEQUENCE], return_tensors="pt", add_special_tokens=False)
    tokenized = {k: v.to(device) for k, v in tokenized.items()}

    with torch.no_grad():
        output = model(**tokenized)

    print(f"Prediction completed in {time.time()-t0:.1f}s")

    # Extract pLDDT scores
    plddt = output["plddt"].squeeze(0).cpu().numpy()  # (L, 37) per-atom
    plddt_ca = plddt[:, 1]  # CA atom pLDDT
    mean_plddt = float(plddt_ca.mean())
    print(f"Mean pLDDT (CA): {mean_plddt:.2f}")

    # Save pLDDT scores
    import numpy as np
    np.save(OUT_DIR / "enterohemolysin_plddt.npy", plddt_ca)

    # Convert to PDB
    from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
    from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

    converted = model.output_to_pdb(output)
    pdb_str = converted[0]
    pdb_path = OUT_DIR / "enterohemolysin_00061_esmfold.pdb"
    with open(pdb_path, "w") as f:
        f.write(pdb_str)
    print(f"PDB saved: {pdb_path}")
    print(f"Mean pLDDT: {mean_plddt:.2f}")
    if mean_plddt > 70:
        print("✅ HIGH CONFIDENCE: pLDDT > 70 — reliable structural prediction")
    elif mean_plddt > 50:
        print("⚠️  MEDIUM CONFIDENCE: pLDDT 50-70")
    else:
        print("❌ LOW CONFIDENCE: pLDDT < 50")

except ImportError as e:
    print(f"Import error: {e}")
    print("Trying ESMFold API instead...")
    # Fallback: use ESMAtlas API
    import requests, json
    url = "https://esmatlas.com/api/predictStructure"
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    data = f"sequence={SEQUENCE}"
    print("Submitting to ESMAtlas API (this may take 2-3 min)...")
    t0 = time.time()
    resp = requests.post(url, headers=headers, data=data, timeout=300)
    print(f"Response status: {resp.status_code} (took {time.time()-t0:.1f}s)")
    if resp.status_code == 200:
        pdb_str = resp.text
        pdb_path = OUT_DIR / "enterohemolysin_00061_esmfold.pdb"
        with open(pdb_path, "w") as f:
            f.write(pdb_str)
        print(f"PDB saved: {pdb_path}")
        # Extract pLDDT from PDB B-factor column
        plddt_vals = []
        for line in pdb_str.split("\n"):
            if line.startswith("ATOM") and " CA " in line:
                try:
                    plddt_vals.append(float(line[60:66]))
                except:
                    pass
        if plddt_vals:
            import numpy as np
            arr = np.array(plddt_vals)
            np.save(OUT_DIR / "enterohemolysin_plddt.npy", arr)
            mean_plddt = float(arr.mean())
            print(f"Mean pLDDT (from PDB B-factors): {mean_plddt:.2f}")
    else:
        print(f"API error: {resp.text[:500]}")
