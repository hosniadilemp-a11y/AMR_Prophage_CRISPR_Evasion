# Execution Plan: Bacteriophage Mobilome Improvements (Paper 2)

This plan outlines the steps, code, and validation process to execute three key improvements for Paper 2 (Group B: The Bacteriophage Mobilome Story) using the draft assembly of *E. coli* QA5221.

## GPU Requirement Assessment
None of the proposed steps require a GPU. They are entirely CPU-bound bioinformatics analyses (parsing GFF/GBK, running sequence lookups, and visual plotting) and can be executed quickly on standard CPU cores.

---

## Proposed Changes & Tasks

### Task 1: Cargo Gene Identification and Annotation
* **Goal:** Extract and annotate all genes located within the prophage region on NODE_1 (approx. 70,000 to 98,000 bp) to identify any functional cargo (toxins, virulence, or metabolic markers).
* **Steps:**
  1. Write a Python script `scripts/extract_prophage_cargo.py` to parse `results/Step3_Annotation/prokka_out/Ecoli_isolate.gff`.
  2. Extract all genes between the coordinates 70,000 and 98,000 bp on NODE_1.
  3. Filter the hypothetical proteins in this region and run a local database search (or web BLAST) to identify their remote homologs and domains.
  4. Catalog the results in a new report file `reports/Step12_Prophage_Cargo_Report.txt`.

### Task 2: Prophage Viability Assessment
* **Goal:** Inventory the prophage region to determine if it encodes all modules required for active lysogeny and cell lysis.
* **Steps:**
  1. Analyze the functional annotations of all genes in the prophage region.
  2. Check for the presence of:
     * *Lysogeny / Regulation:* Integrase (`int` / `ybcK`), repressors, division inhibitors (`ybcS`).
     * *DNA Replication / Modification:* Replicative helicase loader (`dnaC` / `ybcN`), nucleases (`ybcO`).
     * *Packaging:* Phage terminase (`_00091`).
     * *Structural:* Portal, capsid (`_00109`), tail assembly (`_00107`, `_00108`), baseplate J-like (`_00114`).
     * *Lysis:* Endolysin, holin, or tail lysozymes.
  3. Formulate a diagnostic conclusion: Is the prophage a complete, potentially inducible phage, or a decayed pseudoprophage? Add this analysis to `reports/Step12_Prophage_Cargo_Report.txt`.

### Task 3: Bacteriophage Synteny Mapping (clinker)
* **Goal:** Align and visualize synteny between the ST354 prophage region and reference phages (like E. coli bacteriophage DLP12) using `clinker`.
* **Steps:**
  1. Extract the genbank/GFF coordinates of the prophage region from `Ecoli_isolate.gff` and write a sub-genbank file `prophage_QA5221.gbk`.
  2. Download reference genbank files for closely related phages (e.g. NC_001799 for DLP12 or related lambdoid phages).
  3. Run `clinker` in `amr_env` to generate an interactive HTML comparison plot:
     ```bash
     conda run -n amr_env clinker prophage_QA5221.gbk dlp12_ref.gbk --plot clinker_prophage_synteny.html
     ```
  4. Save the results and export the figures to `AMR_Work/manuscript/plots/` for inclusion in the paper.

---

## Verification Plan

### Automated Verification
* Verify that the script `scripts/extract_prophage_cargo.py` runs with exit code 0 and generates `reports/Step12_Prophage_Cargo_Report.txt` containing the annotated cargo list.
* Verify that the compiled `clinker` HTML is created successfully and contains aligned blocks showing percentage sequence identity.

### Manual Verification
* Inspect the prophage cargo table and check if the annotated products make biological sense.
* Open the `clinker_prophage_synteny.html` plot in the web browser and verify that the 5 focus genes and cargo genes align with reference sequences.
