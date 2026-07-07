# Detailed Roadmap: Proposed Improvements for Paper 2 (Group B: The Bacteriophage Mobilome Story)

This document provides a comprehensive analysis of the improvements proposed for Paper 2 (the characterization of a novel active prophage region in the ST354 lineage, focusing on the 5 phage-associated candidates: antirepressor _00084, terminase _00091, tail _00107, baseplate _00114, and structural capsid _00109). For each item, we address the data status, implementation protocol, time/compute requirements, and peer-review impact.

---

### 1. Pangenome Database Expansion (400 Genomes)
* **Status of Data:** **In Progress (90% complete).** FastTree phase is active.
* **How to Implement:** Once the 400-genome phylogeny tree is complete, query the presence and conservation of this prophage locus across all ST354 genomes.
* **Time Required:** ~2.0 hours.
* **Compute Power Needed:** Active usage of 36 threads, ~35 GB RAM.
* **Impact on Paper:** Evaluates the evolutionary conservation of the prophage. It shows whether this prophage is a highly mobile element newly acquired by QA5221, or a stable lineage-wide prophage region that has evolved with the ST354 clonal complex.

---

### 2. Prophage Boundary and Attachment Junction (attL/attR) Characterization
* **Status of Data:** **No data currently.**
* **How to Implement:** Run a dedicated prophage prediction server (e.g. PHASTER or PhageFinder) on the contigs of QA5221 (specifically NODE_1, which contains the 5 candidates). Locate the exact coordinates of the integrase, lysin, tail, and capsid elements. Identify the flanking direct repeats representing the attachment sites (attL and attR).
* **Time Required:** Fast (~30 minutes of tool run).
* **Compute Power Needed:** Low (runs on PHASTER web server or locally via Docker).
* **Impact on Paper:** Crucial for genomic annotation. High-quality virology and genomics journals require precise mapping of prophage boundaries to prove the structural limits of the mobilome element.

---

### 3. Cargo Gene Identification and Annotation
* **Status of Data:** **Available.** We have the full annotation files.
* **How to Implement:** Extract all open reading frames (ORFs) located between the predicted attL and attR attachment sites. Screen any "hypothetical proteins" in this region (cargo genes) against Swiss-Prot and Pfam. Assess if any cargo genes encode virulence factors, toxin-antitoxin systems, or metabolic enzymes that could provide selective advantage to E. coli.
* **Time Required:** Fast (~2 hours of sequence analysis).
* **Compute Power Needed:** Low.
* **Impact on Paper:** Reveals the biological theme. Characterizing the cargo genes shows exactly how this prophage acts as an active vector for genetic plasticity, carrying novel traits into the ST354 host.

---

### 4. Comparative Synteny Analysis against Reference Phages
* **Status of Data:** **Available.**
* **How to Implement:** Retrieve complete phage genomes showing blast similarity. Execute clinker or Easyfig to align and visualize synteny between our ST354 prophage and reference bacteriophages, showing conservation of the structural tail, baseplate, and capsid blocks.
* **Time Required:** Fast (~2 hours).
* **Compute Power Needed:** Low.
* **Impact on Paper:** Places the prophage in taxonomic context. Showing whether it belongs to the Myoviridae, Siphoviridae, or Podoviridae families based on structural homology adds significant quality to the mobilome narrative.

---

### 5. Assessment of Prophage Viability (Active vs Decayed Pseudoprophage)
* **Status of Data:** **Can be analyzed immediately.**
* **How to Implement:** Evaluate the completeness of the prophage structural machinery. Check for the presence of essential phage modules: (1) Lysogeny/regulation (integrase, repressors), (2) Replication, (3) Packaging (terminase), (4) Structural (capsid, tail, baseplate), and (5) Lysis (endolysin, holin).
* **Time Required:** Minimal (~2 hours of database screening).
* **Compute Power Needed:** Low.
* **Impact on Paper:** Determines the biological viability. If the prophage contains all essential structural and lysis genes, it represents a potentially active, inducible prophage capable of lysogenic cell lysis, rather than a silent pseudoprophage remnant.