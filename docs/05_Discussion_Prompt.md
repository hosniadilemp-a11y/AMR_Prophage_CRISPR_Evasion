================================================================================
LLM PROMPT FOR CURATING RELEVANT LITERATURE REFERENCES FOR PAPER 2
================================================================================

Instructions: Copy and paste the entire text below into a high-capacity LLM (e.g. Claude 3, GPT-4o, or Gemini Pro) to find and format supporting scientific literature for our paper.

------------------------------ COPY START ------------------------------

You are an expert scientific literature search specialist and bioinformatician. For our paper (Paper 2) titled: "Genomic Epidemiology, Mobilome Characterisation, and Intact Prophage Analysis of the Clinical Multidrug-Resistant Escherichia coli ST354 Isolate QA5221", please perform a literature search and provide a curated list of relevant reference papers related to our findings.

For each reference paper you find, you must provide:
1. The Title, Authors, and Year of the publication.
2. What specific finding or aspect of our work it supports, or what is similar to our work (e.g. pangenomics, resistome networks, codon amelioration, intact prophage boundaries, CRISPR evasion, Siphoviridae taxonomic markers, or MD simulations).
3. Why it is important to our work (detailed explanation of the scientific connection, methodology, or biological context).

Rules for references selection:
- Prioritize recent papers (published between 2021 and 2026).
- Focus on papers using similar pipelines: pangenomics (Panaroo/Roary), mobilome/resistome (CARD, Jaccard networks), codon usage bias, prophage curation, CRISPR-Cas targeting, structural homology (ESMFold/Foldseek), and explicit-solvent Molecular Dynamics simulations (OpenMM/GROMACS) for validating fold stability.
- At the end of your response, provide the complete, copy-pasteable BibTeX entries for all cited references. Every BibTeX entry MUST include a valid, clickable DOI field (e.g. `doi = {10.1016/j.jmb.2023.168124}`). Do not use placeholder DOIs.

Here are the key results sections and findings of our paper that you need to find references for:

1. BACTERIAL EPIDEMIOLOGY AND PANGENOMICS:
   - Our findings: Comparative pangenomic analysis of E. coli QA5221 placed within an expanded clonal cohort (400 ST354 genomes) and species-wide cohort (800 genomes) shows a highly fluid, open pangenome structure (Heap's law exponent alpha = 0.5716). The core genome constitutes 8.27% of the total 37,366 gene families, while the accessory cloud genome represents the largest partition (85.72%).
   - We need references on: Pangenome openness in ExPEC lineages, Heap's law modeling in E. coli, and pangenomics resolving clonal complexes vs. species-wide diversity.

2. INTEGRATED MOBILOME AND RESISTOME MAPPING:
   - Our findings: Mapped 12 acquired resistance and biocide tolerance genes conferring an 8-class resistance profile. Resistance determinants are physically linked on horizontal mobile backbones (CALIN arrays and a Class 1 integron flanked by active IS3 elements). Jaccard co-occurrence network analysis confirmed tight linkage of sul1, sul3, dfrA7, and biocides (qacL, qacE\Delta1).
   - We need references on: Biocide tolerance driving co-selection of AMR cassettes in hospitals/farms, and network-based analysis of resistance gene co-occurrence.

3. CODON ADAPTATION AND GC/GC3 AMELIORATION SIGNATURES:
   - Our findings: Significant discrepancies in GC content and GC3 codon bias between core (51.86% GC, 55.90% GC3), accessory (48.49% GC), and singleton (46.78% GC) partitions. Kruskal-Wallis H-testing (p < 10^-60) and post-hoc Mann-Whitney U-tests confirm GC depression in singletons, providing sequence-independent mathematical signatures of horizontal gene transfer (HGT) from low-GC donor reservoirs (e.g. Firmicutes, phages).
   - We need references on: Codon adaptation index (CAI) and GC3 amelioration signatures as mathematical indicators of horizontal gene transfer (HGT) in Proteobacteria.

4. INTACT PROPHAGE BOUNDARIES AND CHIMERIC MODULES:
   - Our findings: Identified a complete, 42.3-kb prophage on NODE_1 (coordinates 56,422--98,747 bp) containing all five predicted functional modules (Lysogeny, Replication, Packaging, Structural, Lysis). DLP12 synteny comparison confirms that replication/regulation modules are shared, but the tail and capsid structural core represent a newly acquired, complete horizontal block with homology to HK97.
   - We need references on: Chimeric prophage evolution in Enterobacteriaceae (lambdoid/HK97 hybrids), and mosaic structural module inheritance in Caudoviricetes.

5. CRISPR-CAS IMMUNE EVASION (SPACER-NAIVE GAP):
   - Our findings: Detected a functional Type I-E CRISPR array on NODE_9 (16 repeats, 15 spacers) with a complete cas operon. However, BLASTn spacer targeting analysis revealed only short seed fragments (max 17 bp), far below the biological threshold for interference (20 bp seed match). This naive integration state explains the stable maintenance of the prophage without host immune clearance.
   - We need references on: CRISPR-Cas spacer-naive gaps enabling the stable integration and persistence of temperate bacteriophages and plasmids.

6. PHAGE STRUCTURAL TAXONOMY (STRUCTURAL SINGLETONS):
   - Our findings: Dual BLASTp and Foldseek structural homology search of structural singletons (terminase AC4NUP_00460 and major capsid AC4NUP_00555) classified the prophage as Siphoviridae (order Caudoviricetes) with a restricted E. coli host range (all sequence-level homologs confined to Escherichia).
   - We need references on: Sequence-independent structural taxonomy (using Foldseek/TM-align) to classify highly divergent phage structural components and define host-range constraints.

7. LINEAGE-WIDE PROPHAGE CARRIAGE PREVALENCE (FIGURE 17 COHORTS):
   - Our findings: Prevalence analysis (Figure 17) shows that the prophage cargo (enterohemolysin, rcbA) and switches behave as private singletons in the local 32-genome cohort (3.12% prevalence), but are present in exactly 9 genomes (2.25%) in the 400-genome ST354 cohort, restricted to a single correlated sub-clade. Capsid and tail assembly genes are present in only 2 genomes (0.50% prevalence), indicating highly localized clade-specific prophage mobilization.
   - We need references on: Clade-specific carriage of prophage elements, and transition of private singletons in small cohorts into low-frequency accessory markers in expanded cohorts.

8. MOLECULAR DYNAMICS (MD) FOLD STABILITY VALIDATION:
   - Our findings: Explicit-solvent MD simulations validate the biophysical feasibility of key candidate folds under physiological conditions:
     (i) Enterohemolysin host cargo (AC4NUP_00305) simulated embedded in POPC/POPG (3:1) membrane bilayer achieved structural stability (C\alpha RMSD 5.65 \pm 1.04 \AA, Rg 27.18 \pm 0.20 \AA).
     (ii) Antirepressor regulator (AC4NUP_00420) simulated in aqueous solution achieved rapid equilibration (C\alpha RMSD 2.24 \pm 0.18 \AA, Rg 19.45 \pm 0.12 \AA).
   - We need references on: Explicit membrane bilayer molecular dynamics simulations (in POPC/POPG/POPE) to evaluate the structural integrity and stability of membrane-associated proteins, and aqueous MD validation of transcriptional regulators.

------------------------------- COPY END -------------------------------
