PAPER 2 — FINDINGS SUMMARY
=============================================================================
Title: Computational Characterization of a Novel Genomically Intact Prophage
       Mobilome and Genomic Profiling of Virulence and Resistance Hotspots in
       the Extraintestinal Pathogenic Escherichia coli ST354 Lineage
Journal: FEMS Microbiology Letters, 2026, Volume 373, Issue 1
Isolate: QA5221 (ExPEC ST354) | NCBI Accession: PRJNA1481519 (AC4NUP_XXXXX)
=============================================================================

──────────────────────────────────────────────────────────────────────────────
FINDING 1 ▶ OPEN PANGENOME WITH EXTREME PLASTICITY
──────────────────────────────────────────────────────────────────────────────
• Pangenome constructed from 32 local ST354 genomes + 800-genome global cohort
• Heap's Law exponent α = 0.5716 → Open pangenome (α < 1)
• Total gene families: 37,366 across 800 genomes
  - Core genome (present in ≥ 99%): 3,088 genes = 8.27%
  - Shell genome (accessory, 15–99%): 2,228 genes = 5.97%
  - Cloud genome (singletons/rare, < 15%): 32,050 genes = 85.72%
• QA5221-specific singletons: highly unique, private gene pool
• New gene discovery rate has NOT plateaued → lineage still expanding

──────────────────────────────────────────────────────────────────────────────
FINDING 2 ▶ MULTIDRUG RESISTANCE RESISTOME (12 ACQUIRED GENES, 8 CLASSES)
──────────────────────────────────────────────────────────────────────────────
• 12 acquired resistance/biocide genes detected by ABRicate + AMRFinderPlus
• Classes RESISTANT (R): beta-lactam, aminoglycoside, sulfonamide,
  trimethoprim, phenicol, macrolide, tetracycline, biocide/QAC
• Classes INTERMEDIATE (I) via intrinsic efflux: fluoroquinolone,
  polymyxin, aminocoumarin, rifamycin, glycylcycline
• Carbapenem susceptibility: RETAINED (no acquired carbapenemase genes)
• Key genes: blaTEM-1B, aph(3')-Ia, aadA1, aadA2, sul1, sul3, dfrA7,
  cmlA1, estX, tet(A), qacL, qacEΔ1
• Physical clustering within a class 1 integron on NODE_37 (with IS3
  transposons and CALIN cassette array)
• Jaccard co-occurrence network: strong physical linkage between
  sulfonamide/trimethoprim genes and QAC biocide determinants
  → Biocide co-selection maintains the entire integron cassette array

──────────────────────────────────────────────────────────────────────────────
FINDING 3 ▶ COMPOSITIONAL MOSAICISM — CODON USAGE EVIDENCE OF RECENT HGT
──────────────────────────────────────────────────────────────────────────────
• Singleton GC content: 46.78% vs. core genome: 51.86% (p < 0.001)
• Singleton GC3 content: 44.56% vs. core genome: 55.90% (p < 0.001)
• This GC/GC3 depression is consistent with genome amelioration theory:
  recently acquired genes retain the low-GC signature of their donor
• Codon Adaptation Index (CAI) and Nc (effective number of codons)
  calculations confirm recent horizontal gene transfer of singletons
• Provides a sequence-independent method to identify novel foreign genes

──────────────────────────────────────────────────────────────────────────────
FINDING 4 ▶ INTACT PROPHAGE: 42.3 kb LAMBDOID SIPHOVIRIDAE
──────────────────────────────────────────────────────────────────────────────
• Complete, 42,325-bp prophage resolved on NODE_1 (coordinates 56,422–98,747 bp)
• Integration junctions: attL (AAAAAAACCGCC) and attR (AAAACAACCGCC)
  → Single nucleotide transition confirms site-specific integration
• Automated tool PhiSpy MISSED the boundaries due to host-adapted GC composition
  → Manual curation was required
• Five functional modules confirmed:
  (i)  Lysogeny & Regulation (ydaT transcriptional regulator)
  (ii) DNA Replication & Repair (rcbA — double-strand break reduction)
  (iii) Packaging — large terminase subunit (AC4NUP_00460)
  (iv) Head & Capsid — major capsid protein HK97-like (AC4NUP_00555)
  (v)  Virulence Cargo — enterohemolysin cytolysin (AC4NUP_00305)
• DLP12 synteny analysis: lysogeny/replication modules share lambda-like origin;
  capsid/tail structural core is a newly acquired HK97-homologous block
  → Chimeric mosaic prophage (lambda backbone + HK97 structural module)

──────────────────────────────────────────────────────────────────────────────
FINDING 5 ▶ CRISPR IMMUNE GAP — PROPHAGE EVADES HOST IMMUNITY PASSIVELY
──────────────────────────────────────────────────────────────────────────────
• One CRISPR array (CRISPR-1) on NODE_9 (36,345–37,288 bp)
• 16 direct repeats, 15 spacers (32 bp average)
• Complete Type I-E cas operon: cas3, casA, casB, casC, casD, casE, cas1, cas2
  → Structurally INTACT Cascade system
• BLASTn of all 15 spacers vs. prophage sequence:
  Maximum seed match = 17 bp (< 20 bp biological minimum for Cas3 targeting)
• Result: CRISPR spacer-naive gap → prophage persists without immune clearance
• No anti-CRISPR (Acr) proteins found in AcrDB searches
• A canonical LexA binding box (5'-CTGTATGGGTAAATACAG-3') was identified
  42 bp upstream of the terminase → SOS response governs lytic induction

──────────────────────────────────────────────────────────────────────────────
FINDING 6 ▶ PHAGE TAXONOMY — SEQUENCE-INDEPENDENT STRUCTURAL CLASSIFICATION
──────────────────────────────────────────────────────────────────────────────
• BLASTp — Terminase (AC4NUP_00460, 249 aa): 100% hits to E. coli only
  (99.2–100% identity, E-value < 10⁻¹⁷⁹)
• BLASTp — Major Capsid (AC4NUP_00555, 344 aa): 100% hits to E. coli only
  (99.4–100% identity, E-value = 0)
• Foldseek structural homologs (afdb50):
  - Terminase: Klebsiella aerogenes (89.8%), Salmonella enterica (87.5%),
    Citrobacter portucalensis (71.4%) → Siphoviridae / Lambdoid
  - Capsid: Serratia sp. (70.4%), Jejubacter calystegiae (48.4%),
    Cronobacter turicensis (41.5%) → Enterobacterial Siphovirus
• TAXONOMIC CONCLUSION: QA5221 prophage = lambdoid Siphoviridae (Caudoviricetes)
  with an E. coli-restricted host range

──────────────────────────────────────────────────────────────────────────────
FINDING 7 ▶ PROPHAGE CANDIDATE PREVALENCE — COHORT-SCALE EPIDEMIOLOGY
──────────────────────────────────────────────────────────────────────────────
• Local 32-genome cohort:
  - All prophage genes (cargo + structural): present in 1/32 = 3.12%
  → Private to QA5221
• Expanded 400-genome ST354 clonal cohort:
  - Enterohemolysin (AC4NUP_00305), rcbA, antirepressor: 9/400 = 2.25%
  - Structural capsid/tail genes: 2/400 = 0.50%
  → Restricted to one phylogenetically correlated sub-clade
• Global 800-genome E. coli cohort:
  - Antirepressor switch: 23/800 = 2.88%
  - Terminase: 15/800 = 1.88%
  - Structural capsid/tail: 4/800 = 0.50%
• Pattern: Cargo genes (enterohemolysin) are more widely retained than
  structural phage machinery → consistent with "moron" cargo gene dynamics
  (hosts selectively retain virulence cargo; redundant structural machinery decays)

──────────────────────────────────────────────────────────────────────────────
FINDING 8 ▶ ESMFOLD STRUCTURAL MODELING — HIGH-CONFIDENCE PREDICTIONS
──────────────────────────────────────────────────────────────────────────────
• Enterohemolysin (AC4NUP_00305, 350 aa):
  - Mean pLDDT = 88.4 (very high confidence)
  - 234/350 residues (66.9%) > pLDDT 90
  - 300/350 residues (85.7%) > pLDDT 70
  - Top Foldseek hit: E. coli putative enterohemolysin (98.2% identity)
  - Additional hits: Salmonella enterica (47.1%), Proteus sp. (45.5%)
  - Confirmed cytolysin toxin / pore-forming fold

──────────────────────────────────────────────────────────────────────────────
FINDING 9 ▶ MOLECULAR DYNAMICS — BIOPHYSICAL VALIDATION OF 3D FOLDS
──────────────────────────────────────────────────────────────────────────────
• Enterohemolysin (AC4NUP_00305) — Membrane simulation (100 ns):
  - Bilayer: POPC/POPG (3:1) — models anionic bacterial inner membrane
  - Backbone RMSD: 5.65 ± 1.04 Å (stable equilibrium)
  - Radius of gyration Rg: 27.18 ± 0.20 Å (compact, stable)
  - RMSF: < 2.0 Å in transmembrane helical core
  → CONFIRMS: stable transmembrane topology = active cytolysin fold

• Antirepressor switch (AC4NUP_00420, 203 aa) — Aqueous simulation (100 ns):
  - Backbone RMSD: 2.24 ± 0.18 Å (rapid equilibration)
  - Radius of gyration Rg: 19.45 ± 0.12 Å
  - RMSF: < 1.5 Å in structured helical domains
  → CONFIRMS: stable folded transcriptional regulator conformation

• OVERALL: Both ESMFold models correspond to dynamically stable folds,
  validating them as biologically functional proteins rather than artifacts

=============================================================================
OVERALL CONCLUSION
=============================================================================
The ST354 ExPEC clinical isolate QA5221 carries a novel, genomically intact
lambdoid Siphoviridae prophage that:
  1. Evades CRISPR immunity through a spacer-naive gap (not anti-CRISPR)
  2. Harbors virulence cargo (enterohemolysin cytolysin + rcbA) confirmed by
     structural homology and MD simulation
  3. Is rare at population scale (< 3% prevalence) but exclusively E. coli-
     restricted — suggesting recent, clonal mobilization
  4. Co-exists with an 8-class MDR resistome physically co-localized within
     a class 1 integron — subject to biocide co-selection pressure
  5. Has a chimeric mosaic architecture (lambda regulatory core + HK97 structural
     head) — demonstrating active recombinatorial phage evolution
=============================================================================
