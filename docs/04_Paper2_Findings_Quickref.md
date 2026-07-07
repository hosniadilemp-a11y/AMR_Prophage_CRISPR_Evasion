PAPER 2 — KEY FINDINGS QUICK REFERENCE CARD
=============================================================================
Isolate: QA5221 (ExPEC ST354) | Journal: FEMS Microbiology Letters, 2026
=============================================================================

PANGENOME
  • Open pangenome (α = 0.5716); 37,366 gene families across 800 genomes
  • Core genome = 8.27% | Cloud (singletons) = 85.72% of total gene pool
  • Gene discovery has not plateaued → ST354 lineage is still expanding

RESISTOME
  • 12 acquired resistance/biocide genes spanning 8 antibiotic classes
  • Carbapenem susceptibility is RETAINED (no carbapenemase genes found)
  • All resistance genes cluster inside a class 1 integron on NODE_37
  • Biocide (QAC) genes physically co-selected with sulfonamide/trimethoprim genes

HORIZONTAL GENE TRANSFER SIGNATURES
  • Singleton GC content 5.08% lower than core genome (p < 0.001)
  • Singleton GC3 content 11.34% lower than core genome (p < 0.001)
  • Depression in GC/GC3 = molecular proof of recent, un-ameliorated HGT

PROPHAGE — STRUCTURE
  • Intact 42.3 kb lambdoid prophage on NODE_1 (attL / attR junction confirmed)
  • PhiSpy missed boundaries — manual curation was required
  • All 5 functional modules present: lysogeny, replication, packaging, capsid, cargo
  • Chimeric architecture: lambda-like regulatory core + HK97-like capsid block

PROPHAGE — CRISPR EVASION
  • One Type I-E CRISPR array (16 repeats, 15 spacers) on NODE_9
  • Maximum spacer-prophage match = 17 bp (below the ≥20 bp Cas3 threshold)
  • Prophage evades immunity PASSIVELY — no anti-CRISPR proteins present
  • LexA binding site 42 bp upstream of terminase → SOS-controlled lytic induction

PROPHAGE — PHAGE TAXONOMY
  • Terminase (AC4NUP_00460): 100% BLASTp hits to E. coli only
  • Capsid (AC4NUP_00555): 100% BLASTp hits to E. coli only
  • Foldseek structural taxonomy: Siphoviridae / Lambdoid (Enterobacteriaceae)
  • Conclusion: novel E. coli-restricted lambdoid Siphoviridae

PROPHAGE — POPULATION PREVALENCE
  • Local (32 genomes): 3.12% — private to QA5221
  • ST354 cohort (400 genomes): cargo genes 2.25%, structural genes 0.50%
  • Global (800 genomes): cargo genes 1.62–2.88%, structural genes 0.50%
  • Cargo genes more widely distributed than structural genes (moron dynamics)

ESMFOLD STRUCTURAL PREDICTION
  • Enterohemolysin (AC4NUP_00305): mean pLDDT = 88.4 → cytolysin toxin fold
  • 85.7% of residues exceed pLDDT 70 threshold → very high confidence
  • Top Foldseek hit: E. coli enterohemolysin (98.2% structural identity)

MOLECULAR DYNAMICS VALIDATION
  • Enterohemolysin (membrane sim, 100 ns): RMSD = 5.65 ± 1.04 Å — STABLE
  • Antirepressor switch (aqueous sim, 100 ns): RMSD = 2.24 ± 0.18 Å — STABLE
  • Both proteins maintain compact folds → biologically functional predictions

OVERALL MESSAGE
  • QA5221 ST354 carries a rare but dangerous intact prophage that passively
    evades host CRISPR and delivers a confirmed cytolysin virulence toxin,
    embedded within a multidrug-resistant mobile genetic element background.

=============================================================================
Files: paper2.pdf | supplementary.pdf | Finals/ | Submition/paper2/
Code:  https://github.com/hosniadilemp-a11y/AMR_Novel_Gene_Discovery
Data:  NCBI PRJNA1481519 (AC4NUP_XXXXX annotations)
=============================================================================
