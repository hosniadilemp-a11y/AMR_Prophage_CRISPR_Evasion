# Figures

This directory contains all publication-ready figures in both PDF (vector) and PNG (raster) formats, organized sequentially according to the steps of the computational pipeline.

| Figure File Name | Step | Title / Description | Generating Script |
|---|---|---|---|
| `step0_pipeline_schema` | - | Integrated computational genomics workflow | `scripts/generate_pipeline_schema.py` |
| `step1_assembly_coverage` | Step 1 | Assembly statistics and read coverage verification | `scripts/generate_paper2_plots.py` |
| `step1_prophage_module_map` | Step 1 | Prophage boundary (attL/attR) and functional module map | `scripts/generate_prophage_figures.py` |
| `step2_prophage_composition` | Step 2 | Prophage CDS functional module composition | `scripts/generate_prophage_figures.py` |
| `step2_prophage_synteny` | Step 2 | Comparative synteny map vs. DLP12 reference (clinker) | `scripts/generate_synteny_figure.py` |
| `step2_phage_taxonomy` | Step 2 | Phage taxonomic classification (BLASTp & Foldseek) | `scripts/generate_phage_taxonomy_figure.py` |
| `step3_amr_genetic_context` | Step 3 | AMR genetic context maps and Class 1 integron arrays | `scripts/generate_paper2_plots.py` |
| `step3_virulence_heatmap` | Step 3 | Lineage-wide virulence factor profiling heatmap | `scripts/generate_paper2_plots.py` |
| `step4_crispr_array` | Step 4 | Host Type I-E CRISPR array and spacer alignment profiling | `scripts/generate_prophage_figures.py` |
| `step4_card_phenotype` | Step 4 | Genome-wide CARD antibiotic resistance phenotype prediction | `scripts/generate_exp5_exp6_figures.py` |
| `step4_candidate_prevalence`| Step 4 | Carriage frequency of candidates across ST354 cohorts | `scripts/generate_exp5_exp6_figures.py` |
| `step5_esmfold_enterohemolysin` | Step 5 | ESMFold structure & Foldseek homologs of enterohemolysin | `scripts/generate_esmfold_figure.py` |
| `step5_md_rmsd_rg` | Step 5 | OpenMM explicit-solvent MD trajectory RMSD & Rg profiles | `scripts/generate_exp5_exp6_figures.py` |
| `step5_md_rmsf` | Step 5 | OpenMM explicit-solvent MD trajectory per-residue RMSF | `scripts/generate_exp5_exp6_figures.py` |

## Regenerating All Figures

To regenerate all plots, run the master plotting script:

```bash
conda activate amr_prophage_env
python3 scripts/06_generate_paper2_plots.py --all --output figures/
```
