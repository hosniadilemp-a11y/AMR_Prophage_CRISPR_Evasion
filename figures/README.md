# Figures

This directory contains all publication-ready figures in both PDF (vector) and PNG (raster) formats.

All figures correspond to either the main manuscript or the supplementary materials.

## Main Manuscript Figures

| Figure File | Title | Generating Script |
|---|---|---|
| `figure01_discovery_pipeline_schema` | Integrated computational pipeline schema | `scripts/generate_pipeline_schema.py` |
| `figure12a_prophage_module_map` | Prophage linear functional module map (NODE_1, 42.3 kb) | `scripts/generate_prophage_figures.py` |
| `figure12b_prophage_functional_composition` | Prophage CDS functional composition (donut + bar) | `scripts/generate_prophage_figures.py` |
| `figure12c_prophage_synteny` | Comparative synteny vs. DLP12 pseudoprophage (clinker) | `scripts/generate_synteny_figure.py` |
| `figure13_card_resistance_phenotype` | CARD resistance phenotype matrix (54 hits, 14 drug classes) | `scripts/generate_exp5_exp6_figures.py` |
| `figure14_crispr_array` | CRISPR-1 array architecture and spacer alignment (NODE_9) | `scripts/generate_prophage_figures.py` |
| `figure15_esmfold_enterohemolysin` | ESMFold pLDDT profile + Foldseek structural homologs | `scripts/generate_esmfold_figure.py` |
| `figure16_phage_taxonomy` | Phage taxonomic classification (BLASTp + Foldseek) | `scripts/generate_phage_taxonomy_figure.py` |
| `figure17_candidate_prevalence` | Lineage-wide candidate prevalence (32 / 400 / 800 genomes) | `scripts/generate_exp5_exp6_figures.py` |

## Supplementary Figures

| Figure File | Title | Generating Script |
|---|---|---|
| `figure1_assembly_coverage` | Assembly statistics and read coverage (QA5221 baseline) | `scripts/generate_paper2_plots.py` |
| `figure2_pangenome_partition` | ST354 pangenome partition summary | `scripts/generate_paper2_plots.py` |
| `figure3_phylogeny_tree` | Core-genome ML phylogenetic tree (ST354 cohort) | `scripts/generate_paper2_plots.py` |
| `figure4_is_distribution` | Insertion sequence element distribution (22 IS, 9 families) | `scripts/generate_paper2_plots.py` |
| `figure5_stats_comparison` | Statistical comparison across genomic partitions | `scripts/generate_paper2_plots.py` |
| `figure6_amr_comparison` | AMR gene comparison across ST354 cohort | `scripts/generate_paper2_plots.py` |
| `figure7_amr_genetic_context` | AMR genetic context maps (integron arrays) | `scripts/generate_paper2_plots.py` |
| `figure8_amr_gc3_scatter` | GC vs. GC3 codon bias scatter (resistome amelioration) | `scripts/generate_paper2_plots.py` |
| `figure10_md_rmsd_rg` | MD RMSD + Rg panel (enterohemolysin + antirepressor) | `scripts/generate_exp5_exp6_figures.py` |
| `figure10_md_rmsf` | MD per-residue RMSF panel | `scripts/generate_exp5_exp6_figures.py` |
| `figure_amr_network` | AMR Jaccard co-occurrence network | `scripts/run_amr_cooccurrence.py` |
| `figure_virulence_heatmap` | Lineage-wide virulence factor heatmap (ST354 cohort) | `scripts/generate_paper2_plots.py` |
| `figureX_pangenome_curves` | Pangenome accumulation curves (Heap's Law fit) | `scripts/generate_paper2_plots.py` |

## Regenerating All Figures

```bash
conda activate amr_prophage_env
python3 scripts/06_generate_paper2_plots.py --all --output figures/
```
