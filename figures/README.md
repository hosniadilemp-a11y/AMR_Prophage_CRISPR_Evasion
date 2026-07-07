# Figures

This directory contains all publication-ready figures in both PDF (vector) and PNG (raster) formats.

## Figure–Script Mapping

| Figure File | Title | Generating Script |
|---|---|---|
| `figure01_discovery_pipeline_schema` | Integrated pipeline schema | `scripts/generate_pipeline_schema.py` |
| `figure12a_prophage_module_map` | Prophage linear module map (NODE_1) | `scripts/generate_prophage_figures.py` |
| `figure12b_prophage_functional_composition` | Prophage CDS functional composition | `scripts/generate_prophage_figures.py` |
| `figure12c_prophage_synteny` | Comparative synteny vs. DLP12 (clinker) | `scripts/generate_synteny_figure.py` |
| `figure13_card_resistance_phenotype` | CARD resistance phenotype matrix | `scripts/generate_exp5_exp6_figures.py` |
| `figure14_crispr_array` | CRISPR-1 array architecture (NODE_9) | `scripts/generate_prophage_figures.py` |
| `figure15_esmfold_enterohemolysin` | ESMFold pLDDT + Foldseek homologs | `scripts/generate_esmfold_figure.py` |
| `figure16_phage_taxonomy` | Phage taxonomy (BLASTp + Foldseek) | `scripts/generate_phage_taxonomy_figure.py` |
| `figure17_candidate_prevalence` | Lineage-wide candidate prevalence | `scripts/generate_exp5_exp6_figures.py` |
| `figure10_md_rmsd_rg` | MD RMSD + Rg panel | `scripts/generate_exp5_exp6_figures.py` |
| `figure10_md_rmsf` | MD per-residue RMSF panel | `scripts/generate_exp5_exp6_figures.py` |
| `figure_amr_network` | AMR Jaccard co-occurrence network | `scripts/generate_paper2_plots.py` |
| `figure_virulence_heatmap` | Lineage virulence factor heatmap | `scripts/generate_paper2_plots.py` |

## Regenerating All Figures

```bash
conda activate amr_prophage_env
python3 scripts/06_generate_paper2_plots.py --all --output figures/
```
