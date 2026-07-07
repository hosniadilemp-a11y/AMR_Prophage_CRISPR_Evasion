# Results

This directory contains the pre-computed outputs from each pipeline step.

## Directory Structure

| Subdirectory | Contents |
|---|---|
| `clinker/` | Interactive HTML clinker synteny visualization (open in browser) |
| `candidates/` | Priority candidate summary TSV + advanced stats JSON |
| `md_structures/` | Prepared PDB structure files for MD simulations |
| `defense_systems/` | DefenseFinder system detection results (TSVs) |

## Key Pre-computed Results

- **Clinker synteny:** Open `clinker/clinker_prophage_synteny.html` in a web browser to view the interactive prophage–DLP12 comparative synteny map.
- **CRISPR immune gap:** 15 spacers from CRISPR-1 (NODE_9) were aligned against the prophage region; no spacer achieves ≥20 bp continuous match.
- **CARD phenotype:** 54 high-confidence hits across 14 drug classes; 8-class Resistant, 5-class Intermediate, carbapenem Susceptible.
- **DefenseFinder:** CRISPR-Cas Type I-E system detected; no anti-CRISPR proteins detected in prophage CDS.
