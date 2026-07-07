# Logs

This directory contains execution logs from the pipeline steps.

## Log File Registry

| Log File | Step | Description |
|---|---|---|
| `step1_att_sites.log` | Step 1 | Sliding-window att-site scanner output |
| `step3_abricate_card.log` | Step 3 | ABRicate CARD screening log |
| `step4_minced.log` | Step 4 | MinCED CRISPR detection log |
| `step4_blastn_spacers.log` | Step 4 | BLASTn spacer alignment log |
| `step5_md_ehly.csv` | Step 5 | Enterohemolysin bilayer MD thermodynamic log |
| `step5_md_antirepressor.csv` | Step 5 | Antirepressor aqueous MD thermodynamic log |

> **Note:** Large MD trajectory files (.dcd) are not stored here due to GitHub size limits.
> They are archived at Zenodo: https://doi.org/10.5281/zenodo.21073430
