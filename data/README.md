# Data

This directory contains input sequences, reference files, and pre-computed protein FASTA files.

## Files

| File | Description |
|---|---|
| `prophage_QA5221.gbk` | Prophage region GenBank slice: NODE_1, 56,422–98,747 bp (42,325 bp, 62 CDS) |
| `dlp12_ref.gbk` | DLP12 reference defective pseudoprophage, *E. coli* K-12 MG1655 (NC_001799.1) |
| `enterohemolysin_00061.faa` | Enterohemolysin cargo protein sequence (350 aa, FASTA) |
| `candidates.faa` | All candidate protein sequences (FASTA) |
| `prioritized_candidates.faa` | 5 priority prophage candidate sequences (FASTA) |

## Downloading Raw Reads

Raw Illumina paired-end reads for isolate QA5221:

```bash
fasterq-dump SRR39314025 --outdir data/raw_reads/ --split-files
```

- **SRA Accession:** SRR39314025
- **BioProject:** PRJNA1481519
