# Paper 2: Active Prophage Mobilome & CRISPR Evasion Analysis

This repository contains the analysis scripts, environment configurations, pipelines, and logs for **Paper 2**: *Genomic Epidemiology, Mobilome Dynamics, and Statistical Resistome Characterization of the Extraintestinal Pathogenic Escherichia coli ST354 Lineage*.

The project centers on clinical isolate **QA5221 (ST354)** and characterizes its active prophage mobilome, CRISPR evasion mechanisms, biocide/AMR co-selection dynamics, and biophysical membrane interactions.

---

## Directory Structure

This repository is structured to mirror the layout of the `AMR_Novel_Gene_Discovery` workspace:

```
repo/
├── .gitignore
├── LICENSE
├── CITATION.cff
├── README.md
├── install_and_run.sh
├── install_and_run.bat
├── config/
│   ├── pipeline_config.yaml    # Config parameters for Paper 2 pipeline
│   └── md_config.yaml          # Parameters for OpenMM simulations
├── environment/
│   ├── environment.yml        # Conda environment specifications
│   ├── requirements.txt       # Python package requirements
│   └── INSTALL.md             # Detailed installation instructions
├── workflows/
│   └── run_all.sh             # Pipeline execution orchestration script
├── scripts/                   # Analysis & figure generation scripts
│   ├── 01_find_att_sites.py   # Sliding-window attL/attR junction scanner
│   ├── 02_extract_prophage_cargo.py # Prophage region gene extraction & annotation
│   ├── 03_run_amr_cooccurrence.py   # Jaccard correlation & network analysis
│   ├── 04_run_experiments_2_5_6.py  # CRISPR profiling & CARD phenotype mapping
│   ├── 05_analyze_paper2_mds.py     # Trajectory biophysical analysis (RMSD, Rg, RMSF)
│   ├── 06_generate_paper2_plots.py  # Plots for Paper 2 figures
│   ├── download_reference_phage.py  # Download DLP12 reference from NCBI
│   └── (plotting_utilities)          # Local LaTeX & plot formatting scripts
├── data/                      # Input genomes & annotation reports (placeholders)
├── docs/                      # Reference notes & plans (placeholders)
├── figures/                   # Output PNG/PDF plots (placeholders)
├── logs/                      # Step execution log files (placeholders)
├── results/                   # Pipeline outputs & trajectories (placeholders)
└── supplementary/             # Supplemental data and files (placeholders)
```

---

## Scientific Focus & Topics Treated

1. **Multidrug Resistance Hotspots**: Profiling the atypical Class 1 integron featuring a `sul3`-linked cassette array on `NODE_37`.
2. **Biocide Co-Selection**: Modeling co-occurrence of quaternary ammonium compound (QAC) biocide resistance genes (`qacL`, `qacEΔ1`) alongside clinical antibiotic cassettes using Jaccard correlation networks.
3. **Active Intact Prophage**: Slicing the 42.3-kb prophage on `NODE_1` and auditing its structural/regulatory core modules to assess lytic viability.
4. **CRISPR-Cas Immune Evasion**: Characterizing the Type I-E Cascade operon/CRISPR array on `NODE_9` and mapping spacer-naive gaps against the active prophage.
5. **Molecular Dynamics (MD) Validation**:
   - Pore-forming enterohemolysin cytolysin (`AC4NUP_00305`) in a POPC/POPG (3:1) bilayer (100 ns simulation).
   - Antirepressor lysogeny-regulation switch (`AC4NUP_00420`) in aqueous solvent (100 ns simulation).

---

## Installation & Running

### Prerequisites
* Linux or Windows OS
* [Conda](https://docs.conda.io/en/latest/) or Miniconda installed

### Setup & Run
To automatically create the Conda environment and run the pipeline:

**On Linux:**
```bash
chmod +x install_and_run.sh
./install_and_run.sh
```

**On Windows:**
```cmd
install_and_run.bat
```

For manual installation and customization options, please refer to [environment/INSTALL.md](environment/INSTALL.md).
