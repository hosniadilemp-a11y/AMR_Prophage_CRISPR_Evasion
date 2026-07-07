# Installation Instructions

This document provides steps to set up the execution environment for the AMR Prophage CRISPR Evasion reproducibility pipeline.

## 1. Conda Environment Setup (Recommended)

To install all dependencies (bioinformatics tools, Python library packages, and compilers) in a single step, use the provided `environment.yml`:

```bash
# Navigate to the repository root
cd AMR_Prophage_CRISPR_Evasion

# Create the conda environment
conda env create -f environment/environment.yml

# Activate the environment
conda activate amr_env
```

## 2. Pip Installation (Alternative)

If you already have the external bioinformatics packages (`abricate`, `clinker`, `blast`) installed locally on your system path, you can install the Python dependencies directly via pip:

```bash
pip install -r environment/requirements.txt
```

## 3. GPU Support for OpenMM

If you plan to run long Molecular Dynamics production simulations of the prophage enterohemolysin cytolysin, you will benefit from GPU acceleration. Ensure your CUDA drivers are installed, then verify OpenMM is detecting your GPU:

```bash
python -c "import openmm; print(openmm.pluginFolder)"
python -c "from openmm import Platform; print([Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])"
```

The output should list `CUDA` and `OpenCL` platforms if the GPU setup is working correctly.
