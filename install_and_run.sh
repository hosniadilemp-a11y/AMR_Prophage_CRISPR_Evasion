#!/bin/bash
# =======================================================================
#          Paper 2 Reproducibility Pipeline Automated Linux Launcher
# =======================================================================
set -e

echo "Checking for Conda installation..."

# Source conda environment scripts if conda is not in PATH
if ! command -v conda &> /dev/null; then
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    else
        echo "[ERROR] Conda was not found in your PATH or common installation directories."
        echo "Please make sure conda is installed and initialized."
        exit 1
    fi
fi

# Ensure conda function is available
if ! declare -F conda > /dev/null; then
    eval "$(conda shell.bash hook)"
fi

# Check if amr_env already exists
if conda env list | grep -q "amr_env"; then
    echo "Environment 'amr_env' already exists."
else
    echo "Creating environment 'amr_env' from environment/environment.yml (this may take a few minutes)..."
    conda env create -f environment/environment.yml -y
fi

echo "Activating environment 'amr_env'..."
conda activate amr_env

echo "Running Paper 2 Reproducibility Pipeline..."
bash workflows/run_all.sh

echo ""
echo "Pipeline execution completed successfully."
