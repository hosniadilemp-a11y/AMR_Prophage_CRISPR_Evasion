#!/usr/bin/env python3
"""
Wrapper to generate Figure 4 (and Figure 19 combined) using the authentic 
publication pipeline in EXP_05_Figures_Legibility_AltText/generate_figure4_labeled.py.
"""
import sys
import subprocess
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
EXP_05_SCRIPT = SCRIPT_DIR.parents[3] / "01_Decision_and_Checklists" / "Experiments_Tracking" / "EXP_05_Figures_Legibility_AltText" / "generate_figure4_labeled.py"

if EXP_05_SCRIPT.exists():
    ret = subprocess.run([sys.executable, str(EXP_05_SCRIPT)], check=True)
    sys.exit(ret.returncode)
else:
    print(f'ERROR: {EXP_05_SCRIPT} not found!')
    sys.exit(1)
