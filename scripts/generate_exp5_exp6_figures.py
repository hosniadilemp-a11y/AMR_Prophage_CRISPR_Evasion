#!/usr/bin/env python3
"""
Generate publication figures for Experiments 5 and 6:
 - Figure 13: CARD resistance phenotype heatmap (acquired + intrinsic)
 - Figure 14: CRISPR array diagram
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from io import StringIO
from pathlib import Path
import subprocess

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})

OUT_PLOTS  = Path("/media/adel/Data/Hosni/openmm_windows_setup/AMR_Work/manuscript/paper2/plots")
OUT_REPORT = Path("/media/adel/Data/Hosni/openmm_windows_setup/AMR_Work/reports")

# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENT 6 FIGURE: Resistance phenotype matrix
# ─────────────────────────────────────────────────────────────────────────────

# Load CARD data
raw = open(OUT_REPORT / "Step_CARD_RGI/QA5221_card_results.tsv").readlines()
hi  = next(i for i,l in enumerate(raw) if l.startswith('#FILE'))
keep = [raw[hi]] + [l for l in raw[hi+1:] if l.strip() and
        not any(l.startswith(x) for x in ('Using ','Processing','Found ','Tip:','Done'))]
df = pd.read_csv(StringIO(''.join(keep)), sep='\t')
df.columns = [c.lstrip('#') for c in df.columns]
df_acq = df[(df['%COVERAGE'] >= 80) & (df['%IDENTITY'] >= 85)].copy()

# Acquired resistance genes (from ResFinder/paper2 Table)
ACQUIRED = {'TEM-1':'blaTEM-1B', "APH(3')-Ia":'aph(3′)-Ia', 'aadA':'aadA1',
            'aadA2':'aadA2', 'sul1':'sul1', 'sul3':'sul3', 'dfrA7':'dfrA7',
            'cmlA1':'cmlA1', 'estX':'estX', 'tet(A)':'tet(A)', 'qacL':'qacL',
            'qacEdelta1':'qacEΔ1'}

# Build resistance matrix: rows=antibiotic class, cols=[Acquired Genes, Intrinsic Efflux]
CLASSES = [
    ('Beta-lactam',        ['TEM-1'],              ['acrA','acrB','TolC'],  'R',  '#d62828'),
    ('Aminoglycoside',     ["APH(3')-Ia",'aadA','aadA2'], ['acrD'],          'R',  '#d62828'),
    ('Sulfonamide',        ['sul1','sul3'],         [],                      'R',  '#d62828'),
    ('Trimethoprim',       ['dfrA7'],               [],                      'R',  '#d62828'),
    ('Phenicol',           ['cmlA1'],               ['acrA','acrB'],         'R',  '#d62828'),
    ('Macrolide',          ['estX'],                ['acrA','acrB','mdfA'],  'R',  '#d62828'),
    ('Tetracycline',       ['tet(A)'],              ['acrA','acrB','mdfA'],  'R',  '#d62828'),
    ('Biocide/QAC',        ['qacL','qacEdelta1'],   ['acrA','acrB'],         'R',  '#d62828'),
    ('Fluoroquinolone',    [],                      ['acrA','acrB','acrEF','emrAB'], 'I', '#f4a261'),
    ('Carbapenem',         [],                      ['acrA','acrB','TolC'],  'S',  '#2a9d8f'),
    ('Polymyxin/Peptide',  [],                      ['pmrF','eptA'],         'I',  '#f4a261'),
    ('Aminocoumarin',      [],                      ['TolC','mdtABC'],       'I',  '#f4a261'),
    ('Rifamycin',          [],                      ['acrA','acrB','TolC'],  'I',  '#f4a261'),
    ('Glycylcycline',      [],                      ['acrA','acrB'],         'I',  '#f4a261'),
]

fig, ax = plt.subplots(figsize=(13, 7.5))

n_classes = len(CLASSES)
col_w     = 0.28  # column width per gene category

# Headers
ax.text(0.5, n_classes + 0.6, 'Antibiotic Class', ha='center', va='bottom',
        fontsize=10, fontweight='bold', color='#1a1a2e')
ax.text(3.5, n_classes + 0.6, 'Acquired Resistance Genes', ha='center', va='bottom',
        fontsize=10, fontweight='bold', color='#d62828')
ax.text(7.2, n_classes + 0.6, 'Intrinsic Efflux/Regulatory', ha='center', va='bottom',
        fontsize=10, fontweight='bold', color='#457b9d')
ax.text(9.6, n_classes + 0.6, 'Predicted\nPhenotype', ha='center', va='bottom',
        fontsize=10, fontweight='bold', color='#1a1a2e')

for i, (cls, acq, intr, phen, color) in enumerate(CLASSES):
    y = n_classes - i - 1
    # Row background
    bg = '#fff8f0' if i % 2 == 0 else '#f8f9fa'
    ax.barh(y, 10.5, height=0.85, left=-0.2, color=bg, zorder=0)

    # Antibiotic class label
    ax.text(0.4, y, cls, va='center', ha='left', fontsize=9,
            fontweight='bold' if phen == 'R' else 'normal', color='#1a1a2e')

    # Acquired gene badges
    for j, g in enumerate(acq):
        bx = 2.8 + j * 1.1
        rect = mpatches.FancyBboxPatch((bx - 0.42, y - 0.29), 0.84, 0.58,
                                        boxstyle='round,pad=0.05',
                                        facecolor='#d62828', edgecolor='none', zorder=2)
        ax.add_patch(rect)
        ax.text(bx, y, g.replace("APH(3')-Ia","aph3'Ia"),
                ha='center', va='center', fontsize=7.5, color='white', fontweight='bold')

    # Intrinsic gene badges
    for j, g in enumerate(intr[:4]):
        bx = 5.9 + j * 0.93
        rect = mpatches.FancyBboxPatch((bx - 0.40, y - 0.25), 0.80, 0.50,
                                        boxstyle='round,pad=0.05',
                                        facecolor='#457b9d', edgecolor='none', zorder=2, alpha=0.85)
        ax.add_patch(rect)
        ax.text(bx, y, g, ha='center', va='center', fontsize=7, color='white')

    # Phenotype badge
    PHEN_COL = {'R': '#d62828', 'I': '#f4a261', 'S': '#2a9d8f'}
    PHEN_TXT = {'R': 'RESISTANT', 'I': 'INTERMEDIATE', 'S': 'SUSCEPTIBLE'}
    pc = PHEN_COL[phen]
    rect2 = mpatches.FancyBboxPatch((9.0, y - 0.24), 1.2, 0.48,
                                     boxstyle='round,pad=0.05',
                                     facecolor=pc, edgecolor='none', zorder=2)
    ax.add_patch(rect2)
    ax.text(9.6, y, PHEN_TXT[phen], ha='center', va='center',
            fontsize=7, color='white', fontweight='bold')

# Column dividers
ax.axvline(2.6, color='#dee2e6', lw=1, zorder=1, ymin=0, ymax=1)
ax.axvline(9.0, color='#dee2e6', lw=1, zorder=1, ymin=0, ymax=1)

# Legend
leg_handles = [
    mpatches.Rectangle((0,0),1,1, facecolor='#d62828', label='Acquired resistance gene (mobile element)'),
    mpatches.Rectangle((0,0),1,1, facecolor='#457b9d', label='Intrinsic efflux pump / regulator'),
    mpatches.Rectangle((0,0),1,1, facecolor='#f4a261', label='Intermediate — phenotype depends on expression'),
]
ax.legend(handles=leg_handles, loc='lower right', fontsize=8, framealpha=0.92)

ax.set_xlim(-0.2, 10.8)
ax.set_ylim(-0.6, n_classes + 1.0)
ax.axis('off')
ax.set_title(
    'Predicted Antibiotic Resistance Phenotype of E. coli QA5221 (CARD Database, n=54 genes)',
    fontsize=11, fontweight='bold', pad=10)

plt.tight_layout()
for ext in ('pdf','png'):
    fig.savefig(OUT_PLOTS / f'figure13_card_resistance_phenotype.{ext}', dpi=300, bbox_inches='tight')
print("✓ Figure 13 (CARD phenotype matrix) saved.")
plt.close()

# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENT 5 FIGURE: CRISPR array diagram
# ─────────────────────────────────────────────────────────────────────────────
SPACERS = [
    "GTGACCATAAGGGCATTTCATACTGCCTTGCT",
    "GTGAACGTGACAACGCCCGCGACAGAAACGTC",
    "GGAAGAAGCGGTAAAGCGGGTTCCAGTTCAGC",
    "AAACTCCGTGGTGTAGTGCGTTTGTTGGTGCG",
    "TTGGTATTAACGGGAATGTAGTCACTCCAGGA",
    "GCGATCGTATTTTCCGCTACGTCGCTCAAATG",
    "TCGTTGATGGTGTGATGAGTGTACAAAAAAAA",
    "ACGACTGAATCCGGCGAGGCTATCAACTGGCC",
    "CCGCACACGTCGAGCTGGTGGGGATTAATGCT",
    "CAGTAAACGTTCGCGCGTTCTCGCGCTCACGA",
    "TCTGTGTATTCGCTACCAGCCAGCACCGGTAC",
    "CAGGGTAACTGGCTGCTCCTCAACCATTCCGC",
    "GGGTTTTAGCGACCGGACCAGAGTACGCAACA",
    "TTGCGCCTTAAGTTCCTTAACTTCCTCTGTAT",
    "GGAGTGATGAATTGGAAGTAACGGATAACCGC",
]
REPEAT = "GTGTTCCCTGCGCCAGCGGGGATAAACCG"
N = len(SPACERS)

fig2, ax2 = plt.subplots(figsize=(14, 3.5))
ax2.set_xlim(-1, N * 2 + 0.5)
ax2.set_ylim(-0.8, 2.0)
ax2.axis('off')

R_H, S_H = 0.4, 0.4
Y = 0.8

for i in range(N):
    xr = i * 2
    xs = xr + 1
    # Repeat (diamond-ish rectangle)
    ax2.add_patch(mpatches.FancyBboxPatch((xr, Y - R_H/2), 0.9, R_H,
                   boxstyle='round,pad=0.05', facecolor='#264653', edgecolor='white', lw=0.7))
    ax2.text(xr + 0.45, Y, 'DR', ha='center', va='center', fontsize=6.5,
             color='white', fontweight='bold')
    # Spacer
    ax2.add_patch(mpatches.FancyBboxPatch((xs, Y - S_H/2), 0.9, S_H,
                   boxstyle='round,pad=0.05', facecolor='#e9c46a', edgecolor='white', lw=0.7))
    ax2.text(xs + 0.45, Y, f'Sp{i+1}', ha='center', va='center', fontsize=6.5,
             color='#1a1a2e', fontweight='bold')
    # Spacer sequence below
    ax2.text(xs + 0.45, Y - S_H/2 - 0.12, SPACERS[i][:16]+'…',
             ha='center', va='top', fontsize=5.5, color='#495057', family='monospace')

# Final repeat
xr_last = N * 2
ax2.add_patch(mpatches.FancyBboxPatch((xr_last, Y - R_H/2), 0.9, R_H,
               boxstyle='round,pad=0.05', facecolor='#264653', edgecolor='white', lw=0.7))
ax2.text(xr_last + 0.45, Y, 'DR', ha='center', va='center', fontsize=6.5,
         color='white', fontweight='bold')

# Labels
ax2.text(-0.5, Y, 'CRISPR-1\n(NODE_9)', ha='right', va='center',
         fontsize=9, fontweight='bold', color='#1a1a2e')
ax2.text(N * 2 + 1.0, Y, f'36,345–37,288 bp\n16 repeats × 29 bp\n15 spacers × 32 bp',
         ha='left', va='center', fontsize=8, color='#495057')

# No-match annotation
ax2.annotate('', xy=(N * 2 + 0.5, 1.5), xytext=(0.45, 1.5),
             arrowprops=dict(arrowstyle='-', color='#adb5bd', lw=1.5))
ax2.text(N + 0.45, 1.65, '✗ No spacer targets the prophage on NODE_1',
         ha='center', va='bottom', fontsize=9, color='#d62828', fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#fff3cd', edgecolor='#d62828', lw=1))

ax2.set_title(
    'CRISPR Array Detected in E. coli QA5221 (NODE_9, 36,345–37,288 bp)\n'
    'None of the 15 spacers show sequence identity to the prophage region (NODE_1)',
    fontsize=10, fontweight='bold')

plt.tight_layout()
for ext in ('pdf','png'):
    fig2.savefig(OUT_PLOTS / f'figure14_crispr_array.{ext}', dpi=300, bbox_inches='tight')
print("✓ Figure 14 (CRISPR array diagram) saved.")
plt.close()

print("\nAll figures saved successfully.")
