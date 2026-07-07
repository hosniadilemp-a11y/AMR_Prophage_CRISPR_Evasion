#!/usr/bin/env python3
"""
Figure 16: Phage Taxonomic Classification — BLASTp + Foldseek evidence
Two-panel publication figure showing:
  (A) BLASTp identity distribution for terminase (00091) and capsid (00109)
  (B) Foldseek structural taxonomy across databases
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9,
    'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False,
})

OUT = Path("/media/adel/Data/Hosni/openmm_windows_setup/AMR_Work/manuscript/paper2/plots")

# ── Data ──────────────────────────────────────────────────────────────────────
# BLASTp top hits (from NCBI remote search)
blast_term = {
    "Escherichia coli WP_001118115.1":  99.6,
    "Escherichia coli WP_069723171.1":  99.6,
    "Escherichia coli WP_001118123.1":  99.2,
    "Escherichia coli WP_219788440.1":  99.2,
    "Escherichia coli WP_335102787.1":  99.2,
    "QA5221 (self)":                   100.0,
}
blast_caps = {
    "Escherichia coli WP_000042293.1": 100.0,
    "Escherichia coli WP_318074569.1":  99.7,
    "Escherichia coli O156 EEZ64907":   99.7,
    "Escherichia coli WP_203778530.1":  99.7,
    "Escherichia coli WP_057781017.1":  99.7,
    "QA5221 (self)":                   100.0,
}

# Foldseek structural homologs (afdb50) — phage-annotated hits
foldseek_term = [
    ("Klebsiella aerogenes",           89.8, "Terminase small subunit"),
    ("Salmonella enterica",            87.5, "Putative prophage terminase small subunit"),
    ("Citrobacter portucalensis",      71.4, "Terminase"),
    ("Klebsiella pneumoniae",          65.9, "Putative bacteriophage protein"),
    ("Kosakonia oryziphila",           57.7, "Phage terminase small subunit"),
    ("Klebsiella michiganensis",       57.3, "Terminase (Fragment)"),
    ("Serratia marcescens",            54.2, "Terminase"),
]
foldseek_caps = [
    ("Serratia sp. Ag1",               70.4, "Bacteriophage protein"),
    ("Salmonella enterica Houtenae",   44.7, "Bacteriophage protein"),
    ("Cronobacter turicensis",         41.5, "Bacteriophage protein"),
    ("Enterobacteriaceae bacterium",   38.4, "Bacteriophage protein"),
    ("Sodalis glossinidius",           39.7, "Bacteriophage protein"),
    ("Jejubacter calystegiae",         48.4, "Bacteriophage protein"),
    ("Morganella sp.",                 36.4, "Bacteriophage protein"),
]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(
    'Phage Taxonomic Classification of Prophage Structural Genes in QA5221\n'
    'KNGPFPPJ_00091 (Terminase) and KNGPFPPJ_00109 (Major Capsid Protein)',
    fontsize=12, fontweight='bold', y=0.98
)

COLORS = {'self': '#2a9d8f', 'ecoli': '#457b9d', 'other': '#adb5bd'}
PHAGE_COL = '#e76f51'
SIPHO_COL = '#457b9d'

# ── Panel A: BLASTp terminase ─────────────────────────────────────────────────
ax = axes[0, 0]
labels_t = list(blast_term.keys())
vals_t   = list(blast_term.values())
y_pos = np.arange(len(labels_t))[::-1]
cols_t = ['#2a9d8f' if 'self' in l else '#457b9d' for l in labels_t]
bars = ax.barh(y_pos, vals_t, color=cols_t, edgecolor='white', height=0.6)
for bar, val in zip(bars, vals_t):
    ax.text(bar.get_width() - 1.5, bar.get_y() + bar.get_height()/2,
            f'{val:.1f}%', va='center', ha='right', fontsize=8,
            color='white', fontweight='bold')
ax.set_yticks(y_pos)
ax.set_yticklabels([l[:35] for l in labels_t], fontsize=8)
ax.set_xlim(96, 101.5)
ax.set_xlabel('Sequence Identity (%)', fontsize=9)
ax.set_title('A. BLASTp — Terminase (KNGPFPPJ_00091, 249 aa)\nNCBI nr database | All hits: Escherichia coli',
             fontsize=9.5, fontweight='bold')
ax.text(0.02, 0.05, '100% E. coli-specific\n→ Intra-species prophage\nconservation',
        transform=ax.transAxes, fontsize=8,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#e8f4f8', edgecolor='#457b9d'))

# ── Panel B: BLASTp capsid ────────────────────────────────────────────────────
ax = axes[0, 1]
labels_c = list(blast_caps.keys())
vals_c   = list(blast_caps.values())
y_pos = np.arange(len(labels_c))[::-1]
cols_c = ['#2a9d8f' if 'self' in l else '#457b9d' for l in labels_c]
bars = ax.barh(y_pos, vals_c, color=cols_c, edgecolor='white', height=0.6)
for bar, val in zip(bars, vals_c):
    ax.text(bar.get_width() - 1.5, bar.get_y() + bar.get_height()/2,
            f'{val:.1f}%', va='center', ha='right', fontsize=8,
            color='white', fontweight='bold')
ax.set_yticks(y_pos)
ax.set_yticklabels([l[:35] for l in labels_c], fontsize=8)
ax.set_xlim(97.5, 101.5)
ax.set_xlabel('Sequence Identity (%)', fontsize=9)
ax.set_title('B. BLASTp — Major Capsid Protein (KNGPFPPJ_00109, 344 aa)\nNCBI nr database | All hits: Escherichia coli',
             fontsize=9.5, fontweight='bold')
ax.text(0.02, 0.05, '100% E. coli-specific\n→ Confirms E. coli\nprophage origin',
        transform=ax.transAxes, fontsize=8,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#e8f4f8', edgecolor='#457b9d'))

# ── Panel C: Foldseek terminase ───────────────────────────────────────────────
ax = axes[1, 0]
labels_ft = [f"{h[0][:28]}\n{h[2][:30]}" for h in foldseek_term]
vals_ft   = [h[1] for h in foldseek_term]
y_pos = np.arange(len(labels_ft))[::-1]
cols_ft = ['#e76f51' if v > 80 else '#f4a261' if v > 60 else '#adb5bd' for v in vals_ft]
bars = ax.barh(y_pos, vals_ft, color=cols_ft, edgecolor='white', height=0.55)
for bar, val in zip(bars, vals_ft):
    ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
            f'{val:.1f}%', va='center', fontsize=8, fontweight='bold')
ax.axvline(50, color='#adb5bd', linewidth=1, linestyle='--', alpha=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(labels_ft, fontsize=7.5)
ax.set_xlim(0, 105)
ax.set_xlabel('Sequence Identity (%)', fontsize=9)
ax.set_title('C. Foldseek — Terminase Structural Homologs (afdb50)\nAll named as "Terminase" / "Prophage terminase"',
             fontsize=9.5, fontweight='bold')
# Taxonomy box
ax.text(0.98, 0.05,
        'Host range:\nEnterobacteriaceae\n(Klebsiella, Salmonella,\nCitrobacter, Serratia)\n→ Siphoviridae / Lambdoid',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#fff3e0', edgecolor='#e76f51'))

# ── Panel D: Foldseek capsid ──────────────────────────────────────────────────
ax = axes[1, 1]
labels_fc = [f"{h[0][:28]}\n{h[2][:30]}" for h in foldseek_caps]
vals_fc   = [h[1] for h in foldseek_caps]
y_pos = np.arange(len(labels_fc))[::-1]
cols_fc = ['#e76f51' if v > 60 else '#f4a261' if v > 40 else '#adb5bd' for v in vals_fc]
bars = ax.barh(y_pos, vals_fc, color=cols_fc, edgecolor='white', height=0.55)
for bar, val in zip(bars, vals_fc):
    ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
            f'{val:.1f}%', va='center', fontsize=8, fontweight='bold')
ax.axvline(40, color='#adb5bd', linewidth=1, linestyle='--', alpha=0.7,
           label='40% identity threshold')
ax.set_yticks(y_pos)
ax.set_yticklabels(labels_fc, fontsize=7.5)
ax.set_xlim(0, 90)
ax.set_xlabel('Sequence Identity (%)', fontsize=9)
ax.set_title('D. Foldseek — Capsid Structural Homologs (afdb50)\nAll annotated as "Bacteriophage protein"',
             fontsize=9.5, fontweight='bold')
ax.legend(fontsize=8, loc='lower right')
ax.text(0.98, 0.05,
        'Host range:\nEnterobacteriaceae\n(Serratia, Salmonella,\nCronobacter, Morganella)\n→ Enterobacterial Siphovirus',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#fff3e0', edgecolor='#e76f51'))

plt.tight_layout(rect=[0, 0, 1, 0.96])

# ── Add overall taxonomy verdict ──────────────────────────────────────────────
fig.text(0.5, 0.01,
         '▶  Taxonomic Conclusion: QA5221 prophage belongs to a lambdoid Siphoviridae family '
         'infecting Enterobacteriaceae (E. coli, Salmonella, Klebsiella, Serratia, Citrobacter)',
         ha='center', fontsize=10, fontweight='bold', color='#e76f51',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='#fff3e0', edgecolor='#e76f51', alpha=0.9))

for ext in ('pdf', 'png'):
    fig.savefig(OUT / f'figure16_phage_taxonomy.{ext}', dpi=300, bbox_inches='tight')
print("✓ Figure 16 (phage taxonomy) saved.")
plt.close()
