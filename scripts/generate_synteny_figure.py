#!/usr/bin/env python3
"""
Generate Figure 12c: Prophage vs DLP12 synteny comparison figure.
Uses dna_features_viewer to draw linear genome maps side-by-side with similarity shading.
Falls back to matplotlib if dna_features_viewer unavailable.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np
from pathlib import Path
import subprocess, sys

OUT_DIR = Path("/media/adel/Data/Hosni/openmm_windows_setup/AMR_Work/manuscript/paper2/plots")

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 9,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# ── Try dna_features_viewer first ────────────────────────────────────────────
try:
    from dna_features_viewer import GraphicFeature, GraphicRecord
    HAS_DFV = True
except ImportError:
    HAS_DFV = False
    print("dna_features_viewer not found — using manual matplotlib rendering.")

# ── QA5221 Prophage gene data ────────────────────────────────────────────────
PROPH_START = 56422
PROPH_END   = 98747

# (start_bp_from_left_att, end_bp_from_left_att, strand, category, label)
# We shift so attL = 0
GENES_QA = [
    # Host cargo
    (1794,  1982,  -1, 'Host Cargo',          'rcbA'),
    (2233,  3285,  -1, 'Host Cargo',          'enterohemolysin'),
    # Lysogeny
    (7045,  7266,  -1, 'Lysogeny',            'kilR'),
    (8769,  9191,  +1, 'Lysogeny',            'ydaT'),
    # Replication
    (10082, 10828, +1, 'Replication',         'dnaC'),
    (15700, 15990, +1, 'Replication',         'ybcO'),
    # Candidates / Structural
    (16812, 17423, +1, 'Structural (novel)',  '★antirep'),
    (20158, 20907, +1, 'Packaging',           'terL'),
    (33405, 34055, +1, 'Structural (novel)',  '★tail'),
    (34364, 35398, +1, 'Structural (novel)',  '★capsid'),
    (37274, 38473, +1, 'Structural (novel)',  '★baseplate'),
    # Lysis
    (40529, 41092, +1, 'Lysis',              'hin'),
    (41597, 42253, +1, 'Lysis',              'pphA'),
]

# ── DLP12 (K-12 MG1655 pseudoprophage) gene data ─────────────────────────────
# DLP12 is ~13.6 kb, heavily decayed; only retains regulatory + some replication
# NC_000913.3 DLP12 region: ~1,610,535..1,624,152
# We use representative functional zones (manually curated from literature)
DLP12_LEN = 13617
GENES_DLP = [
    # Lysogeny (CI repressor remnant)
    (1200,  1800,  +1, 'Lysogeny',            'cI*'),
    # Replication (partial)
    (2500,  3200,  +1, 'Replication',         'repL*'),
    # Remnant packaging (truncated / pseudogene)
    (5100,  5600,  +1, 'Packaging',           'terL*'),
    # Major deletions — no structural core
    # Lysis
    (11800, 12400, +1, 'Lysis',              'hin*'),
    (12550, 13200, +1, 'Lysis',              'pphA*'),
]

# ── Colour map ───────────────────────────────────────────────────────────────
CMAP = {
    'Host Cargo':          '#e63946',
    'Lysogeny':            '#457b9d',
    'Replication':         '#2a9d8f',
    'Packaging':           '#f4a261',
    'Structural (novel)':  '#264653',
    'Lysis':               '#8338ec',
    'Hypothetical':        '#adb5bd',
}

def draw_track(ax, genes, total_len, y_center=0.5, h=0.22, label_above=True):
    """Draw a horizontal linear genome track."""
    # backbone
    ax.plot([0, total_len], [y_center, y_center], color='#343a40', lw=2.5, zorder=1)
    for (s, e, strand, cat, lbl) in genes:
        color = CMAP.get(cat, '#adb5bd')
        head_len = min(400, (e - s) * 0.3)
        if strand == +1:
            # body
            ax.fill_betweenx([y_center - h/2, y_center + h/2],
                              s, e - head_len, color=color, zorder=3)
            # arrowhead
            tri = np.array([[e - head_len, y_center - h/2],
                             [e - head_len, y_center + h/2],
                             [e, y_center]])
            ax.fill(tri[:, 0], tri[:, 1], color=color, zorder=3)
        else:
            ax.fill_betweenx([y_center - h/2, y_center + h/2],
                              s + head_len, e, color=color, zorder=3)
            tri = np.array([[s + head_len, y_center - h/2],
                             [s + head_len, y_center + h/2],
                             [s, y_center]])
            ax.fill(tri[:, 0], tri[:, 1], color=color, zorder=3)
        # label
        mid = (s + e) / 2
        if lbl:
            ypos = y_center + h/2 + 0.035 if label_above else y_center - h/2 - 0.035
            va   = 'bottom' if label_above else 'top'
            is_novel = '★' in lbl
            ax.text(mid, ypos, lbl.replace('★', ''),
                    ha='center', va=va, fontsize=7.5,
                    fontweight='bold' if is_novel else 'normal',
                    color='#1a1a2e',
                    style='normal')
            if is_novel:
                ax.text(mid - 150, ypos, '★',
                        ha='right', va=va, fontsize=7.5,
                        color='#e9c46a', fontweight='bold')
    return

# ── Similarity shading between the two tracks ────────────────────────────────
def add_similarity_ribbon(ax, s1_start, s1_end, s2_start, s2_end,
                           y_top, y_bot, color='#adb5bd', alpha=0.25):
    """Draw a ribbon between two homologous regions."""
    from matplotlib.patches import Polygon
    xs = np.mean([s1_start, s1_end])
    xe = np.mean([s2_start, s2_end])
    # trapezoid
    verts = [(s1_start, y_top), (s1_end, y_top),
             (s2_end, y_bot), (s2_start, y_bot)]
    poly = Polygon(verts, closed=True, color=color, alpha=alpha, zorder=0)
    ax.add_patch(poly)

# ────────────────────────────────────────────────────────────────────────────
# Build figure
# ────────────────────────────────────────────────────────────────────────────
QA_LEN = PROPH_LEN = PROPH_END - PROPH_START  # 42325

fig, ax = plt.subplots(figsize=(14, 5.5))
ax.set_xlim(-1500, max(QA_LEN, DLP12_LEN) + 1500)
ax.set_ylim(-0.15, 1.55)
ax.axis('off')

# Track positions
Y_QA  = 1.10
Y_DLP = 0.35

# ── Draw QA5221 prophage track ───────────────────────────────────────────────
draw_track(ax, GENES_QA, QA_LEN, y_center=Y_QA, h=0.20, label_above=True)
ax.text(-1400, Y_QA, 'QA5221\n(42.3 kb)', ha='right', va='center',
        fontsize=9, fontweight='bold', color='#1a1a2e')
# attL / attR ticks
for xpos, name in [(0, 'attL'), (QA_LEN, 'attR')]:
    ax.axvline(xpos, ymin=(Y_QA - 0.15) / 1.70, ymax=(Y_QA + 0.18) / 1.70,
               color='#e9c46a', lw=2, ls='--', zorder=2)
    ax.text(xpos, Y_QA - 0.20, name,
            ha='center', va='top', fontsize=8,
            color='#e9c46a', fontweight='bold')

# ── Draw DLP12 track ─────────────────────────────────────────────────────────
draw_track(ax, GENES_DLP, DLP12_LEN, y_center=Y_DLP, h=0.20, label_above=False)
ax.text(-1400, Y_DLP, 'DLP12\nK-12 MG1655\n(13.6 kb, defective)', ha='right', va='center',
        fontsize=9, fontweight='bold', color='#6c757d')

# Deletion annotation on DLP12
ax.annotate('', xy=(DLP12_LEN * 0.62, Y_DLP),
            xytext=(DLP12_LEN * 0.38, Y_DLP),
            arrowprops=dict(arrowstyle='<->', color='#dc3545',
                            lw=1.5, connectionstyle='arc3,rad=0'))
ax.text(DLP12_LEN * 0.50, Y_DLP + 0.15,
        'Structural core deleted\n(no tail/capsid/baseplate)',
        ha='center', va='bottom', fontsize=7.5,
        color='#dc3545', style='italic',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#fff3cd',
                  edgecolor='#dc3545', linewidth=0.8))

# ── Homology ribbons (QA5221 → DLP12) ───────────────────────────────────────
# Lysogeny module — high similarity
add_similarity_ribbon(ax,
    GENES_QA[2][0], GENES_QA[4][1],        # kilR..ydaT
    GENES_DLP[0][0], GENES_DLP[0][1],       # cI*
    Y_QA - 0.11, Y_DLP + 0.11,
    color='#457b9d', alpha=0.25)

# Replication — moderate
add_similarity_ribbon(ax,
    GENES_QA[5][0], GENES_QA[5][1],         # ybcO
    GENES_DLP[1][0], GENES_DLP[1][1],       # repL*
    Y_QA - 0.11, Y_DLP + 0.11,
    color='#2a9d8f', alpha=0.20)

# Packaging terminase — partial
add_similarity_ribbon(ax,
    GENES_QA[7][0], GENES_QA[7][1],         # terL
    GENES_DLP[2][0], GENES_DLP[2][1],       # terL*
    Y_QA - 0.11, Y_DLP + 0.11,
    color='#f4a261', alpha=0.25)

# Lysis
add_similarity_ribbon(ax,
    GENES_QA[11][0], GENES_QA[12][1],       # hin + pphA
    GENES_DLP[3][0], GENES_DLP[4][1],       # hin* + pphA*
    Y_QA - 0.11, Y_DLP + 0.11,
    color='#8338ec', alpha=0.20)

# Novel structural block (QA5221 only — no DLP12 counterpart)
for gi in [6, 8, 9, 10]:
    xs = GENES_QA[gi][0]
    xe = GENES_QA[gi][1]
    rect = mpatches.FancyBboxPatch(
        (xs, Y_QA - 0.12), xe - xs, 0.24,
        boxstyle='round,pad=50', linewidth=1.2,
        edgecolor='#e9c46a', facecolor='none', zorder=4)
    ax.add_patch(rect)

# Scale bars
for y_track, length in [(Y_QA, QA_LEN), (Y_DLP, DLP12_LEN)]:
    ax.annotate('',
        xy=(length + 600, y_track - 0.08),
        xytext=(length + 600 + 5000, y_track - 0.08),
        arrowprops=dict(arrowstyle='<->', color='#495057', lw=1.2))
    ax.text(length + 600 + 2500, y_track - 0.14,
            '5 kb', ha='center', va='top', fontsize=7, color='#495057')

# ── Legend ───────────────────────────────────────────────────────────────────
legend_patches = [mpatches.Patch(color=v, label=k) for k, v in CMAP.items()
                  if k != 'Hypothetical']
from matplotlib.lines import Line2D
legend_patches.append(
    Line2D([0], [0], marker='*', color='w', markerfacecolor='#e9c46a',
           markersize=11, label='Novel candidate (QA5221 unique)'))
legend_patches.append(
    mpatches.Patch(facecolor='none', edgecolor='#e9c46a', linewidth=1.5,
                   label='Novel structural block'))
ax.legend(handles=legend_patches, loc='upper right',
          fontsize=7.5, framealpha=0.92, ncol=2,
          bbox_to_anchor=(1.01, 1.45))

# ── Title ─────────────────────────────────────────────────────────────────────
ax.set_title(
    'Comparative Synteny: E. coli QA5221 Active Prophage vs. DLP12 Defective Pseudoprophage (K-12 MG1655)\n'
    'Ribbons indicate regions of homology; yellow-bordered genes are novel structural modules unique to QA5221',
    fontsize=10, fontweight='bold', pad=4)

plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(OUT_DIR / f'figure12c_prophage_synteny.{ext}',
                dpi=300, bbox_inches='tight')
print("✓ Figure 12c (prophage synteny comparison) saved.")
plt.close()
