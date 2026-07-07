#!/usr/bin/env python3
"""
Generate publication-quality prophage figures for paper2:
  Figure 12a: Prophage module map (linear genome map, NODE_1 56,422-98,747 bp)
  Figure 12b: Prophage gene inventory / functional module bar chart
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import matplotlib.patheffects as pe
import numpy as np
from pathlib import Path

# ── Output directory ────────────────────────────────────────────────────────
OUT_DIR = Path("/media/adel/Data/Hosni/openmm_windows_setup/AMR_Work/manuscript/paper2/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 9,
    'axes.linewidth': 0.8,
    'pdf.fonttype': 42,  # TrueType in PDF
    'ps.fonttype': 42,
})

# ── Prophage boundaries ──────────────────────────────────────────────────────
PROPH_START = 56422
PROPH_END   = 98747
PROPH_LEN   = PROPH_END - PROPH_START   # 42325 bp

# ── Colour palette (colourblind-friendly Mako-inspired) ─────────────────────
COLORS = {
    'Host Cargo':           '#e63946',   # red
    'Lysogeny/Regulation':  '#457b9d',   # steel blue
    'DNA Replication':      '#2a9d8f',   # teal
    'Packaging':            '#f4a261',   # orange
    'Structural Core':      '#264653',   # dark teal
    'Recombination/Lysis':  '#8338ec',   # violet
    'Hypothetical':         '#adb5bd',   # grey
    'attL/attR':            '#e9c46a',   # gold
}

# ── Gene data from Step12 report + paper2.tex table ────────────────────────
# (locus_tag, start, end, strand, category, short_label)
GENES = [
    # Host Cargo
    ('KNGPFPPJ_00060', 58216, 58404, '-', 'Host Cargo',           'rcbA'),
    ('KNGPFPPJ_00061', 58655, 59707, '-', 'Host Cargo',           'hlyE'),
    # Lysogeny/Regulation
    ('KNGPFPPJ_00065', 63467, 63688, '-', 'Lysogeny/Regulation',  'kilR'),
    ('KNGPFPPJ_00069', 65191, 65613, '+', 'Lysogeny/Regulation',  'ydaT'),
    # DNA Replication
    ('KNGPFPPJ_00071', 66504, 67250, '+', 'DNA Replication',      'dnaC'),
    ('KNGPFPPJ_00082', 72122, 72412, '+', 'DNA Replication',      'ybcO'),
    # Packaging
    ('KNGPFPPJ_00091', 76580, 77329, '+', 'Packaging',            'terL'),
    # Structural Core — 5 priority candidates
    ('KNGPFPPJ_00084', 73234, 73845, '+', 'Structural Core',      'antirep'),
    ('KNGPFPPJ_00107', 89827, 90477, '+', 'Structural Core',      'tail'),
    ('KNGPFPPJ_00109', 90786, 91820, '+', 'Structural Core',      'capsid'),
    ('KNGPFPPJ_00114', 93696, 94895, '+', 'Structural Core',      'baseplate'),
    # Other hypothetical structural (selected representative)
    ('KNGPFPPJ_00092', 77268, 78698, '+', 'Hypothetical',         ''),
    ('KNGPFPPJ_00093', 78698, 80164, '+', 'Hypothetical',         ''),
    ('KNGPFPPJ_00095', 80804, 82024, '+', 'Hypothetical',         ''),
    ('KNGPFPPJ_00097', 82544, 83485, '+', 'Hypothetical',         ''),
    ('KNGPFPPJ_00100', 84119, 84673, '+', 'Hypothetical',         ''),
    ('KNGPFPPJ_00103', 85591, 86736, '+', 'Hypothetical',         ''),
    ('KNGPFPPJ_00106', 87821, 89827, '+', 'Hypothetical',         ''),
    ('KNGPFPPJ_00112', 92591, 93343, '+', 'Hypothetical',         ''),
    ('KNGPFPPJ_00115', 94892, 95572, '+', 'Hypothetical',         ''),
    ('KNGPFPPJ_00116', 95572, 96228, '+', 'Hypothetical',         ''),
    # Recombination/Lysis
    ('KNGPFPPJ_00118', 96951, 97514, '+', 'Recombination/Lysis',  'hin'),
    ('KNGPFPPJ_00119', 98019, 98675, '+', 'Recombination/Lysis',  'pphA'),
]

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1 — Linear prophage module map
# ─────────────────────────────────────────────────────────────────────────────
def draw_arrow(ax, start, end, strand, y, h, color, label='', important=False):
    """Draw a gene arrow on y at height h."""
    xs = (start - PROPH_START) / PROPH_LEN
    xe = (end   - PROPH_START) / PROPH_LEN
    width = xe - xs
    head = min(0.012, width * 0.35)

    if strand == '+':
        ax.annotate('',
            xy=(xe, y), xytext=(xs, y),
            arrowprops=dict(arrowstyle=f'-|>, head_width={h*0.7}, head_length={head}',
                            color=color, lw=0, shrinkA=0, shrinkB=0),
            annotation_clip=False)
        ax.fill_between([xs, xe - head], [y - h/2]*2, [y + h/2]*2,
                         color=color, linewidth=0, zorder=3)
        # arrowhead triangle
        tri_x = [xe - head, xe - head, xe]
        tri_y = [y - h/2, y + h/2, y]
        ax.fill(tri_x, tri_y, color=color, linewidth=0, zorder=3)
    else:
        ax.fill_between([xs + head, xe], [y - h/2]*2, [y + h/2]*2,
                         color=color, linewidth=0, zorder=3)
        tri_x = [xs + head, xs + head, xs]
        tri_y = [y - h/2, y + h/2, y]
        ax.fill(tri_x, tri_y, color=color, linewidth=0, zorder=3)

    if label and important:
        mid_x = (xs + xe) / 2
        ax.text(mid_x, y + h / 2 + 0.025, label,
                ha='center', va='bottom', fontsize=7,
                fontweight='bold', color='#1a1a2e',
                path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])


fig, ax = plt.subplots(figsize=(14, 3.8))

# Backbone line
ax.plot([0, 1], [0.5, 0.5], color='#343a40', lw=2, zorder=1)

# attL / attR markers
for site_pos, site_name in [(PROPH_START, 'attL'), (PROPH_END, 'attR')]:
    xp = (site_pos - PROPH_START) / PROPH_LEN
    ax.axvline(xp, color=COLORS['attL/attR'], lw=2, ls='--', zorder=2, alpha=0.9)
    ax.text(xp, 0.08, site_name,
            ha='center', va='bottom', fontsize=8, color=COLORS['attL/attR'],
            fontweight='bold')

# Genomic scale ticks (every 5 kb)
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(0.0, 1.0)
tick_positions = range(0, PROPH_LEN + 1, 5000)
for tp in tick_positions:
    xp = tp / PROPH_LEN
    ax.axvline(xp, color='#dee2e6', lw=0.5, zorder=0)
    kb_label = f"{(tp + PROPH_START)//1000}kb"
    ax.text(xp, 0.03, kb_label,
            ha='center', va='bottom', fontsize=6.5, color='#6c757d')

# Module span shading
module_spans = [
    (PROPH_START, 62000,   'Host Cargo',          0.36),
    (63000, 68000,         'Lysogeny/Regulation',  0.36),
    (66000, 72500,         'DNA Replication',       0.36),
    (73000, 78000,         'Packaging',             0.36),
    (73000, 97000,         'Structural Core',       0.15),
    (96800, 98750,         'Recombination/Lysis',  0.36),
]
for ms, me, mcat, alpha in module_spans:
    xs = (ms - PROPH_START) / PROPH_LEN
    xe = (me - PROPH_START) / PROPH_LEN
    ax.axvspan(xs, xe, ymin=0.15, ymax=0.85,
               alpha=alpha * 0.35, color=COLORS[mcat], zorder=0)

# Draw genes
GENE_Y   = 0.50   # centre of backbone
GENE_H   = 0.16
PRIORITY_TAGS = {'KNGPFPPJ_00084', 'KNGPFPPJ_00107', 'KNGPFPPJ_00109',
                 'KNGPFPPJ_00114', 'KNGPFPPJ_00091'}

for tag, s, e, strand, cat, lbl in GENES:
    important = (tag in PRIORITY_TAGS) or (cat == 'Host Cargo') or bool(lbl and cat != 'Hypothetical')
    draw_arrow(ax, s, e, strand, GENE_Y, GENE_H,
               color=COLORS[cat], label=lbl, important=important)

# Priority candidate star markers
for tag, s, e, strand, cat, lbl in GENES:
    if tag in PRIORITY_TAGS:
        mid_x = ((s + e) / 2 - PROPH_START) / PROPH_LEN
        ax.plot(mid_x, GENE_Y + GENE_H / 2 + 0.10, '*',
                color='#e9c46a', markersize=9, zorder=5)

# Legend
legend_patches = [
    mpatches.Patch(color=v, label=k) for k, v in COLORS.items()
    if k not in ('attL/attR',)
]
# add a star for candidates
from matplotlib.lines import Line2D
legend_patches.append(
    Line2D([0], [0], marker='*', color='w', markerfacecolor='#e9c46a',
           markersize=10, label='Priority candidate'))
ax.legend(handles=legend_patches, loc='upper right',
          fontsize=7, framealpha=0.9, ncol=2,
          bbox_to_anchor=(1.0, 1.0))

ax.set_xlabel(f'Coordinate on NODE_1  (attL: {PROPH_START:,} bp  →  attR: {PROPH_END:,} bp,  Δ {PROPH_LEN:,} bp)',
              fontsize=8.5)
ax.set_yticks([])
ax.set_xticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

ax.set_title('Prophage Boundary Map and Functional Module Organisation (NODE_1, QA5221)',
             fontsize=10, fontweight='bold', pad=8)

plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(OUT_DIR / f'figure12a_prophage_module_map.{ext}',
                dpi=300, bbox_inches='tight')
print("✓ Figure 12a (prophage module map) saved.")
plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2 — Prophage functional module composition (stacked bar / donut)
# ─────────────────────────────────────────────────────────────────────────────
# Count genes per category
from collections import Counter
cat_counts = Counter(g[4] for g in GENES)

# Order for display
CAT_ORDER = [
    'Host Cargo', 'Lysogeny/Regulation', 'DNA Replication',
    'Packaging', 'Structural Core', 'Recombination/Lysis', 'Hypothetical'
]
counts = [cat_counts.get(c, 0) for c in CAT_ORDER]
colors = [COLORS[c] for c in CAT_ORDER]

fig2, (ax_donut, ax_bar) = plt.subplots(1, 2, figsize=(12, 5))

# ── Donut chart ──────────────────────────────────────────────────────────────
total_genes = sum(counts)
wedges, texts, autotexts = ax_donut.pie(
    counts,
    labels=None,
    colors=colors,
    autopct=lambda p: f'{p:.1f}%' if p > 3 else '',
    startangle=90,
    wedgeprops=dict(width=0.55, edgecolor='white', linewidth=1.5),
    pctdistance=0.80,
    textprops={'fontsize': 8}
)
for at in autotexts:
    at.set_fontweight('bold')
    at.set_color('white')

# Centre text
ax_donut.text(0, 0.05, str(total_genes), ha='center', va='center',
              fontsize=22, fontweight='bold', color='#1a1a2e')
ax_donut.text(0, -0.20, 'CDS in\nprophage', ha='center', va='center',
              fontsize=9, color='#495057')
ax_donut.set_title('Prophage CDS Functional Composition\n(62 coding sequences, 42.3 kb)',
                   fontsize=10, fontweight='bold', pad=10)

# Legend
patches = [mpatches.Patch(color=c, label=f'{n}  {l}')
           for c, l, n in zip(colors, CAT_ORDER, counts)]
ax_donut.legend(handles=patches, loc='lower left',
                bbox_to_anchor=(-0.35, -0.15), fontsize=8, framealpha=0.9)

# ── Horizontal bar chart ─────────────────────────────────────────────────────
bars = ax_bar.barh(CAT_ORDER[::-1], counts[::-1],
                   color=colors[::-1], edgecolor='white', linewidth=0.8)
for bar, cnt in zip(bars, counts[::-1]):
    ax_bar.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2,
                str(cnt), va='center', ha='left', fontsize=9, fontweight='bold',
                color='#1a1a2e')

ax_bar.set_xlabel('Number of CDS', fontsize=9)
ax_bar.set_title('CDS Count by Functional Module', fontsize=10, fontweight='bold')
ax_bar.spines[['top', 'right']].set_visible(False)
ax_bar.tick_params(axis='y', labelsize=8.5)
ax_bar.set_xlim(0, max(counts) + 2.5)

plt.tight_layout()
for ext in ('pdf', 'png'):
    fig2.savefig(OUT_DIR / f'figure12b_prophage_functional_composition.{ext}',
                 dpi=300, bbox_inches='tight')
print("✓ Figure 12b (functional composition) saved.")
plt.close()

print("\nAll prophage figures generated successfully.")
