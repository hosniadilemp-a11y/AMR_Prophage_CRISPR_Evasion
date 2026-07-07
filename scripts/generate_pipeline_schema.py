#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, Circle, FancyArrowPatch
from matplotlib.path import Path
import matplotlib.patches as patches

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR   = os.path.dirname(SCRIPT_DIR)
OUTPUT_DIR = os.path.join(REPO_DIR, "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Set styling params
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['text.color'] = '#1E293B'

# Create canvas (tighter canvas size)
fig, ax = plt.subplots(figsize=(15, 6.0), facecolor='white')
ax.set_xlim(-0.6, 4.2)
ax.set_ylim(-0.65, 2.05)
ax.axis('off')

# Design colors
COLOR_STEP_BG = '#F8FAFC'
COLOR_BORDER = '#CBD5E1'
COLOR_ACCENT = '#2563EB'      # Blue
COLOR_INPUT = '#047857'       # Darker Emerald for readability
COLOR_OUTPUT = '#B91C1C'      # Darker Red for readability
COLOR_MUTED = '#64748B'

# Step-specific accent colors for headers
STEP_COLORS = {
    1: '#2563EB', # Blue
    2: '#7C3AED', # Purple
    3: '#059669', # Emerald
    4: '#D97706', # Amber
    5: '#DB2777', # Pink
    6: '#0891B2', # Cyan
    7: '#4F46E5', # Indigo
    8: '#DC2626'  # Red
}

# Coordinates of the steps in the 2x4 winding layout
# Row 1: 1 -> 2 -> 3 -> 4
# Row 2: 8 <- 7 <- 6 <- 5
coords = {
    1: (0.0, 1.4),
    2: (1.20, 1.4),
    3: (2.40, 1.4),
    4: (3.60, 1.4),
    5: (3.60, 0.0),
    6: (2.40, 0.0),
    7: (1.20, 0.0),
    8: (0.0, 0.0)
}

step_details = {
    1: {
        "title_short": "Genome Assembly & QC",
        "tools": "FastQC, SPAdes, mosdepth",
        "input": "Raw FASTQ Reads",
        "output": "Clean Contigs FASTA",
        "desc": "Trim reads, assemble contigs,\nverify continuity & coverage."
    },
    2: {
        "title_short": "Core Annotation",
        "tools": "Prokka, Bakta, VFDB",
        "input": "Clean Contigs FASTA",
        "output": "Hypothetical ORFs",
        "desc": "Predict coding ORFs;\nannotate known elements."
    },
    3: {
        "title_short": "Pangenome Matrix",
        "tools": "Panaroo, IQ-TREE2",
        "input": "400 public E. coli GFFs",
        "output": "QA5221 Cloud Singletons",
        "desc": "Cluster orthologous genes;\nidentify private singletons."
    },
    4: {
        "title_short": "Mobilome & Context",
        "tools": "ISEScan, clinker",
        "input": "Candidate flanking regions",
        "output": "Synteny & Context Table",
        "desc": "Analyze gene neighbourhoods\nfor transposases & MGEs."
    },
    5: {
        "title_short": "Compositional QC",
        "tools": "CodonW, SciPy (ANOVA)",
        "input": "Candidate sequences",
        "output": "Prioritized ORFs",
        "desc": "Filter candidates by length,\nGC% & GC3 codon bias."
    },
    6: {
        "title_short": "3D Protein Folding",
        "tools": "ESMFold, Foldseek, PyMOL",
        "input": "Protein sequences",
        "output": "3D Folds & TM-alignments",
        "desc": "Fold 3D structures;\nalign active site clefts."
    },
    7: {
        "title_short": "PLM Latent Proximity",
        "tools": "ESM-2 Embeddings, UMAP",
        "input": "Acquired AMR database",
        "output": "Functional Clustering",
        "desc": "Project sequence embeddings;\nconfirm functional clustering."
    },
    8: {
        "title_short": "Molecular Dynamics",
        "tools": "OpenMM (127 ns production)",
        "input": "3D Structures & Ligands",
        "output": "Backbone RMSD / RMSF",
        "desc": "Explicit solvent MD;\nverify catalytic rigidity."
    }
}

# Helper to draw a box with shadow and pills (more compact dimensions)
def draw_step_box(ax, x, y, step_num, info):
    w, h = 0.85, 1.08
    
    # 1. Draw Drop Shadow
    shadow = FancyBboxPatch((x - w/2 + 0.02, y - h/2 - 0.02), w, h, boxstyle="round,pad=0.03", 
                           linewidth=0, facecolor='#E2E8F0', zorder=1)
    ax.add_patch(shadow)
    
    # 2. Main card box
    box = FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle="round,pad=0.03", 
                         linewidth=1.0, edgecolor='#CBD5E1', facecolor='white', zorder=2)
    ax.add_patch(box)
    
    # 3. Step Header circle with number (top-left aligned)
    circle_x = x - w/2 + 0.08
    circle_y = y + h/2 - 0.08
    circle = Circle((circle_x, circle_y), 0.038, facecolor=STEP_COLORS[step_num], edgecolor='none', zorder=3)
    ax.add_patch(circle)
    ax.text(circle_x, circle_y, str(step_num), ha='center', va='center', 
            fontsize=7.0, fontweight='bold', color='white', zorder=4)
    
    # 4. Header title text next to circle
    ax.text(x - w/2 + 0.15, y + h/2 - 0.08, info["title_short"], ha='left', va='center', 
            fontsize=7.2, fontweight='bold', color='#0F172A', zorder=4)
    
    # 5. Separator line under header
    ax.plot([x - w/2 + 0.05, x + w/2 - 0.05], [y + h/2 - 0.16, y + h/2 - 0.16], color='#E2E8F0', lw=0.8, zorder=3)
    
    # 6. Description (centered)
    ax.text(x, y - 0.09, info["desc"], ha='center', va='center', fontsize=7.0, color='#334155', linespacing=1.2, zorder=3)
    
    # 7. Tools box
    ax.text(x, y - 0.23, f"Tools: {info['tools']}", ha='center', va='center', fontsize=6.8, fontstyle='italic', color=COLOR_MUTED, zorder=3)
    
    # 8. Input / Output Pills at the bottom (tighter positioning)
    # Input Pill
    ax.add_patch(FancyBboxPatch((x - 0.37, y - 0.38), 0.74, 0.06, boxstyle="round,pad=0.01", 
                                linewidth=0, facecolor='#ECFDF5', zorder=3))
    ax.text(x, y - 0.35, f"IN: {info['input']}", ha='center', va='center', 
            fontsize=6.2, fontweight='bold', color=COLOR_INPUT, zorder=4)
            
    # Output Pill
    ax.add_patch(FancyBboxPatch((x - 0.37, y - 0.47), 0.74, 0.06, boxstyle="round,pad=0.01", 
                                linewidth=0, facecolor='#FEF2F2', zorder=3))
    ax.text(x, y - 0.44, f"OUT: {info['output']}", ha='center', va='center', 
            fontsize=6.2, fontweight='bold', color=COLOR_OUTPUT, zorder=4)

# Draw custom vector pictograms inside each box
def draw_pictogram(ax, step_num, x, y):
    px, py = x, y + 0.14
    color = STEP_COLORS[step_num]
    
    if step_num == 1:
        # Step 1: Assembly (Contig line and reads mapping)
        ax.plot([px - 0.16, px + 0.16], [py, py], color='#475569', lw=2.2, zorder=4)
        ax.plot([px - 0.12, px - 0.02], [py + 0.04, py + 0.04], color=color, lw=1.2, zorder=4)
        ax.plot([px - 0.07, px + 0.05], [py + 0.08, py + 0.08], color=color, lw=1.2, zorder=4)
        ax.plot([px + 0.02, px + 0.12], [py + 0.04, py + 0.04], color=color, lw=1.2, zorder=4)
        
    elif step_num == 2:
        # Step 2: DNA double helix
        t = np.linspace(-np.pi, np.pi, 40)
        xs = px + 0.16 * t / np.pi
        ys1 = py + 0.04 * np.sin(2.5 * t)
        ys2 = py - 0.04 * np.sin(2.5 * t)
        ax.plot(xs, ys1, color=color, lw=1.2, zorder=4)
        ax.plot(xs, ys2, color='#A78BFA', lw=1.2, zorder=4)
        for i in range(0, len(t), 4):
            ax.plot([xs[i], xs[i]], [ys1[i], ys2[i]], color='#E2E8F0', lw=0.7, zorder=3)
            
    elif step_num == 3:
        # Step 3: Phylogenetic Cladogram Tree
        # Root
        ax.plot([px - 0.12, px - 0.07], [py, py], color='#475569', lw=1.2, zorder=4)
        # First branch split
        ax.plot([px - 0.07, px - 0.07], [py - 0.05, py + 0.05], color='#475569', lw=1.2, zorder=4)
        # Upper branch
        ax.plot([px - 0.07, px + 0.02], [py + 0.05, py + 0.05], color='#475569', lw=1.2, zorder=4)
        # Lower branch
        ax.plot([px - 0.07, px + 0.02], [py - 0.05, py - 0.05], color='#475569', lw=1.2, zorder=4)
        # Upper sub-branches
        ax.plot([px + 0.02, px + 0.02], [py + 0.02, py + 0.08], color='#475569', lw=1.2, zorder=4)
        ax.plot([px + 0.02, px + 0.12], [py + 0.08, py + 0.08], color=color, lw=1.2, zorder=4)
        ax.plot([px + 0.02, px + 0.12], [py + 0.02, py + 0.02], color=color, lw=1.2, zorder=4)
        # Lower branch end
        ax.plot([px + 0.02, px + 0.12], [py - 0.05, py - 0.05], color='#64748B', lw=1.2, zorder=4)
        
    elif step_num == 4:
        # Step 4: MGE flanking arrows
        def draw_arrow_icon(ax, x1, x2, y_c, fill_c):
            h_val = 0.04
            arrow = Polygon([(x1, y_c - h_val/2), (x2 - 0.03, y_c - h_val/2), (x2 - 0.03, y_c - h_val), (x2, y_c), (x2 - 0.03, y_c + h_val), (x2 - 0.03, y_c + h_val/2), (x1, y_c + h_val/2)], facecolor=fill_c, edgecolor='none', zorder=4)
            ax.add_patch(arrow)
            
        draw_arrow_icon(ax, px - 0.16, px - 0.05, py, '#94A3B8') 
        draw_arrow_icon(ax, px - 0.03, px + 0.05, py, color)     
        draw_arrow_icon(ax, px + 0.07, px + 0.16, py, '#3B82F6') 
        
    elif step_num == 5:
        # Step 5: Funnel or filtering diagram
        poly = Polygon([(px - 0.14, py + 0.07), (px + 0.14, py + 0.07), (px + 0.04, py - 0.07), (px - 0.04, py - 0.07)], facecolor='#FCE7F3', edgecolor=color, lw=1.0, zorder=4)
        ax.add_patch(poly)
        # Filtered lines inside
        ax.plot([px - 0.10, px + 0.10], [py + 0.02, py + 0.02], color='#F472B6', lw=0.7, ls='-', zorder=4)
        ax.plot([px - 0.06, px + 0.06], [py - 0.02, py - 0.02], color='#F472B6', lw=0.7, ls='-', zorder=4)
        # Drip arrow
        ax.plot([px, px], [py - 0.07, py - 0.12], color=color, lw=2.0, zorder=5)
        
    elif step_num == 6:
        # Step 6: 3D Protein fold spline ribbon
        t = np.linspace(0, 1, 100)
        x_curve = px + 0.14 * np.sin(2 * np.pi * t)
        y_curve = py + 0.06 * np.cos(3 * np.pi * t)
        ax.plot(x_curve, y_curve, color=color, lw=2.0, zorder=4, alpha=0.9)
        c_lig = Circle((px + 0.04, py - 0.02), 0.02, facecolor='#F59E0B', zorder=5, edgecolor='black', lw=0.4)
        ax.add_patch(c_lig)
        
    elif step_num == 7:
        # Step 7: PLM Latent space embeddings clusters
        np.random.seed(42)
        ax.scatter(px - 0.08 + 0.03*np.random.randn(4), py + 0.02 + 0.02*np.random.randn(4), s=12, color='#818CF8', edgecolors='none', zorder=4)
        ax.scatter(px + 0.08 + 0.03*np.random.randn(4), py - 0.02 + 0.02*np.random.randn(4), s=12, color='#22D3EE', edgecolors='none', zorder=4)
        ax.plot(px, py, marker='*', markersize=9, color=color, markeredgecolor='black', markeredgewidth=0.4, zorder=5)
        
    elif step_num == 8:
        # Step 8: MD simulation solvent box with water molecules
        box = plt.Rectangle((px - 0.14, py - 0.08), 0.28, 0.16, fill=True, facecolor='#FEF2F2', edgecolor=color, lw=0.8, zorder=3)
        ax.add_patch(box)
        for wx, wy in [(px-0.10, py-0.04), (px-0.07, py+0.04), (px+0.05, py-0.03), (px+0.10, py+0.04), (px+0.02, py+0.05), (px-0.02, py-0.05)]:
            ax.plot(wx, wy, marker='o', markersize=1.5, color='#FCA5A5', zorder=4)
        t = np.linspace(-0.10, 0.10, 50)
        y_fluc = py + 0.015 * np.sin(10 * t) + 0.004 * np.random.randn(len(t))
        ax.plot(px + t, y_fluc, color=color, lw=1.0, zorder=5)

# Draw connections between steps (adjusted for w=0.85, h=0.95)
def draw_connections(ax):
    w, h = 0.85, 0.95
    # Horizontal links (1->2->3->4)
    for s in [1, 2, 3]:
        x1, y1 = coords[s]
        x2, y2 = coords[s+1]
        arrow = FancyArrowPatch((x1 + (w/2 + 0.03), y1), (x2 - (w/2 + 0.03), y2),
                                arrowstyle='-|>', mutation_scale=12,
                                linewidth=2.0, color='#94A3B8', zorder=1)
        ax.add_patch(arrow)
        
    # Vertical link down (4->5)
    x1, y1 = coords[4]
    x2, y2 = coords[5]
    arrow = FancyArrowPatch((x1, y1 - (h/2 + 0.03)), (x2, y2 + (h/2 + 0.03)),
                            arrowstyle='-|>', mutation_scale=12,
                            linewidth=2.0, color='#94A3B8', zorder=1)
    ax.add_patch(arrow)
    
    # Horizontal links back (5->6->7->8)
    for s in [5, 6, 7]:
        x1, y1 = coords[s]
        x2, y2 = coords[s+1]
        arrow = FancyArrowPatch((x1 - (w/2 + 0.03), y1), (x2 + (w/2 + 0.03), y2),
                                arrowstyle='-|>', mutation_scale=12,
                                linewidth=2.0, color='#94A3B8', zorder=1)
        ax.add_patch(arrow)

# Execute drawing
for s in range(1, 9):
    cx, cy = coords[s]
    draw_step_box(ax, cx, cy, s, step_details[s])
    draw_pictogram(ax, s, cx, cy)

draw_connections(ax)

# Save figure in high-res PNG, vector PDF and SVG
png_out = os.path.join(OUTPUT_DIR, "step0_pipeline_schema.png")
pdf_out = os.path.join(OUTPUT_DIR, "step0_pipeline_schema.pdf")
svg_out = os.path.join(OUTPUT_DIR, "step0_pipeline_schema.svg")

plt.savefig(png_out, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(pdf_out, bbox_inches='tight', facecolor='white')
plt.savefig(svg_out, bbox_inches='tight', facecolor='white')
plt.close()

print("Pipeline schema regenerated successfully with compact drop-shadow card layout and no canvas titles!")
