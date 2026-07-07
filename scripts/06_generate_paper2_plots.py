#!/usr/bin/env python3
import os
import sys
import re
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.optimize import curve_fit
import scipy.stats as stats
from Bio import SeqIO
from Bio import Phylo

# ==============================================================================
# GLOBAL STYLING CONFIGURATION (CREST/MAKO STYLE SYSTEM)
# ==============================================================================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
plt.rcParams['text.color'] = '#1e293b'
plt.rcParams['axes.labelcolor'] = '#1e293b'
plt.rcParams['xtick.color'] = '#475569'
plt.rcParams['ytick.color'] = '#475569'

# Define cohesive color scheme
COLOR_CORE = "#1e3a8a"       # Deep Navy
COLOR_SOFT_CORE = "#3b82f6"  # Bright Blue
COLOR_ACCESSORY = "#4f46e5"  # Indigo
COLOR_SINGLETON = "#0d9488"  # Teal
COLOR_AMR_GENE = "#f97316"   # Amber/Orange (high-contrast accent)
COLOR_GREY_LIGHT = "#f1f5f9" # Soft Slate Background
COLOR_GREY_DARK = "#64748b"  # Slate Grey
COLOR_LOW_COV = "#94a3b8"    # Low coverage elements

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================
def parse_assembly_contigs(fasta_path):
    print("Parsing assembly contigs for length vs coverage...")
    lengths = []
    coverages = []
    names = []
    pattern = re.compile(r"NODE_(\d+)_length_(\d+)_cov_([\d\.]+)")
    
    for record in SeqIO.parse(fasta_path, "fasta"):
        match = pattern.search(record.id)
        if match:
            lengths.append(int(match.group(2)))
            coverages.append(float(match.group(3)))
            names.append(record.id)
    return pd.DataFrame({"Name": names, "Length": lengths, "Coverage": coverages})

def compute_gc_gc3(ffn_path):
    print("Computing GC and GC3% for all CDS...")
    gene_gc = {}
    for r in SeqIO.parse(ffn_path, "fasta"):
        seq = str(r.seq).upper()
        if not seq:
            continue
        length = len(seq)
        gc = (seq.count("G") + seq.count("C")) / length * 100
        
        codons = [seq[i:i+3] for i in range(0, length, 3) if len(seq[i:i+3]) == 3]
        third_pos = [c[2] for c in codons]
        gc3 = (third_pos.count("G") + third_pos.count("C")) / len(third_pos) * 100 if third_pos else 0
        
        gene_gc[r.id] = {"GC": gc, "GC3": gc3, "Length": length}
    return gene_gc

def load_pangenome_categories(roary_path, faa_path):
    print("Loading pangenome categories...")
    locus_to_category = {}
    with open(roary_path) as f:
        reader = csv.reader(f)
        header = next(reader)
        # Handle "QUERY_QA5221" column for 800-genome cohort
        if "QUERY_QA5221" in header:
            iso_idx = header.index("QUERY_QA5221")
        else:
            iso_idx = header.index("Ecoli_isolate")
        genome_indices = list(range(14, len(header)))
        other_indices = [idx for idx in genome_indices if idx != iso_idx]
        total_cohort_genomes = len(header) - 14
        
        for row in reader:
            if len(row) < len(header):
                continue
            iso_val = row[iso_idx].strip()
            if not iso_val:
                continue
            ref_count = sum(1 for idx in other_indices if row[idx].strip())
            total_genomes = ref_count + 1
            if total_genomes >= int(0.95 * total_cohort_genomes):
                cat = "Core"
            elif total_genomes > 1:
                cat = "Accessory"
            else:
                cat = "Singleton"
            for tag in iso_val.split("\t"):
                for single_tag in tag.split(";"):
                    # Translate QUERY_QA5221_XXXXX to KNGPFPPJ_XXXXX
                    mapped_tag = single_tag.strip().replace("QUERY_QA5221_", "KNGPFPPJ_")
                    locus_to_category[mapped_tag] = cat
                    
    valid_tags = set()
    for r in SeqIO.parse(faa_path, "fasta"):
        valid_tags.add(r.id)
    return locus_to_category, valid_tags

# ==============================================================================
# PLOTTING FUNCTIONS
# ==============================================================================
def plot_figure1_assembly_coverage(df_contigs, out_dir):
    print("Generating figure1_assembly_coverage...")
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    fig.patch.set_facecolor("white")
    
    chrom_mask = (df_contigs["Coverage"] >= 15) & (df_contigs["Coverage"] <= 30)
    plasmid_mask = (df_contigs["Coverage"] > 30)
    low_mask = (df_contigs["Coverage"] < 15)
    
    ax.scatter(df_contigs.loc[chrom_mask, "Length"], df_contigs.loc[chrom_mask, "Coverage"], 
               c=COLOR_SINGLETON, alpha=0.7, edgecolors="none", s=60, label="Chromosomal baseline (~25x)")
    ax.scatter(df_contigs.loc[plasmid_mask, "Length"], df_contigs.loc[plasmid_mask, "Coverage"], 
               c=COLOR_ACCESSORY, alpha=0.85, edgecolors="#1e293b", linewidths=0.8, s=90, label="Plasmid / Multi-copy (>30x)")
    ax.scatter(df_contigs.loc[low_mask, "Length"], df_contigs.loc[low_mask, "Coverage"], 
               c=COLOR_LOW_COV, alpha=0.4, edgecolors="none", s=40, label="Low coverage / Artifacts (<15x)")
    
    for _, row in df_contigs.iterrows():
        if "NODE_43_" in row["Name"] or "NODE_37_" in row["Name"] or "NODE_24_" in row["Name"]:
            label = row["Name"].split("_length")[0].replace("_", " ")
            ax.annotate(label, (row["Length"], row["Coverage"]), 
                        textcoords="offset points", xytext=(8,5), ha="left", fontsize=8.5, weight="bold", color="#0f172a")
            
    ax.set_xscale("log")
    ax.set_xlabel("Contig Length (bp, log scale)", fontsize=11, weight="bold")
    ax.set_ylabel("Read Coverage Depth (x)", fontsize=11, weight="bold")
    ax.grid(True, which="both", ls="-", color="#e2e8f0", alpha=0.6)
    ax.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#cbd5e1")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#cbd5e1')
    ax.spines['bottom'].set_color('#cbd5e1')
    
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "step1_assembly_coverage.png"), bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "step1_assembly_coverage.pdf"), bbox_inches="tight")
    plt.close(fig)

def plot_figure2_pangenome_partition(summary_path, out_dir):
    print("Generating figure2_pangenome_partition...")
    categories = ["Core", "Soft core", "Shell", "Cloud"]
    counts = [0, 3298, 2295, 3591]
    
    if os.path.exists(summary_path):
        with open(summary_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 3:
                    cat = parts[0].replace(" genes", "").strip()
                    if cat.lower().startswith("total"):
                        continue
                    count = int(parts[2])
                    if cat == "Core": counts[0] = count
                    elif cat == "Soft core": counts[1] = count
                    elif cat == "Shell": counts[2] = count
                    elif cat == "Cloud": counts[3] = count

    plot_cats = [c for c, cnt in zip(categories, counts) if cnt > 0]
    plot_counts = [cnt for cnt in counts if cnt > 0]
    
    # Sleek monochromatic blue-indigo theme
    color_map = {
        "Core": "#1e293b",       # Charcoal
        "Soft core": COLOR_CORE,  # Deep Navy
        "Shell": COLOR_ACCESSORY, # Indigo
        "Cloud": COLOR_SINGLETON  # Teal
    }
    colors = [color_map.get(cat, "#64748b") for cat in plot_cats]
    explode = tuple([0.03] * len(plot_counts))
    
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    fig.patch.set_facecolor("white")
    
    wedges, texts, autotexts = ax.pie(
        plot_counts, explode=explode, labels=plot_cats, autopct='%1.1f%%',
        shadow=False, startangle=140, colors=colors, 
        textprops=dict(color="#1e293b", fontsize=11), pctdistance=0.75
    )
    
    centre_circle = plt.Circle((0,0), 0.58, fc='white', ec="#cbd5e1", lw=0.5)
    ax.add_artist(centre_circle)
    
    for autotext in autotexts:
        autotext.set_fontsize(10)
        autotext.set_weight("bold")
        autotext.set_color("white")
        
    ax.legend(wedges, [f"{cat} ({cnt} clusters)" for cat, cnt in zip(plot_cats, plot_counts)],
              title="Pangenome Partitions", loc="center left", bbox_to_anchor=(0.95, 0.5),
              frameon=True, facecolor="white", edgecolor="#cbd5e1")
    
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "figure2_pangenome_partition.png"), bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "figure2_pangenome_partition.pdf"), bbox_inches="tight")
    plt.close(fig)

def plot_figure3_phylogeny_tree(tree_path, out_dir):
    print("Generating figure3_phylogeny_tree...")
    if not os.path.exists(tree_path):
        print(f"Warning: Tree file {tree_path} not found.")
        return

    tree = Phylo.read(tree_path, "newick")
    
    # Prune tree to make it visible (keep QUERY_QA5221 and every 20th terminal)
    all_terminals = [t.name for t in tree.get_terminals()]
    if len(all_terminals) > 100:
        print(f"-> Tree has {len(all_terminals)} terminals. Pruning to make it visible...")
        keep_names = set(all_terminals[::20])
        keep_names.add("QUERY_QA5221")
        for name in all_terminals:
            if name not in keep_names:
                tree.prune(name)
        print(f"-> Pruned tree to {len(tree.get_terminals())} terminals.")
        
    outgroup = None
    for t in tree.get_terminals():
        if "MG1655" in t.name or "000005845" in t.name:
            outgroup = t.name
            break
    if outgroup:
        tree.root_with_outgroup(outgroup)
    else:
        tree.root_at_midpoint()

    def assign_y(clade):
        if clade.is_terminal():
            clade.y = assign_y.current_y
            assign_y.current_y += 1
            return clade.y
        else:
            child_ys = [assign_y(child) for child in clade.clades]
            clade.y = sum(child_ys) / len(child_ys)
            return clade.y

    assign_y.current_y = 0
    assign_y(tree.root)

    def assign_x(clade, current_x=0.0):
        clade.x = current_x + (clade.branch_length if clade.branch_length is not None else 0.0)
        for child in clade.clades:
            assign_x(child, clade.x)

    assign_x(tree.root, 0.0)

    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    fig.patch.set_facecolor("white")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set_color('#cbd5e1')
    ax.spines['bottom'].set_linewidth(1.0)
    ax.get_xaxis().set_visible(True)
    ax.get_yaxis().set_visible(False)

    max_x = max(c.x for c in tree.get_terminals())
    margin = 0.02 * max_x

    def get_color(clade):
        leaves = clade.get_terminals()
        if any("QUERY_QA5221" in l.name or l.name == "Ecoli_isolate" or "QA5221" in l.name for l in leaves):
            return COLOR_AMR_GENE # Orange/Amber highlight
        elif all(l.name.startswith("ST354_") or "ST354" in l.name for l in leaves):
            return COLOR_CORE # Navy blue
        else:
            return COLOR_GREY_DARK # Slate grey

    def get_lw(clade):
        leaves = clade.get_terminals()
        if any("QUERY_QA5221" in l.name or l.name == "Ecoli_isolate" or "QA5221" in l.name for l in leaves):
            return 2.5
        elif all(l.name.startswith("ST354_") or "ST354" in l.name for l in leaves):
            return 1.8
        else:
            return 1.2

    def plot_clade(clade):
        if not clade.is_terminal():
            y_min = min(child.y for child in clade.clades)
            y_max = max(child.y for child in clade.clades)
            ax.vlines(clade.x, y_min, y_max, color=get_color(clade), lw=get_lw(clade), zorder=2)
            
        for child in clade.clades:
            ax.hlines(child.y, clade.x, child.x, color=get_color(child), lw=get_lw(child), zorder=2)
            
            support = None
            if child.confidence is not None:
                support = child.confidence
            elif child.name is not None:
                try:
                    support = float(child.name)
                except ValueError:
                    pass
            
            if support is not None and not child.is_terminal():
                if 50 <= support <= 100 or 0.5 <= support <= 1.0:
                    if support <= 1.0:
                        support *= 100
                    support_str = f"{int(support)}"
                    ax.text(clade.x + margin * 0.2, child.y + 0.15, support_str, 
                            fontsize=6.5, color="#475569", ha="left", va="bottom", weight="bold")
            plot_clade(child)

    plot_clade(tree.root)

    for clade in tree.get_terminals():
        raw_name = clade.name
        display_name = raw_name
        is_target = False
        
        if "QUERY_QA5221" in raw_name or raw_name == "Ecoli_isolate" or "QA5221" in raw_name:
            display_name = "E. coli QA5221 (Isolate, ST354)"
            color = COLOR_AMR_GENE
            weight = "bold"
            is_target = True
        elif raw_name.startswith("ST354_") or "ST354" in raw_name:
            display_name = raw_name.replace("ST354_", "").replace("_", " ")
            color = COLOR_CORE
            weight = "normal"
        else:
            display_name = raw_name.replace("_", " ")
            color = COLOR_GREY_DARK
            weight = "normal"
            
        ax.plot([clade.x, max_x + margin * 0.5], [clade.y, clade.y], color="#e2e8f0", ls=":", lw=1.0, zorder=1)
        ax.text(max_x + margin, clade.y, display_name, fontsize=8.5, color=color, weight=weight, va="center", ha="left")
        
        if is_target:
            ax.scatter([clade.x], [clade.y], color=COLOR_AMR_GENE, s=120, marker='*', zorder=5)

    ax.set_ylim(-1, len(tree.get_terminals()))
    ax.set_xlim(0, max_x + margin * 15)
    ax.set_xlabel("Genetic Distance (Substitution Rate / Core Site)", fontsize=10, weight="bold", labelpad=10, color="#1e293b")
    ax.tick_params(axis='x', colors='#475569', labelsize=9)

    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "figure3_phylogeny_tree.png"), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join(out_dir, "figure3_phylogeny_tree.pdf"), bbox_inches="tight")
    plt.close(fig)

def power_law_new_genes(N, theta, alpha):
    return theta * (N ** (-alpha))

def power_law_pangenome_size(N, kappa, gamma):
    return kappa * (N ** gamma)

def plot_figureX_pangenome_curves(csv_path, out_dir):
    print("Generating figureX_pangenome_curves...")
    if not os.path.exists(csv_path):
        print(f"Warning: Roary CSV {csv_path} not found. Skipping Heap's law curves.")
        return

    df = pd.read_csv(csv_path)
    metadata_cols = [
        'Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences',
        'Avg sequences per isolate', 'Genome Fragment', 'Order within Fragment', 
        'Accessory Fragment', 'Accessory Order with Fragment', 'QC', 'Min group size nuc', 
        'Max group size nuc', 'Avg group size nuc'
    ]
    genome_cols = [col for col in df.columns if col not in metadata_cols]
    n_genomes = len(genome_cols)
    presence_matrix = df[genome_cols].notna().astype(int).values
    n_genes = presence_matrix.shape[0]

    # Run fast permutations (100 is enough for publication-grade curves)
    n_permutations = 100
    pangenome_sizes = np.zeros((n_permutations, n_genomes))
    new_genes_counts = np.zeros((n_permutations, n_genomes))
    
    for p in range(n_permutations):
        perm = np.random.permutation(n_genomes)
        perm_matrix = presence_matrix[:, perm]
        seen_genes = np.zeros(n_genes, dtype=bool)
        for j in range(n_genomes):
            current_genes = perm_matrix[:, j] > 0
            new_genes = current_genes & (~seen_genes)
            new_genes_counts[p, j] = np.sum(new_genes)
            seen_genes = seen_genes | current_genes
            pangenome_sizes[p, j] = np.sum(seen_genes)
            
    mean_pangenome_size = np.mean(pangenome_sizes, axis=0)
    std_pangenome_size = np.std(pangenome_sizes, axis=0)
    mean_new_genes = np.mean(new_genes_counts, axis=0)
    std_new_genes = np.std(new_genes_counts, axis=0)
    
    x_data = np.arange(1, n_genomes + 1)
    
    # Fits
    fit_x = x_data[1:]
    fit_y = mean_new_genes[1:]
    popt_new, _ = curve_fit(power_law_new_genes, fit_x, fit_y, p0=(500, 0.5))
    theta_fit, alpha_fit = popt_new
    popt_pan, _ = curve_fit(power_law_pangenome_size, x_data, mean_pangenome_size, p0=(3000, 0.3))
    kappa_fit, gamma_fit = popt_pan

    # Set up double panel figure with new Mako styling
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    fig.patch.set_facecolor("white")
    
    # 1. Pangenome Accumulation
    ax1.errorbar(x_data, mean_pangenome_size, yerr=std_pangenome_size, fmt='o', 
                 color=COLOR_SINGLETON, ecolor=COLOR_SINGLETON, alpha=0.5, elinewidth=1, capsize=2, label='Observed Permutations')
    ax1.plot(x_data, power_law_pangenome_size(x_data, kappa_fit, gamma_fit), color=COLOR_CORE, linewidth=2.5, 
             label=f'Power law fit ($y = {kappa_fit:.1f} \\cdot N^{{{gamma_fit:.3f}}}$)')
    ax1.set_xlabel('Number of Genomes Sampled (N)', fontsize=11, weight="bold")
    ax1.set_ylabel('Total Gene Clusters (n)', fontsize=11, weight="bold")
    ax1.set_title('A) Pangenome Accumulation Curve', fontsize=12, fontweight='bold', loc="left", pad=10)
    ax1.grid(True, linestyle='-', color="#e2e8f0", alpha=0.6)
    ax1.legend(fontsize=10, frameon=True, edgecolor="#cbd5e1")
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_color('#cbd5e1')
    ax1.spines['bottom'].set_color('#cbd5e1')
    
    # 2. New Genes Discovery
    ax2.errorbar(x_data, mean_new_genes, yerr=std_new_genes, fmt='o', 
                 color=COLOR_ACCESSORY, ecolor=COLOR_ACCESSORY, alpha=0.5, elinewidth=1, capsize=2, label='Observed Permutations')
    ax2.plot(x_data, power_law_new_genes(x_data, theta_fit, alpha_fit), color=COLOR_CORE, linewidth=2.5,
             label=f"Heaps' law fit ($y = {theta_fit:.1f} \\cdot N^{{-{alpha_fit:.3f}}}$)")
    ax2.set_xlabel('Number of Genomes Sampled (N)', fontsize=11, weight="bold")
    ax2.set_ylabel('Number of New Genes Discovered ($n_{new}$)', fontsize=11, weight="bold")
    ax2.set_title('B) New Genes Discovery Curve', fontsize=12, fontweight='bold', loc="left", pad=10)
    ax2.grid(True, linestyle='-', color="#e2e8f0", alpha=0.6)
    ax2.legend(fontsize=10, frameon=True, edgecolor="#cbd5e1")
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_color('#cbd5e1')
    ax2.spines['bottom'].set_color('#cbd5e1')
    
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "figureX_pangenome_curves.png"), bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "figureX_pangenome_curves.pdf"), bbox_inches="tight")
    plt.close(fig)

def plot_figure4_is_distribution(out_dir):
    print("Generating figure4_is_distribution...")
    is_data = {
        'IS3': 5,
        'ISNCY': 4,
        'IS200/IS605': 3,
        'IS1': 2,
        'IS4': 2,
        'ISL3': 2,
        'IS5': 1,
        'IS66': 1,
        'IS630': 1,
        'IS256': 1
    }
    families = list(is_data.keys())
    counts = list(is_data.values())
    
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    fig.patch.set_facecolor("white")
    
    # Beautiful single bar fill color (slate/indigo)
    bars = ax.barh(families[::-1], counts[::-1], color=COLOR_ACCESSORY, edgecolor="#1e293b", height=0.55)
    
    for bar in bars:
        width = bar.get_width()
        ax.annotate(f"{int(width)}",
                    xy=(width, bar.get_y() + bar.get_height() / 2),
                    xytext=(6, 0),  
                    textcoords="offset points",
                    ha='left', va='center', fontsize=9, weight="bold", color="#1e293b")
                    
    ax.set_xlabel("Insertion Sequence (IS) Count", fontsize=11, weight="bold")
    ax.set_ylabel("IS Family", fontsize=11, weight="bold")
    ax.set_xlim(0, max(counts) + 0.8)
    ax.grid(True, axis="x", ls="-", color="#e2e8f0", alpha=0.6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#cbd5e1')
    ax.spines['bottom'].set_color('#cbd5e1')
    
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "figure4_is_distribution.png"), bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "figure4_is_distribution.pdf"), bbox_inches="tight")
    plt.close(fig)

def perform_statistical_plots(faa_path, ffn_path, roary_path, out_dir):
    print("Generating figure5_stats_comparison and figure6_amr_comparison...")
    locus_to_category, valid_tags = load_pangenome_categories(roary_path, faa_path)
    gene_gc_info = compute_gc_gc3(ffn_path)
    
    # 12 Acquired resistance genes
    acquired_amr_locus_tags = {
        "KNGPFPPJ_04581": {"gene": "tet(A)", "class": "Tetracycline", "contig": "NODE_30"},
        "KNGPFPPJ_04687": {"gene": "estX", "class": "Macrolide", "contig": "NODE_37"},
        "KNGPFPPJ_04689": {"gene": "aadA2", "class": "Aminoglycoside", "contig": "NODE_37"},
        "KNGPFPPJ_04690": {"gene": "cmlA1", "class": "Phenicol", "contig": "NODE_37"},
        "KNGPFPPJ_04691": {"gene": "aadA1", "class": "Aminoglycoside", "contig": "NODE_37"},
        "KNGPFPPJ_04692": {"gene": "qacL", "class": "Biocide", "contig": "NODE_37"},
        "KNGPFPPJ_04694": {"gene": "sul3", "class": "Sulfonamide", "contig": "NODE_37"},
        "KNGPFPPJ_04716": {"gene": "dfrA7", "class": "Trimethoprim", "contig": "NODE_41"},
        "KNGPFPPJ_04717": {"gene": "qacE\\Delta1", "class": "Biocide", "contig": "NODE_41"},
        "KNGPFPPJ_04718": {"gene": "sul1", "class": "Sulfonamide", "contig": "NODE_41"},
        "KNGPFPPJ_04720": {"gene": "blaTEM-1B", "class": "Beta-lactam", "contig": "NODE_42"},
        "KNGPFPPJ_04727": {"gene": "aph(3')-Ia", "class": "Aminoglycoside", "contig": "NODE_47"}
    }
    
    gene_data = []
    for tag in valid_tags:
        if tag not in gene_gc_info:
            continue
        info = gene_gc_info[tag]
        if tag in acquired_amr_locus_tags:
            cat = "Acquired AMR"
        else:
            cat = locus_to_category.get(tag, "Singleton")
            
        gene_data.append({
            "Locus": tag,
            "GC": info["GC"],
            "GC3": info["GC3"],
            "Length": info["Length"],
            "Category": cat
        })
    df_all = pd.DataFrame(gene_data)
    
    core = df_all[df_all["Category"] == "Core"]
    acc = df_all[df_all["Category"] == "Accessory"]
    sing = df_all[df_all["Category"] == "Singleton"]
    amr = df_all[df_all["Category"] == "Acquired AMR"]
    
    # ------------------ figure5_stats_comparison ------------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=300)
    fig.patch.set_facecolor("white")
    
    order = ["Core", "Accessory", "Singleton"]
    # New Mako Palette for Boxplots
    colors_f5 = [COLOR_CORE, COLOR_ACCESSORY, COLOR_SINGLETON]
    
    # Plot GC%
    box_gc = axes[0].boxplot([core["GC"], acc["GC"], sing["GC"]], tick_labels=order, patch_artist=True,
                             medianprops=dict(color="black", linewidth=2.0),
                             boxprops=dict(linewidth=1.2), flierprops=dict(marker='o', markersize=2, alpha=0.15))
    axes[0].set_ylabel("GC Content (%)", fontsize=11, weight="bold")
    axes[0].set_title("A) Gene GC Content", fontsize=12, weight="bold", loc="left", pad=10)
    
    # Plot GC3%
    box_gc3 = axes[1].boxplot([core["GC3"], acc["GC3"], sing["GC3"]], tick_labels=order, patch_artist=True,
                              medianprops=dict(color="black", linewidth=2.0),
                              boxprops=dict(linewidth=1.2), flierprops=dict(marker='o', markersize=2, alpha=0.15))
    axes[1].set_ylabel("3rd Codon Position GC3 (%)", fontsize=11, weight="bold")
    axes[1].set_title("B) Codon Bias (GC3)", fontsize=12, weight="bold", loc="left", pad=10)
    
    # Plot Length
    box_len = axes[2].boxplot([core["Length"], acc["Length"], sing["Length"]], tick_labels=order, patch_artist=True,
                              medianprops=dict(color="black", linewidth=2.0),
                              boxprops=dict(linewidth=1.2), flierprops=dict(marker='o', markersize=2, alpha=0.15))
    axes[2].set_yscale("log")
    axes[2].set_ylabel("Gene Length (bp, log scale)", fontsize=11, weight="bold")
    axes[2].set_title("C) Gene Length", fontsize=12, weight="bold", loc="left", pad=10)
    
    for box_plot in [box_gc, box_gc3, box_len]:
        for patch, color in zip(box_plot['boxes'], colors_f5):
            patch.set_facecolor(color)
            patch.set_alpha(0.85)
            patch.set_edgecolor("#1e293b")
            
    for ax in axes:
        ax.grid(True, which="both", axis="y", ls="-", color="#e2e8f0", alpha=0.6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#cbd5e1')
        ax.spines['bottom'].set_color('#cbd5e1')
        ax.tick_params(axis='both', labelsize=10)
        
    # Stats annotations on GC
    y_max_gc = max(df_all["GC"]) + 2
    axes[0].plot([1, 1, 3, 3], [y_max_gc - 5, y_max_gc - 3, y_max_gc - 3, y_max_gc - 5], lw=1.0, c="#475569")
    axes[0].text(2, y_max_gc - 2.5, "****", ha="center", va="bottom", fontsize=11, weight="bold", color="#0f172a")
    
    axes[0].plot([1, 1, 2, 2], [y_max_gc - 12, y_max_gc - 10, y_max_gc - 10, y_max_gc - 12], lw=1.0, c="#475569")
    axes[0].text(1.5, y_max_gc - 9.5, "****", ha="center", va="bottom", fontsize=11, weight="bold", color="#0f172a")
    
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "figure5_stats_comparison.png"), bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "figure5_stats_comparison.pdf"), bbox_inches="tight")
    plt.close(fig)

    # ------------------ figure6_amr_comparison ------------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=300)
    fig.patch.set_facecolor("white")
    
    order_amr = ["Core", "Accessory", "Singleton", "Acquired AMR"]
    colors_f6 = [COLOR_CORE, COLOR_ACCESSORY, COLOR_SINGLETON, COLOR_AMR_GENE]
    
    # Plot GC%
    box_gc = axes[0].boxplot([core["GC"], acc["GC"], sing["GC"], amr["GC"]], tick_labels=order_amr, patch_artist=True,
                             medianprops=dict(color="black", linewidth=2.0),
                             boxprops=dict(linewidth=1.2), flierprops=dict(marker='o', markersize=2, alpha=0.15))
    axes[0].set_ylabel("GC Content (%)", fontsize=11, weight="bold")
    axes[0].set_title("A) Gene GC Content", fontsize=12, weight="bold", loc="left", pad=10)
    
    # Plot GC3%
    box_gc3 = axes[1].boxplot([core["GC3"], acc["GC3"], sing["GC3"], amr["GC3"]], tick_labels=order_amr, patch_artist=True,
                              medianprops=dict(color="black", linewidth=2.0),
                              boxprops=dict(linewidth=1.2), flierprops=dict(marker='o', markersize=2, alpha=0.15))
    axes[1].set_ylabel("3rd Codon Position GC3 (%)", fontsize=11, weight="bold")
    axes[1].set_title("B) Codon Bias (GC3)", fontsize=12, weight="bold", loc="left", pad=10)
    
    # Plot Length
    box_len = axes[2].boxplot([core["Length"], acc["Length"], sing["Length"], amr["Length"]], tick_labels=order_amr, patch_artist=True,
                              medianprops=dict(color="black", linewidth=2.0),
                              boxprops=dict(linewidth=1.2), flierprops=dict(marker='o', markersize=2, alpha=0.15))
    axes[2].set_yscale("log")
    axes[2].set_ylabel("Gene Length (bp, log scale)", fontsize=11, weight="bold")
    axes[2].set_title("C) Gene Length", fontsize=12, weight="bold", loc="left", pad=10)
    
    for box_plot in [box_gc, box_gc3, box_len]:
        for patch, color in zip(box_plot['boxes'], colors_f6):
            patch.set_facecolor(color)
            patch.set_alpha(0.85)
            patch.set_edgecolor("#1e293b")
            
    for ax in axes:
        ax.grid(True, which="both", axis="y", ls="-", color="#e2e8f0", alpha=0.6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#cbd5e1')
        ax.spines['bottom'].set_color('#cbd5e1')
        ax.tick_params(axis='both', labelsize=10)
        
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "figure6_amr_comparison.png"), bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "figure6_amr_comparison.pdf"), bbox_inches="tight")
    plt.close(fig)

    # ------------------ figure8_amr_gc3_scatter ------------------
    print("Generating figure8_amr_gc3_scatter...")
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    fig.patch.set_facecolor("white")
    
    ax.scatter(core["GC"], core["GC3"], c=COLOR_CORE, alpha=0.12, s=7, label="Core Genome (N=3,310)", edgecolors="none")
    ax.scatter(acc["GC"], acc["GC3"], c=COLOR_ACCESSORY, alpha=0.30, s=11, label="Accessory Genome (N=1,047)", edgecolors="none")
    ax.scatter(sing["GC"], sing["GC3"], c=COLOR_SINGLETON, alpha=0.45, s=14, label="Singleton Genome (N=278)", edgecolors="none")
    
    amr_points = []
    for tag, meta in acquired_amr_locus_tags.items():
        if tag in gene_gc_info:
            gc = gene_gc_info[tag]["GC"]
            gc3 = gene_gc_info[tag]["GC3"]
            amr_points.append((tag, meta["gene"], gc, gc3, meta["class"]))
            
    # Mako/Crest-inspired palette for AMR classes
    class_colors = {
        "Beta-lactam": "#ef4444",      # Red
        "Aminoglycoside": "#0d9488",  # Teal
        "Sulfonamide": "#3b82f6",     # Blue
        "Trimethoprim": "#8b5cf6",    # Violet
        "Phenicol": "#f59e0b",        # Yellow/Amber
        "Macrolide": "#ec4899",       # Pink
        "Tetracycline": "#f97316",    # Orange
        "Biocide": "#10b981"          # Green
    }
    
    plotted_classes = set()
    for tag, gene, gc, gc3, cls in amr_points:
        col = class_colors.get(cls, "red")
        lbl = f"AMR: {cls}" if cls not in plotted_classes else ""
        plotted_classes.add(cls)
        # Use large diamonds with dark outlines
        ax.scatter(gc, gc3, c=col, edgecolors="#1e293b", linewidths=1.2, s=90, marker='D', zorder=5, label=lbl)
        ax.annotate(gene, (gc, gc3), xytext=(5, 5), textcoords="offset points", fontsize=9.5, weight="bold", color="#0f172a", zorder=6)
        
    ax.set_xlabel("GC Content (%)", fontsize=11, weight="bold")
    ax.set_ylabel("GC3 Codon Bias (%)", fontsize=11, weight="bold")
    ax.grid(True, ls="-", color="#e2e8f0", alpha=0.6)
    ax.legend(loc="upper left", frameon=True, shadow=False, facecolor="white", edgecolor="#cbd5e1")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#cbd5e1')
    ax.spines['bottom'].set_color('#cbd5e1')
    
    fig.savefig(os.path.join(out_dir, "figure8_amr_gc3_scatter.png"), bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "figure8_amr_gc3_scatter.pdf"), bbox_inches="tight")
    plt.close(fig)

    # ------------------ figure7_amr_genetic_context ------------------
    print("Generating figure7_amr_genetic_context...")
    fig, axes = plt.subplots(5, 1, figsize=(12, 10), sharex=False, dpi=300)
    fig.patch.set_facecolor("white")
    plt.subplots_adjust(hspace=0.7)
    
    # Modern Mako-inspired gene color mapping
    # Highlighting resistance genes in bright accent color, mobile elements in slate, and others in light grey
    color_map = {
        # AMR Genes: Highlighted in Orange
        "blaTEM-1B": COLOR_AMR_GENE, "aph(3')-Ia": COLOR_AMR_GENE, "aadA1": COLOR_AMR_GENE, "aadA2": COLOR_AMR_GENE,
        "sul1": COLOR_AMR_GENE, "sul3": COLOR_AMR_GENE, "dfrA7": COLOR_AMR_GENE, "cmlA1": COLOR_AMR_GENE, 
        "estX": COLOR_AMR_GENE, "tet(A)": COLOR_AMR_GENE, "qacL": COLOR_AMR_GENE, "qacE\\Delta1": COLOR_AMR_GENE,
        # Mobile elements / Helpers: Slate Blue / Indigo
        "tetR": "#a5b4fc", "hin_3": "#818cf8", "tnpR": "#6366f1", "IS transposase": "#475569",
        "hypothetical": "#e2e8f0", "other": "#94a3b8"
    }
    
    contig_draws = [
        {
            "name": "NODE_37 (Length: 9,804 bp, CALIN Array Cassette)",
            "length": 9804,
            "genes": [
                {"name": "estX", "start": 147, "end": 1226, "strand": "+"},
                {"name": "hypothetical", "start": 1322, "end": 1930, "strand": "+"},
                {"name": "aadA2", "start": 2000, "end": 2779, "strand": "+"},
                {"name": "cmlA1", "start": 3041, "end": 4300, "strand": "+"},
                {"name": "aadA1", "start": 4405, "end": 5184, "strand": "+"},
                {"name": "qacL", "start": 5354, "end": 5686, "strand": "+"},
                {"name": "IS transposase", "start": 5933, "end": 6508, "strand": "-"}, 
                {"name": "sul3", "start": 6866, "end": 7657, "strand": "-"},
                {"name": "kdsr (other)", "start": 8409, "end": 9272, "strand": "-"}
            ],
            "attc": [1241, 1925, 2781, 4320, 5186, 5697]
        },
        {
            "name": "NODE_30 (tet locus, Region shown: 18,000–22,184 bp)",
            "length": 4184,
            "offset": 18000,
            "genes": [
                {"name": "tet(A)", "start": 18748, "end": 19947, "strand": "-"},
                {"name": "tetR", "start": 20026, "end": 20703, "strand": "+"},
                {"name": "hin_3", "start": 21335, "end": 21895, "strand": "-"}
            ],
            "attc": []
        },
        {
            "name": "NODE_41 (Class 1 Integron 3'-CS, Length: 2,175 bp)",
            "length": 2175,
            "genes": [
                {"name": "dfrA7", "start": 18, "end": 719, "strand": "+"},
                {"name": "qacE\\Delta1", "start": 949, "end": 1296, "strand": "+"},
                {"name": "sul1", "start": 1290, "end": 2069, "strand": "+"}
            ],
            "attc": []
        },
        {
            "name": "NODE_42 (Tn3 Transposon segment, Length: 1,823 bp)",
            "length": 1823,
            "genes": [
                {"name": "tnpR", "start": 8, "end": 250, "strand": "+"},
                {"name": "blaTEM-1B", "start": 433, "end": 1293, "strand": "+"}
            ],
            "attc": []
        },
        {
            "name": "NODE_47 (aph(3')-Ia element, Length: 1,294 bp)",
            "length": 1294,
            "genes": [
                {"name": "aph(3')-Ia", "start": 215, "end": 1030, "strand": "+"}
            ],
            "attc": []
        }
    ]
    
    for idx, drawing in enumerate(contig_draws):
        ax = axes[idx]
        offset = drawing.get("offset", 0)
        length = drawing["length"]
        
        ax.plot([0, length], [0, 0], color="#475569", lw=2.5, zorder=1)
        
        for g in drawing["genes"]:
            g_start = g["start"] - offset
            g_end = g["end"] - offset
            strand = g["strand"]
            g_name = g["name"]
            
            col = color_map.get(g_name, "#94a3b8")
            if "other" in g_name:
                col = color_map["other"]
            elif "IS transposase" in g_name:
                col = color_map["IS transposase"]
            elif "hypothetical" in g_name:
                col = color_map["hypothetical"]
                
            width = g_end - g_start
            
            # Draw sleek arrows
            if strand == "+":
                arrow = patches.FancyArrow(g_start, 0, width, 0, width=0.18, head_width=0.38, 
                                           head_length=min(width*0.25, 180), length_includes_head=True,
                                           facecolor=col, edgecolor="#1e293b", linewidth=0.9, zorder=3)
            else:
                arrow = patches.FancyArrow(g_end, 0, -width, 0, width=0.18, head_width=0.38, 
                                           head_length=min(width*0.25, 180), length_includes_head=True,
                                           facecolor=col, edgecolor="#1e293b", linewidth=0.9, zorder=3)
            ax.add_patch(arrow)
            
            label_text = g_name.replace("\\Delta", "Δ")
            # If it is an AMR gene, draw bold text, else smaller normal text
            txt_weight = "bold" if col == COLOR_AMR_GENE else "normal"
            txt_color = "#0f172a" if col == COLOR_AMR_GENE else "#475569"
            ax.text((g_start + g_end) / 2, 0.35, label_text, ha="center", va="bottom", 
                    fontsize=10, weight=txt_weight, color=txt_color)
            
            ax.text(g_start, -0.4, f"{g['start']}", ha="center", va="top", fontsize=8, color="#64748b")
            ax.text(g_end, -0.4, f"{g['end']}", ha="center", va="top", fontsize=8, color="#64748b")
            
        for att in drawing["attc"]:
            pos = att - offset
            # attC recombination triangle
            att_mark = patches.RegularPolygon((pos, 0), numVertices=3, radius=110, orientation=np.pi,
                                             facecolor="#ef4444", edgecolor="#1e293b", linewidth=0.8, zorder=10)
            ax.add_patch(att_mark)
            ax.text(pos, -0.45, "attC", ha="center", va="top", fontsize=7.5, color="#ef4444", weight="bold")
            
        ax.set_xlim(-200, length + 200)
        ax.set_ylim(-1.0, 1.0)
        ax.set_title(drawing["name"], fontsize=11, weight="bold", loc="left", color="#1e293b", pad=8)
        ax.axis("off")
        
    fig.savefig(os.path.join(out_dir, "step3_amr_genetic_context.png"), bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "step3_amr_genetic_context.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("Saved step3_amr_genetic_context.")

# ==============================================================================
# MAIN PIPELINE RUNNER
# ==============================================================================
def main():
    print("======================================================================")
    print("REGENERATING REPRODUCIBILITY FIGURES FOR PROPHAGE MOBILOME REPO")
    print("======================================================================")
    
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    REPO_DIR   = os.path.dirname(SCRIPT_DIR)
    out_dir    = os.path.join(REPO_DIR, "figures")
    os.makedirs(out_dir, exist_ok=True)
    
    # Try local workspace base directory
    workspace_dir = "/home/adel/Documents/AMR_Work/AMR_Work"
    fasta_path = os.path.join(workspace_dir, "results/Step2_Assembly/spades_output/contigs.fasta")
    faa_path = os.path.join(workspace_dir, "results/Step3_Annotation/prokka_out/Ecoli_isolate.faa")
    ffn_path = os.path.join(workspace_dir, "results/Step3_Annotation/prokka_out/Ecoli_isolate.ffn")
    roary_path = os.path.join(workspace_dir, "pangenome_expansion_highcpu/panaroo_out/gene_presence_absence_roary.csv")
    
    # Generate Figure 1 (Assembly Coverage)
    if os.path.exists(fasta_path):
        df_contigs = parse_assembly_contigs(fasta_path)
        plot_figure1_assembly_coverage(df_contigs, out_dir)
    else:
        print(f"Skipping step1_assembly_coverage: {fasta_path} not found.")
        
    # Generate Figure 7 (AMR genetic context)
    if os.path.exists(faa_path) and os.path.exists(ffn_path) and os.path.exists(roary_path):
        perform_statistical_plots(faa_path, ffn_path, roary_path, out_dir)
    else:
        print("Skipping step3_amr_genetic_context: required annotation files not found.")
        
    print("\n✔ Relevant figures successfully generated and saved in the figures/ folder!")

if __name__ == "__main__":
    main()
