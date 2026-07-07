#!/usr/bin/env python3
import os
import subprocess
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# Set consistent fonts
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
plt.rcParams['text.color'] = '#1e293b'
plt.rcParams['axes.labelcolor'] = '#1e293b'

# Paths
AMR_WORK_DIR = "/media/adel/Data/Hosni/openmm_windows_setup/AMR_Work"
REF_GENOMES_DIR = os.path.join(AMR_WORK_DIR, "results/Step3_Annotation/reference_genomes")
QA5221_FASTA = os.path.join(AMR_WORK_DIR, "pangenome_expansion/QA5221.fasta")
OUTPUT_DIR = os.path.join(AMR_WORK_DIR, "results/Step3_Annotation/amr_abricate")
PLOTS_DIR = os.path.join(AMR_WORK_DIR, "manuscript/paper2/plots")

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

def run_abricate_on_genomes():
    # Find all reference genomes starting with ST354_
    ref_fastas = glob.glob(os.path.join(REF_GENOMES_DIR, "ST354_*.fna"))
    print(f"Found {len(ref_fastas)} ST354 reference fasta files.")
    
    all_targets = [("QA5221", QA5221_FASTA)]
    for f in ref_fastas:
        name = os.path.basename(f).replace(".fna", "")
        all_targets.append((name, f))
        
    print(f"Total genomes to analyze: {len(all_targets)}")
    
    # Run abricate on each genome
    for name, filepath in all_targets:
        out_tsv = os.path.join(OUTPUT_DIR, f"{name}_resfinder.tsv")
        if not os.path.exists(out_tsv):
            print(f"Running ABRicate on {name}...")
            cmd = f"conda run -n amr_env abricate --db resfinder {filepath} > {out_tsv}"
            subprocess.run(cmd, shell=True, check=True)
        else:
            print(f"ABRicate report already exists for {name}")
            
    return all_targets

def parse_and_network(genomes):
    # Parse results and aggregate
    gene_profiles = {}
    
    for name, _ in genomes:
        out_tsv = os.path.join(OUTPUT_DIR, f"{name}_resfinder.tsv")
        if not os.path.exists(out_tsv):
            continue
            
        df = pd.read_csv(out_tsv, sep='\t')
        # Filter for high quality hits: identity >= 80%, coverage >= 80%
        df_filtered = df[(df['%IDENTITY'] >= 80.0) & (df['%COVERAGE'] >= 80.0)]
        
        genes = df_filtered['GENE'].unique()
        gene_profiles[name] = {gene: 1 for gene in genes}
        
    # Create DataFrame
    matrix_df = pd.DataFrame(gene_profiles).fillna(0).astype(int)
    print(f"AMR matrix built: {matrix_df.shape[0]} AMR genes x {matrix_df.shape[1]} genomes")
    
    # Save matrix
    matrix_path = os.path.join(OUTPUT_DIR, "amr_matrix.tsv")
    matrix_df.to_csv(matrix_path, sep='\t')
    print(f"AMR matrix saved to {matrix_path}")
    
    # Compute Jaccard correlation between AMR genes
    genes = matrix_df.index.tolist()
    n_genes = len(genes)
    
    jaccard_matrix = np.zeros((n_genes, n_genes))
    for i in range(n_genes):
        for j in range(n_genes):
            g_i = matrix_df.iloc[i].values
            g_j = matrix_df.iloc[j].values
            
            intersection = np.sum((g_i == 1) & (g_j == 1))
            union = np.sum((g_i == 1) | (g_j == 1))
            
            jaccard_matrix[i, j] = (intersection / union) if union > 0 else 0.0
            
    # Build NetworkX Graph
    G = nx.Graph()
    for g in genes:
        # Node weight based on general frequency in the 28 genomes
        freq = float(matrix_df.loc[g].mean())
        G.add_node(g, frequency=freq)
        
    # Add edges for high co-occurrence (Jaccard >= 0.70)
    for i in range(n_genes):
        for j in range(i+1, n_genes):
            weight = jaccard_matrix[i, j]
            if weight >= 0.70:
                G.add_edge(genes[i], genes[j], weight=weight)
                
    # Plot network
    plt.figure(figsize=(10, 8), dpi=300)
    
    # Position nodes using spring layout with larger spacing (k=0.8) and more iterations
    pos = nx.spring_layout(G, k=0.8, seed=42, iterations=100)
    
    # Node colors based on frequency
    node_freqs = [nx.get_node_attributes(G, 'frequency')[node] for node in G.nodes()]
    node_sizes = [400 + freq * 1600 for freq in node_freqs]
    
    # Custom LinearSegmentedColormap (Soft Blue to Deep Navy)
    from matplotlib.colors import LinearSegmentedColormap
    cmap_node = LinearSegmentedColormap.from_list('amr_theme', ['#93c5fd', '#1e3a8a'])
    
    # Draw nodes with dark outlines
    nx.draw_networkx_nodes(
        G, pos, 
        node_size=node_sizes, 
        node_color=node_freqs, 
        cmap=cmap_node, 
        edgecolors='#0f172a', 
        linewidths=1.5
    )
    
    # Draw edges with varying thickness and transparency
    edges = G.edges(data=True)
    weights = [data['weight'] * 5 for u, v, data in edges]
    nx.draw_networkx_edges(G, pos, width=weights, edge_color='#94a3b8', alpha=0.5)
    
    # Draw labels with white outline for perfect legibility
    import matplotlib.patheffects as path_effects
    labels_dict = nx.draw_networkx_labels(
        G, pos, 
        font_size=8, 
        font_weight='bold', 
        font_color='#0f172a',
        font_family='DejaVu Sans'
    )
    for text in labels_dict.values():
        text.set_path_effects([
            path_effects.Stroke(linewidth=2.5, foreground='white'),
            path_effects.Normal()
        ])
    
    plt.axis('off')
    plt.tight_layout()
    
    # Paths for paper2
    pdf_path_p2 = os.path.join(PLOTS_DIR, "figure_amr_network.pdf")
    png_path_p2 = os.path.join(PLOTS_DIR, "figure_amr_network.png")
    
    # Paths for paper1
    PLOTS_DIR_P1 = os.path.join(AMR_WORK_DIR, "manuscript/paper1/plots/final_figures/manuscript/phase4_mge_synteny")
    os.makedirs(PLOTS_DIR_P1, exist_ok=True)
    pdf_path_p1 = os.path.join(PLOTS_DIR_P1, "Fig06-b_AMR_Network.pdf")
    png_path_p1 = os.path.join(PLOTS_DIR_P1, "Fig06-b_AMR_Network.png")
    
    # Save files
    plt.savefig(pdf_path_p2, bbox_inches='tight')
    plt.savefig(png_path_p2, bbox_inches='tight')
    plt.savefig(pdf_path_p1, bbox_inches='tight')
    plt.savefig(png_path_p1, bbox_inches='tight')
    plt.close()
    
    print(f"Network figures saved successfully:")
    print(f"  Paper2: {pdf_path_p2} & {png_path_p2}")
    print(f"  Paper1: {pdf_path_p1} & {png_path_p1}")

if __name__ == "__main__":
    genomes = run_abricate_on_genomes()
    parse_and_network(genomes)
