#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
import pandas as pd
import mdtraj as md

def main():
    parser = argparse.ArgumentParser(description="Analyze Membrane MD Trajectories using MDTraj (PBC-aware)")
    parser.add_argument("--traj", required=True, help="Path to trajectory DCD file")
    parser.add_argument("--top", required=True, help="Path to topology solvated PDB file")
    parser.add_argument("--outdir", required=True, help="Directory to save output CSVs")
    args = parser.parse_args()

    if not os.path.exists(args.traj):
        print(f"[ERROR] Trajectory file '{args.traj}' not found.", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.top):
        print(f"[ERROR] Topology file '{args.top}' not found.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading topology and selecting protein atoms...")
    top = md.load_topology(args.top)
    
    # Select protein atoms and backbone
    protein_indices = top.select("protein")
    backbone_indices = top.select("protein and backbone")
    ca_indices = top.select("protein and name CA")
    
    n_residues = len(ca_indices)
    print(f"Detected {len(protein_indices)} protein atoms, {len(backbone_indices)} backbone atoms, and {n_residues} residues.")

    print("Loading full trajectory (including solvent/membrane for correct PBC wrapping)...")
    # Load full trajectory
    traj = md.load(args.traj, top=args.top)
    print(f"Loaded trajectory with {traj.n_frames} frames.")

    print("Wrapping/imaging molecules across periodic boundaries...")
    # Wrap coordinates to keep the protein contiguous and centered
    traj.image_molecules(inplace=True)

    print("Slicing out protein atoms for analysis...")
    protein_indices = traj.topology.select("protein")
    traj = traj.atom_slice(protein_indices)

    # Extract backbone and CA indices relative to the sliced trajectory
    backbone_indices = traj.topology.select("backbone")
    ca_indices = traj.topology.select("name CA")

    print("Superposing trajectory on backbone (first frame as reference)...")
    traj.superpose(traj, 0, atom_indices=backbone_indices)

    print("Calculating backbone RMSD...")
    rmsd = md.rmsd(traj, traj, 0, atom_indices=backbone_indices) * 10.0 # nm to Angstrom

    print("Calculating Radius of Gyration...")
    rg = md.compute_rg(traj) * 10.0 # nm to Angstrom

    print("Calculating C-alpha RMSF...")
    rmsf = md.rmsf(traj, traj, 0, atom_indices=ca_indices) * 10.0 # nm to Angstrom

    # Save to CSV
    timeseries_csv = os.path.join(args.outdir, "md_analysis_timeseries_wrapped.csv")
    rmsf_csv = os.path.join(args.outdir, "md_analysis_rmsf_wrapped.csv")

    df_timeseries = pd.DataFrame({
        'Frame': np.arange(len(rmsd)),
        'Time_ns': np.arange(len(rmsd)) * 0.02,  # 20 ps interval = 0.02 ns
        'RMSD_Angstrom': rmsd,
        'Rg_Angstrom': rg
    })
    df_timeseries.to_csv(timeseries_csv, index=False)

    df_rmsf = pd.DataFrame({
        'Residue': np.arange(1, len(rmsf) + 1),
        'RMSF_Angstrom': rmsf
    })
    df_rmsf.to_csv(rmsf_csv, index=False)

    print("\n" + "="*40)
    print("   PBC-WRAPPED MD SIMULATION SUMMARY STATISTICS")
    print("="*40)
    print(f"Total Simulation Time : {df_timeseries['Time_ns'].iloc[-1]:.2f} ns ({traj.n_frames} frames)")
    print(f"Backbone RMSD (overall): {rmsd.mean():.4f} Å ± {rmsd.std():.4f} Å (Min: {rmsd.min():.2f} Å, Max: {rmsd.max():.2f} Å)")
    print(f"Backbone RMSD (last 20ns): {rmsd[int(len(rmsd)*0.6):].mean():.4f} Å ± {rmsd[int(len(rmsd)*0.6):].std():.4f} Å")
    print(f"Radius of Gyration (Rg): {rg.mean():.4f} Å ± {rg.std():.4f} Å")
    print(f"C-alpha RMSF (Average) : {rmsf.mean():.4f} Å (Max: {rmsf.max():.2f} Å at Residue {rmsf.argmax()+1})")
    print("="*40 + "\n")

if __name__ == "__main__":
    main()
