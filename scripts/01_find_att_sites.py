#!/usr/bin/env python3
import os
import yaml
from Bio import SeqIO

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def find_flexible_repeats():
    # Load configuration if available
    fasta_path = "data/QA5221.fasta"
    config_path = "config/pipeline_config.yaml"
    
    if os.path.exists(config_path):
        try:
            with open(config_path, "r") as f:
                config = yaml.safe_load(f)
                fasta_path = config.get("assembly_fasta", fasta_path)
        except Exception as e:
            print(f"Warning: Could not parse config, using default assembly path. ({e})")
            
    if not os.path.exists(fasta_path):
        # Fall back to parent directory workspace path if running in parent Cwd
        workspace_fallback = "../pangenome_expansion/QA5221.fasta"
        if os.path.exists(workspace_fallback):
            fasta_path = workspace_fallback
        else:
            print(f"Error: {fasta_path} not found. Please place QA5221.fasta in the data/ directory.")
            return

    node1_seq = None
    for record in SeqIO.parse(fasta_path, "fasta"):
        if "NODE_1_" in record.id or record.id == "NODE_1":
            node1_seq = str(record.seq)
            break

    if not node1_seq:
        print("Error: NODE_1 not found in assembly.")
        return

    # Windows around boundaries
    left_start = 56000
    left_end = 57000
    left_win = node1_seq[left_start:left_end]

    right_start = 98000
    right_end = 99300
    right_win = node1_seq[right_start:right_end]

    print(f"Searching for repeats of length K between Left ({left_start}-{left_end}) and Right ({right_start}-{right_end}) windows...")

    min_len = 8
    max_len = 15
    matches = []

    # Search with up to 1 mismatch
    for k in range(max_len, min_len - 1, -1):
        for i in range(len(left_win) - k + 1):
            sub_left = left_win[i:i+k]
            for j in range(len(right_win) - k + 1):
                sub_right = right_win[j:j+k]
                
                dist = hamming_distance(sub_left, sub_right)
                if dist <= 1:
                    abs_left = left_start + i
                    abs_right = right_start + j
                    matches.append({
                        "left_seq": sub_left,
                        "right_seq": sub_right,
                        "mismatches": dist,
                        "len": k,
                        "left_pos": abs_left,
                        "right_pos": abs_right
                    })
        if matches:
            # Sort by mismatches first, then length
            matches.sort(key=lambda x: (x["mismatches"], -x["len"]))
            break

    os.makedirs("results", exist_ok=True)
    out_file = "results/att_sites_output.txt"
    
    with open(out_file, "w") as f:
        f.write(f"Found {len(matches)} direct repeat candidates of length {k} (with <= 1 mismatch):\n\n")
        for m in matches[:15]:
            f.write(f"Left Seq  : {m['left_seq']}\n")
            f.write(f"Right Seq : {m['right_seq']}\n")
            f.write(f"Mismatches: {m['mismatches']}\n")
            f.write(f"Length    : {m['len']} bp\n")
            f.write(f"attL (left) coordinate  : {m['left_pos']} to {m['left_pos'] + m['len']} bp\n")
            f.write(f"attR (right) coordinate : {m['right_pos']} to {m['right_pos'] + m['len']} bp\n")
            f.write(f"Distance between repeats : {m['right_pos'] - m['left_pos']} bp\n")
            f.write("-" * 50 + "\n")
            
    print(f"Results written to {out_file}")

if __name__ == "__main__":
    find_flexible_repeats()
