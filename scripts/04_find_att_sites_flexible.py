#!/usr/bin/env python3
import os
from Bio import SeqIO

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def find_flexible_repeats():
    fasta_path = "results/Step2_Assembly/spades_output/contigs.fasta"
    if not os.path.exists(fasta_path):
        print(f"Error: {fasta_path} not found.")
        return

    node1_seq = None
    for record in SeqIO.parse(fasta_path, "fasta"):
        if "NODE_1_" in record.id:
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

    print(f"\nFound {len(matches)} direct repeat candidates of length {k} (with <= 1 mismatch):")
    # Print top 15 candidates
    for m in matches[:15]:
        print(f"Left Seq  : {m['left_seq']}")
        print(f"Right Seq : {m['right_seq']}")
        print(f"Mismatches: {m['mismatches']}")
        print(f"Length    : {m['len']} bp")
        print(f"attL (left) coordinate  : {m['left_pos']} to {m['left_pos'] + m['len']} bp")
        print(f"attR (right) coordinate : {m['right_pos']} to {m['right_pos'] + m['len']} bp")
        print(f"Distance between repeats : {m['right_pos'] - m['left_pos']} bp")
        print("-" * 50)

if __name__ == "__main__":
    find_flexible_repeats()
