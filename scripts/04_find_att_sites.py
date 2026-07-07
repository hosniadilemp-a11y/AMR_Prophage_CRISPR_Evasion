#!/usr/bin/env python3
import os
from Bio import SeqIO

def find_direct_repeats_wider():
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

    print(f"NODE_1 length: {len(node1_seq)} bp")

    # Wider search windows
    # Left window: around the transition from yebY (56,393 bp) to first hypothetical gene (56,786 bp)
    left_start = 56000
    left_end = 57000
    left_win = node1_seq[left_start:left_end]

    # Right window: around pphA/yebV/rsmF boundary (98,500 to 99,500 bp)
    right_start = 98500
    right_end = 99500
    right_win = node1_seq[right_start:right_end]

    print(f"Left window: {left_start} to {left_end} bp")
    print(f"Right window: {right_start} to {right_end} bp")

    # Search for matching substrings of length K (10 to 30 bp)
    min_len = 10
    max_len = 30
    matches = []

    for k in range(max_len, min_len - 1, -1):
        for i in range(len(left_win) - k + 1):
            sub = left_win[i:i+k]
            pos_in_right_win = right_win.find(sub)
            if pos_in_right_win != -1:
                abs_left_pos = left_start + i
                abs_right_pos = right_start + pos_in_right_win
                matches.append({
                    "seq": sub,
                    "len": k,
                    "left_pos": abs_left_pos,
                    "right_pos": abs_right_pos
                })
        if matches:
            break

    print(f"\nFound {len(matches)} direct repeat candidates of length {k}:")
    for m in matches:
        print(f"Sequence  : {m['seq']}")
        print(f"Length    : {m['len']} bp")
        print(f"attL (left) coordinate  : {m['left_pos']} to {m['left_pos'] + m['len']} bp")
        print(f"attR (right) coordinate : {m['right_pos']} to {m['right_pos'] + m['len']} bp")
        print(f"Distance between repeats : {m['right_pos'] - m['left_pos']} bp")
        print("-" * 50)

if __name__ == "__main__":
    find_direct_repeats_wider()
