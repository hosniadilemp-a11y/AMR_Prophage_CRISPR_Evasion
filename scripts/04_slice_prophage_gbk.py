#!/usr/bin/env python3
import os
from Bio import SeqIO

def slice_prophage():
    gbk_path = "results/Step3_Annotation/prokka_out/Ecoli_isolate.gbk"
    out_path = "manuscript/paper2/plots/prophage_QA5221.gbk"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    if not os.path.exists(gbk_path):
        print(f"Error: {gbk_path} not found.")
        return
        
    print("Loading GenBank record...")
    target_record = None
    for record in SeqIO.parse(gbk_path, "genbank"):
        if "NODE_1_" in record.id:
            target_record = record
            break
            
    if not target_record:
        print("Error: NODE_1 record not found.")
        return
        
    print(f"Original record ID: {target_record.id}, length: {len(target_record)} bp")
    
    # Coordinates of prophage
    attL = 56422
    attR = 98747
    
    # Slice the record (0-indexed, half-open interval)
    # We slice from attL-1 to attR
    sliced_record = target_record[attL - 1 : attR]
    sliced_record.id = "QA5221_prophage"
    sliced_record.name = "QA5221_phage"
    sliced_record.description = "Escherichia coli QA5221 active prophage region"
    
    # Adjust feature positions
    # Biopython handles feature position shifting automatically during slicing!
    
    # Write to file
    with open(out_path, "w") as out_f:
        SeqIO.write(sliced_record, out_f, "genbank")
        
    print(f"Successfully sliced prophage record (length: {len(sliced_record)} bp) and saved to {out_path}")

if __name__ == "__main__":
    slice_prophage()
