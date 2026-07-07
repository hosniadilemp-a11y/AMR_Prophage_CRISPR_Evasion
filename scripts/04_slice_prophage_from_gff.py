#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def main():
    fasta_path = "results/Step2_Assembly/spades_output/contigs.fasta"
    gff_path = "results/Step3_Annotation/prokka_out/Ecoli_isolate.gff"
    out_path = "manuscript/paper2/plots/prophage_QA5221.gbk"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    if not os.path.exists(fasta_path) or not os.path.exists(gff_path):
        print("Error: FASTA or GFF file is missing.")
        sys.exit(1)
        
    print("Loading NODE_1 sequence...")
    node1_seq = None
    for record in SeqIO.parse(fasta_path, "fasta"):
        if "NODE_1_" in record.id:
            node1_seq = str(record.seq)
            break
            
    if not node1_seq:
        print("Error: NODE_1 sequence not found.")
        sys.exit(1)
        
    # Coordinates of prophage
    attL = 56422
    attR = 98747
    prophage_seq = node1_seq[attL - 1 : attR]
    
    # Create SeqRecord
    record = SeqRecord(
        Seq(prophage_seq),
        id="QA5221_prophage",
        name="QA5221_phage",
        description="Escherichia coli QA5221 active prophage region"
    )
    # Add molecule_type annotation required by modern Biopython
    record.annotations["molecule_type"] = "DNA"
    
    print("Parsing features from GFF...")
    features = []
    
    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom = parts[0]
            feature_type = parts[2]
            if "NODE_1_" not in chrom:
                continue
            if feature_type not in ["CDS", "tRNA"]:
                continue
                
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            
            # Check if feature is within the prophage region
            if start >= attL and end <= attR:
                attrs = parts[8]
                qualifiers = {}
                for attr in attrs.split(";"):
                    if "=" in attr:
                        k, v = attr.split("=", 1)
                        # GenBank format qualifiers must be lists
                        qualifiers[k] = [v]
                
                # Shift coordinates
                # GFF is 1-indexed, inclusive. 
                # Biopython FeatureLocation is 0-indexed, half-open [start, end)
                # Sliced sequence starts at attL - 1
                shifted_start = start - attL
                shifted_end = end - attL + 1
                
                biopython_strand = 1 if strand == "+" else -1
                
                loc = FeatureLocation(shifted_start, shifted_end, strand=biopython_strand)
                feat = SeqFeature(
                    loc,
                    type=feature_type,
                    qualifiers=qualifiers
                )
                features.append(feat)
                
    # Add translation to CDS features if possible
    # We can translate the sequence directly!
    for feat in features:
        if feat.type == "CDS":
            feat_seq = feat.extract(record.seq)
            # Translate sequence
            try:
                # If length is not multiple of 3, truncate
                rem = len(feat_seq) % 3
                if rem != 0:
                    feat_seq = feat_seq[:-rem]
                trans = str(feat_seq.translate(table=11, cds=False))
                # Strip trailing stop codon character '*' if present
                if trans.endswith("*"):
                    trans = trans[:-1]
                feat.qualifiers["translation"] = [trans]
            except Exception as e:
                print(f"Warning: could not translate feature at {feat.location}: {e}")
                
    record.features = features
    
    # Save GenBank
    with open(out_path, "w") as out_f:
        SeqIO.write(record, out_f, "genbank")
        
    print(f"Successfully created prophage GenBank file with {len(features)} features at '{out_path}'")

if __name__ == "__main__":
    main()
