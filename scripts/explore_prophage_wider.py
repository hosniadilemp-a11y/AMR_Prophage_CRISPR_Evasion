#!/usr/bin/env python3
import os

gff_path = "results/Step3_Annotation/prokka_out/Ecoli_isolate.gff"

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
        if feature_type != "CDS":
            continue
        start = int(parts[3])
        end = int(parts[4])
        if start >= 55000 and end <= 105000:
            attrs = parts[8]
            locus = ""
            product = ""
            gene_name = ""
            for attr in attrs.split(";"):
                if attr.startswith("locus_tag="):
                    locus = attr.split("=")[1]
                elif attr.startswith("product="):
                    product = attr.split("=")[1]
                elif attr.startswith("Name="):
                    gene_name = attr.split("=")[1]
            print(f"{start}\t{end}\t{locus}\t{gene_name}\t{product}")
