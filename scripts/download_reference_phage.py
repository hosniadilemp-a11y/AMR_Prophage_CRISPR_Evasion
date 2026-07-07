#!/usr/bin/env python3
import os
import sys
from Bio import Entrez

def download_ref():
    # Set email for NCBI requests
    Entrez.email = "researcher@ens-kouba.dz"
    accession = "NC_001799.1" # E. coli bacteriophage DLP12
    out_path = "data/dlp12_ref.gbk"
    
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    print(f"Downloading reference phage {accession} from NCBI Entrez...")
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gbwithparts", retmode="text")
        data = handle.read()
        handle.close()
        
        with open(out_path, "w") as f:
            f.write(data)
            
        print(f"Successfully downloaded and saved reference to {out_path}")
    except Exception as e:
        print(f"Error downloading reference: {e}")
        sys.exit(1)

if __name__ == "__main__":
    download_ref()
