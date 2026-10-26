#!/usr/bin/env python3
"""
Generate Supplementary Table S4: Genome-Wide Prophage Repertoire Audit for E. coli QA5221
Addresses Reviewer 2 (Domain 2.5), Issue I10 (21-row expansion), and Issue I29 (clean tabularx formatting).
"""
import os
import pandas as pd

prophage_rows = [
    {
        "Prophage ID": "PP-01 (Focal)",
        "Contig": "NODE_1",
        "Coordinates (bp)": "56,422--98,747",
        "Length (kb)": "42.3",
        "CDS Count": "62",
        "Integrity Status": "Genomically Intact",
        "Modules Present": "Lysogeny, Recombination, Replication, Packaging, Structural, Lysis",
        "Virulence / AMR Cargo": "Enterohaemolysin (\\textit{AC4NUP_00305}), \\textit{rcbA}",
        "Priority Rationale": "\\textbf{Prioritised for deep characterisation:} Sole complete, fully modular prophage in genome with intact virulence cargo and \\textit{attL/attR} direct repeats."
    },
    {
        "Prophage ID": "PP-02",
        "Contig": "NODE_1",
        "Coordinates (bp)": "375,146--400,509",
        "Length (kb)": "25.4",
        "CDS Count": "35",
        "Integrity Status": "Incomplete / Degraded",
        "Modules Present": "Partial replication, putative recombination",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Degenerate remnant; lacks structural capsid, packaging terminase, and lysis cassettes."
    },
    {
        "Prophage ID": "PP-03",
        "Contig": "NODE_28",
        "Coordinates (bp)": "13,809--26,361",
        "Length (kb)": "12.5",
        "CDS Count": "26",
        "Integrity Status": "Incomplete / Cryptic",
        "Modules Present": "Integrase, partial structural fragments",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Cryptic prophage fragment; lacking tail assembly, baseplate, and lysis modules."
    },
    {
        "Prophage ID": "PP-04",
        "Contig": "NODE_7",
        "Coordinates (bp)": "42,638--44,585",
        "Length (kb)": "1.9",
        "CDS Count": "2",
        "Integrity Status": "Defective Remnant",
        "Modules Present": "Integrase remnant",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Decayed horizontal insertion; insufficient coding capacity."
    },
    {
        "Prophage ID": "PP-05",
        "Contig": "NODE_9",
        "Coordinates (bp)": "187,885--189,069",
        "Length (kb)": "1.2",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Integrase fragment",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Decayed horizontal scar; solitary integrase lacking replication and structural cassettes."
    },
    {
        "Prophage ID": "PP-06",
        "Contig": "NODE_9",
        "Coordinates (bp)": "155,788--156,777",
        "Length (kb)": "1.0",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Transposase fragment",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Insertion element remnant; no modular viral integrity."
    },
    {
        "Prophage ID": "PP-07",
        "Contig": "NODE_8",
        "Coordinates (bp)": "159,454--161,712",
        "Length (kb)": "2.3",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Phage-related protein",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Degraded prophage scar; lacks capsid, packaging, and regulatory machinery."
    },
    {
        "Prophage ID": "PP-08",
        "Contig": "NODE_5",
        "Coordinates (bp)": "220,328--221,539",
        "Length (kb)": "1.2",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Integrase remnant",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Defective horizontal scar; single ORF lacking replication or lysogeny modules."
    },
    {
        "Prophage ID": "PP-09",
        "Contig": "NODE_4",
        "Coordinates (bp)": "316,174--319,146",
        "Length (kb)": "3.0",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Hypothetical protein",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Decayed remnant; solitary ORF with no phage structural genes."
    },
    {
        "Prophage ID": "PP-10",
        "Contig": "NODE_3",
        "Coordinates (bp)": "305,473--306,150",
        "Length (kb)": "0.7",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Phage tail fragment",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Cryptic structural remnant without assembly genes or virion structure."
    },
    {
        "Prophage ID": "PP-11",
        "Contig": "NODE_3",
        "Coordinates (bp)": "292,532--293,791",
        "Length (kb)": "1.3",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Integrase fragment",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Decayed insertion scar; lacks replication, capsid, and tail modules."
    },
    {
        "Prophage ID": "PP-12",
        "Contig": "NODE_23",
        "Coordinates (bp)": "40,362--41,357",
        "Length (kb)": "1.0",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Transposase remnant",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Mobile element scar; solitary ORF without viral cassettes."
    },
    {
        "Prophage ID": "PP-13",
        "Contig": "NODE_23",
        "Coordinates (bp)": "14,741--15,931",
        "Length (kb)": "1.2",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Integrase remnant",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Defective horizontal insertion; single isolated ORF."
    },
    {
        "Prophage ID": "PP-14",
        "Contig": "NODE_22",
        "Coordinates (bp)": "7,020--8,204",
        "Length (kb)": "1.2",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Integrase fragment",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Degraded integrase scar; lacks structural and lysis cassettes."
    },
    {
        "Prophage ID": "PP-15",
        "Contig": "NODE_2",
        "Coordinates (bp)": "123--1,067",
        "Length (kb)": "0.9",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Hypothetical viral ORF",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Contig edge artifact; solitary fragment without viral machinery."
    },
    {
        "Prophage ID": "PP-16",
        "Contig": "NODE_18",
        "Coordinates (bp)": "49,452--50,348",
        "Length (kb)": "0.9",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Phage protein remnant",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Degenerate viral scar; lacks packaging and capsid machinery."
    },
    {
        "Prophage ID": "PP-17",
        "Contig": "NODE_17",
        "Coordinates (bp)": "90,468--91,595",
        "Length (kb)": "1.1",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Integrase fragment",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Defective insertion scar; isolated ORF without modular viral integrity."
    },
    {
        "Prophage ID": "PP-18",
        "Contig": "NODE_17",
        "Coordinates (bp)": "81,513--82,712",
        "Length (kb)": "1.2",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Hypothetical protein",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Cryptic remnant lacking functional modules or virulence cargo."
    },
    {
        "Prophage ID": "PP-19",
        "Contig": "NODE_13",
        "Coordinates (bp)": "82,401--83,297",
        "Length (kb)": "0.9",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Transposase fragment",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Mobile element remnant; solitary transposase scar."
    },
    {
        "Prophage ID": "PP-20",
        "Contig": "NODE_13",
        "Coordinates (bp)": "4,848--6,038",
        "Length (kb)": "1.2",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Integrase remnant",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Decayed horizontal insertion scar; single ORF."
    },
    {
        "Prophage ID": "PP-21",
        "Contig": "NODE_1",
        "Coordinates (bp)": "106,172--108,220",
        "Length (kb)": "2.0",
        "CDS Count": "1",
        "Integrity Status": "Solitary Remnant",
        "Modules Present": "Phage-related ORF",
        "Virulence / AMR Cargo": "None",
        "Priority Rationale": "Degraded remnant adjacent to ribosomal locus; lacks replication, structural, and packaging modules."
    }
]

df = pd.DataFrame(prophage_rows)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Save TSV
tsv_path = os.path.join(BASE_DIR, "TableS4_Prophage_Repertoire.tsv")
df.to_csv(tsv_path, sep="\t", index=False)

# Save Markdown
md_path = os.path.join(BASE_DIR, "TableS4_Prophage_Repertoire.md")
with open(md_path, "w") as f:
    f.write("# Supplementary Table S4: Repertoire and Genomic Integrity Audit of Predicted Prophage Elements Across E. coli QA5221\n\n")
    f.write("> **Reviewer 2 (Domain 2.5 / Issues I10 & I29) Resolution:** This table documents the comprehensive genome-wide prophage audit across all 37 contigs evaluated by PhiSpy, PHASTEST, and manual curation. Out of 21 candidate regions detected, **only the prophage on NODE_1 (56,422–98,747 bp, 42.3 kb, 62 CDSs) possesses complete genomic and modular integrity alongside confirmed pathogenic virulence cargo (enterohaemolysin).** All other candidate loci (PP-02 through PP-21) are defective pseudoprophages, decayed remnants, or isolated mobile element scars.\n\n")
    clean_df = df.copy()
    for col in clean_df.columns:
        clean_df[col] = clean_df[col].str.replace(r'\\textit\{([^}]+)\}', r'*\1*', regex=True)
        clean_df[col] = clean_df[col].str.replace(r'\\textbf\{([^}]+)\}', r'**\1**', regex=True)
    headers = list(clean_df.columns)
    f.write("| " + " | ".join(headers) + " |\n")
    f.write("| " + " | ".join(["---"] * len(headers)) + " |\n")
    for _, row in clean_df.iterrows():
        f.write("| " + " | ".join(str(val) for val in row) + " |\n")
    f.write("\n")

# Save LaTeX table
tex_path = os.path.join(BASE_DIR, "TableS4_Prophage_Repertoire.tex")
with open(tex_path, "w") as f:
    f.write("% Supplementary Table S4: Genome-wide Prophage Repertoire Audit (Reviewer 2, Domain 2.5, Issues I10 & I29)\n")
    f.write("\\begin{table*}[p]\n\\centering\\footnotesize\n")
    f.write("\\caption{Genome-wide audit of candidate prophage elements and defective pseudoprophage remnants identified in \\textit{E. coli} QA5221 across all 37 contigs ($\\ge 5,000$~bp), providing scientific justification for the focal characterisation of the intact prophage on contig NODE\\_1 (Table S4).}\n")
    f.write("\\label{tab:prophage_repertoire_audit}\n")
    f.write("\\vspace{0.2cm}\n")
    f.write("\\begin{tabularx}{\\linewidth}{lp{1.6cm}rccp{2.6cm}p{3.2cm}p{2.8cm}X}\n")
    f.write("\\toprule\n")
    f.write("\\textbf{ID} & \\textbf{Contig} & \\textbf{Coordinates (bp)} & \\textbf{Size (kb)} & \\textbf{CDS} & \\textbf{Integrity Status} & \\textbf{Modules Present} & \\textbf{Virulence / Cargo} & \\textbf{Prioritisation Rationale} \\\\\n")
    f.write("\\midrule\n")
    for r in prophage_rows:
        contig_escaped = r['Contig'].replace('_', '\\_')
        f.write(f"{r['Prophage ID']} & {contig_escaped} & {r['Coordinates (bp)']} & {r['Length (kb)']} & {r['CDS Count']} & {r['Integrity Status']} & {r['Modules Present']} & {r['Virulence / AMR Cargo']} & {r['Priority Rationale']} \\\\\n")
    f.write("\\bottomrule\n")
    f.write("\\end{tabularx}\n\\end{table*}\n")

print(f"Generated {tsv_path}, {md_path}, and {tex_path} with all 21 prophages successfully.")
