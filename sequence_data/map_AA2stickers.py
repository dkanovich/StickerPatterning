#!/usr/bin/env python3
"""
mask_generator.py

Generate amino-acid-level binary masks for protein sequences.
Any residue matching the user-specified criteria is marked '1', others '2'.
Supports input from FASTA, CSV, and Excel files with a 'Sequence' column.
Options are set directly in the script (no CLI parsing).
"""

import os
import pandas as pd

# === USER CONFIG ===
# List your input files here (FASTA, CSV, or Excel)
input_files = [
    "PLCD_Sequences/Cajal_PLCD.fasta",
    "PLCD_Sequences/Cajal_PLCD.fasta",                
    "PLCD_Sequences/Neuronal_Granule_PLCD.fasta",             
    "PLCD_Sequences/PCG_PLCD.fasta",                          
    "PLCD_Sequences/P_Body_PLCD.fasta",                       
    "PLCD_Sequences/Speckle_PLCD.fasta",                      
    "PLCD_Sequences/Stress_Granule_PLCD.fasta",        
    "PLCD_Sequences/Centrosome_PLCD.fasta",           
    "PLCD_Sequences/Nucleolus_PLCD.fasta",                    
    "PLCD_Sequences/PML_PLCD.fasta",                         
    "PLCD_Sequences/Paraspeckle_PLCD.fasta",                 
    "PLCD_Sequences/Spindle_PLCD.fasta"
]
# input_files = [
#     "dataset_July2025/IDR_Lists_MattKing.xlsx"]

# Output text file for masks
output_file = "binary_masks_PLCD.txt"
# Specify individual amino acids to include (e.g. ['A', 'C', 'D'])
amino_acids = ['Y', 'F', 'R']
# Specify named groups to include (choose from: polar, aromatic, hydrophobic, positive, negative, charged, all)
groups = []  # example: ['polar', 'aromatic']

# Amino acid group definitions
aa_groups = {
    'polar': set('STNQY'),
    'aromatic': set('FWYH'),
    'hydrophobic': set('AVILMFWY'),
    'positive': set('KRH'),
    'negative': set('DE'),
    'charged': set('KRHDE'),
    'all': set('ARNDCQEGHILKMFPSTWYV'),
}

# Build target set based on user config
target_set = set()
# Add from groups
for g in groups:
    g_low = g.lower()
    if g_low not in aa_groups:
        raise ValueError(f"Unknown group: {g}")
    target_set |= aa_groups[g_low]
# Add from individual amino acids
for aa in amino_acids:
    aa = aa.upper()
    if len(aa) != 1 or aa not in aa_groups['all']:
        raise ValueError(f"Invalid amino acid code: {aa}")
    target_set.add(aa)

if not target_set:
    raise ValueError("No amino acids specified in 'amino_acids' or 'groups'.")

# FASTA parser
def parse_fasta(path):
    seqs = []
    with open(path, 'r') as f:
        current = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current:
                    seqs.append(''.join(current).upper())
                    current = []
            else:
                current.append(line)
        if current:
            seqs.append(''.join(current).upper())
    return seqs

# Table parser for CSV/Excel
def parse_table(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == '.csv':
        df = pd.read_csv(path)
    elif ext in ('.xls', '.xlsx'):
        df = pd.read_excel(path, engine='openpyxl')
    else:
        raise ValueError(f"Unknown table format: {path}")
    cols = {c.lower(): c for c in df.columns}
    if 'sequence' not in cols:
        raise KeyError(f"No 'Sequence' column in {path}")
    seq_col = cols['sequence']
    return df[seq_col].dropna().astype(str).str.upper().tolist()

# Generate masks and write output
masks = []
for path in input_files:
    if not os.path.isfile(path):
        print(f"Warning: file not found: {path}")
        continue
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.fa', '.fasta'):
        seqs = parse_fasta(path)
    elif ext in ('.csv', '.xls', '.xlsx'):
        seqs = parse_table(path)
    else:
        print(f"Skipping unsupported file type: {path}")
        continue

    for seq in seqs:
        mask = ''.join('1' if aa in target_set else '2' for aa in seq)
        masks.append(mask)

with open(output_file, 'w') as fout:
    for m in masks:
        fout.write(m + '\n')

print(f"Wrote {len(masks)} masks to {output_file}")
