#!/bin/bash

# Inputs a PDB file of a protein structure and extracts the amino acid sequence
# from a specified chain and residue range, storing it in the variable "target_seq".

# Input PDB file (update this path as needed)
target_pdb="$1"
# target_pdb="/home/ubuntu/protein-binder-design/input/target_A60_90.pdb"

# Create a mapping of amino acid abbreviations to their single-letter codes
declare -A amino_acid_map=(
    [ALA]="A" [ARG]="R" [ASN]="N" [ASP]="D" [CYS]="C"
    [GLN]="Q" [GLU]="E" [GLY]="G" [HIS]="H" [ILE]="I"
    [LEU]="L" [LYS]="K" [MET]="M" [PHE]="F" [PRO]="P"
    [SER]="S" [THR]="T" [TRP]="W" [TYR]="Y" [VAL]="V"
)

# Export the amino acid map as a space-separated string to pass to `awk`
aa_map_string=""
for key in "${!amino_acid_map[@]}"; do
    aa_map_string+="$key=${amino_acid_map[$key]} "
done

# Extract the sequence from the PDB file
target_seq=$(awk -v aa_map_string="$aa_map_string" '
    BEGIN {
        # Parse the amino acid map string into an array
        split(aa_map_string, map_array, " ");
        for (i in map_array) {
            split(map_array[i], pair, "=");
            aa_map[pair[1]] = pair[2];
        }

        # Initialize variables
        seq = "";
    }
    $1 == "ATOM" {
        # Extract the residue abbreviation (column 4) and residue number (column 6)
        res_name = $4;
        res_number = $6;

        # If this residue number has not been seen before, add its letter to the sequence
        if (!(res_number in seen_residues)) {
            seen_residues[res_number] = 1;
            aa_letter = (res_name in aa_map) ? aa_map[res_name] : "X";  # Use "X" for unknown residues
            seq = seq aa_letter;
        }
    }
    END {
        # Print the final sequence
        print seq;
    }
' "$target_pdb")

# Print and store the sequence in the variable "target_seq"
echo "${target_seq}"