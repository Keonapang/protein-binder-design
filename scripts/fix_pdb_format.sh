#!/bin/bash

# Input PDB file
input_file=$1
output_file=$2

# Function to check and fix formatting
fix_pdb() {
  awk '
    # Match lines with "ATOM" or "HETATM" (PDB records) and residue numbers >= 1000
    /^ATOM|^HETATM/ {
      chain_and_residue = substr($0, 22, 5) # Extract chain and residue columns (22-26)
      chain = substr(chain_and_residue, 1, 1) # Extract chain (first character)
      residue = substr(chain_and_residue, 2)  # Extract residue number (remaining characters)

      # Check if residue number starts with a digit and fix the formatting
      if (residue ~ /^[0-9]+$/) {
        # Correct formatting by inserting a space between chain and residue
        $0 = substr($0, 1, 21) chain " " residue substr($0, 27)
      }
    }
    # Print all lines (modified or unmodified)
    { print }
  ' "$input_file" > "$output_file"

  # Check if formatting issues were fixed
  if diff -q "$input_file" "$output_file" > /dev/null; then
    echo "This script is formatted correctly."
    rm "$output_file" # Remove the output file if no changes were made
  else
    echo "Formatting issues detected and fixed. Updated file: $output_file"
  fi
}

# Run the function
fix_pdb