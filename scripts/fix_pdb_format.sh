#!/bin/bash
# Oct 2025

# Fix PDB format col 5 and col 6 spacing (if the chain ID and residue number are not properly spaced)
# How: If the chain and residue number are merged (e.g., A1000), it separates them correctly (e.g., A 1000).

# Output:
#   Modifies the original input file, but saves the original file with a "_old" extension to it
#   If no formatting issues are found, it prints This script is formatted correctly. and exits without making edits.

# Example:
    # before:
        # ATOM   7519  N   ARG A1000     229.060 228.734 335.474  1.00172.02           N
        # ATOM   7520  CA  ARG A1000     229.625 227.596 334.766  1.00172.02           C
    # after:
        # ATOM   7519  N   ARG A 1000     229.060 228.734 335.474  1.00172.02           N
        # ATOM   7520  CA  ARG A 1000     229.625 227.596 334.766  1.00172.02           C

# Usage:
# /fix_pdb_format.sh input.pdb output.pdb

#!/bin/bash

# Input file
input_file="$1"
output_file="${input_file}"   # Overwrite original file if fixed
backup_file="${input_file}_old"  # Backup file if modifications are made

# Ensure the script takes exactly one input argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <input_file>"
fi

# Check if input file exists
if [ ! -f "$input_file" ]; then
  echo "Error: File '$input_file' not found."
fi

# Function to check and fix formatting
fix_pdb() {
  awk '
    # Match ATOM or HETATM lines
    /^ATOM|^HETATM/ {
      chain_and_residue = substr($0, 22, 5)  # Extract chain and residue columns (22-26)
      chain = substr(chain_and_residue, 1, 1)  # Extract chain (first character)
      residue = substr(chain_and_residue, 2)   # Extract residue number (remaining characters)

      # Check if residue number starts with a digit (merged column issue)
      if (residue ~ /^[0-9]+$/) {
        # Fix formatting by inserting a space between chain and residue
        $0 = substr($0, 1, 21) chain " " residue substr($0, 27)
        fixed = 1  # Mark that a fix has been made
      }
    }
    { print }  # Print all lines (modified or unmodified)
    END {
      # Output whether fixes were made
      if (fixed) exit 2  # Exit with code 2 if fixes were made
    }
  ' "$input_file"
}

# Run the fix function and capture its output
fix_output=$(fix_pdb 2>/dev/null)
fix_status=$?

# If no fixes were needed
if [ "$fix_status" -eq 0 ]; then
  echo "This script is formatted correctly.!!"
fi

# If fixes were made
if [ "$fix_status" -eq 2 ]; then
  # Back up the original file
  mv "$input_file" "$backup_file"

  # Write fixed content to the original file path
  echo "$fix_output" > "$output_file"

  echo "Formatting issues detected and fixed."
  echo "Original file backed up as: $backup_file"
  echo "Modified file saved as: $output_file"
fi

# If something else went wrong
echo "An unexpected error occurred."
