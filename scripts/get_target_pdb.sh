#!/usr/bin/env bash
# Command-line tool for PRODIGY

# chmod +x /home/ubuntu/protein-binder-design/scripts/get_target_pdb.sh
# get_target_pdb ${input_file} ${output_file} ${chain} ${start_pos} ${end_pos}

# start_pos=60
# end_pos=90
# chain="A"
# input_file="/home/ubuntu/protein-binder-design/input/pdb2e7a.pdb"
# output_file="/home/ubuntu/protein-binder-design/input/target_${chain}${start_pos}_${end_pos}.pdb"

get_target_pdb() {
    # Arguments: input file, output file, chain, start residue position, end residue position
    input_file="$1"
    output_file="$2"
    chain="$3"
    start_pos="$4"
    end_pos="$5"

    # Check if start and end positions exist in column 6 of the input PDB file for the given chain
    valid_positions=$(awk -v chain="$chain" -v start="$start_pos" -v end="$end_pos" '
    $1 == "ATOM" && $5 == chain && ($6 == start || $6 == end) { print $6 }
    ' "$input_file" | sort -nu)

    # Check if start and end positions are valid
    if ! echo "$valid_positions" | grep -q "$start_pos"; then
        echo "Error: Start position $start_pos not found in column 6 for chain $chain in $input_file."
        return 1
    fi

    if ! echo "$valid_positions" | grep -q "$end_pos"; then
        echo "Error: End position $end_pos not found in column 6 for chain $chain in $input_file."
        return 1
    fi

    # Create the output file and write the "MODEL 1" header (space-delimited)
    printf "MODEL        1\n" > "$output_file"

    # Process the input PDB file
    awk -v chain="$chain" -v start="$start_pos" -v end="$end_pos" '
    BEGIN {
        atom_counter = 1;  # Counter for atom indices
    }
    # Only process rows starting with "ATOM" and matching the specified chain
    $1 == "ATOM" && $5 == chain {
        if ($6 >= start && $6 <= end) {
            # Save the current row details for later use
            last_row = $0;
            last_atom_number = atom_counter;
            last_residue = $4;  # Residue name
            last_residue_number = $6;  # Residue number

            # Print each ATOM line with proper space-delimited formatting
            printf "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f  %5.2f %5.2f           %s\n",
                atom_counter, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12;
            atom_counter++;
        }
    }
    END {
        # Print the TER line with the last residue details
        printf "TER   %5d      %3s %1s%4d\n", last_atom_number + 1, last_residue, chain, last_residue_number;
    }
    ' "$input_file" >> "$output_file"

    # Add the "ENDMDL" and "END" lines (space-delimited)
    printf "ENDMDL\n" >> "$output_file"
    printf "END\n" >> "$output_file"

    echo " "
    echo "Target seq PDB: $output_file"
}

# get_target_pdb() {
#     # Arguments: input file, output file, chain, start residue position, end residue position
#     input_file="$1"
#     output_file="$2"
#     chain="$3"
#     start_pos="$4"
#     end_pos="$5"

#     # Create the output file and write the "MODEL 1" header (space-delimited)
#     printf "MODEL        1\n" > "$output_file"

#     # Process the input PDB file
#     awk -v chain="$chain" -v start="$start_pos" -v end="$end_pos" '
#     BEGIN {
#         atom_counter = 1;  # Counter for atom indices
#     }
#     # Only process rows starting with "ATOM" and matching the specified chain
#     $1 == "ATOM" && $5 == chain {
#         if ($6 >= start && $6 <= end) {
#             # Save the current row details for later use
#             last_row = $0;
#             last_atom_number = atom_counter;
#             last_residue = $4;  # Residue name
#             last_residue_number = $6;  # Residue number

#             # Print each ATOM line with proper space-delimited formatting
#             printf "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f  %5.2f %5.2f           %s\n",
#                 atom_counter, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12;
#             atom_counter++;
#         }
#     }
#     END {
#         # Print the TER line with the last residue details
#         printf "TER   %5d      %3s %1s%4d\n", last_atom_number + 1, last_residue, chain, last_residue_number;
#     }
#     ' "$input_file" >> "$output_file"

#     # Add the "ENDMDL" and "END" lines (space-delimited)
#     printf "ENDMDL\n" >> "$output_file"
#     printf "END\n" >> "$output_file"

#     echo " "
#     echo "Target seq PDB: $output_file"
# }

# get_target_pdb() {
#     # Arguments: input file, output file, chain, start residue position, end residue position
#     input_file="$1"
#     output_file="$2"
#     chain="$3"
#     start_pos="$4"
#     end_pos="$5"

#     # Create the output file and write the "MODEL 1" header (space-delimited)
#     printf "MODEL        1\n" > "$output_file"

#     # Process the input PDB file
#     awk -v chain="$chain" -v start="$start_pos" -v end="$end_pos" '
#     BEGIN {
#         atom_counter = 1;  # Counter for atom indices
#     }
#     # Only process rows starting with "ATOM" and matching the specified chain
#     $1 == "ATOM" && $5 == chain {
#         if ($6 >= start && $6 <= end) {
#             # Print each ATOM line with proper space-delimited formatting
#             printf "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f  %5.2f %5.2f           %s\n",
#                 atom_counter, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12;
#             atom_counter++;
#         }
#     }
#     ' "$input_file" >> "$output_file"

#     # Add the "TER" line (space-delimited)
#     last_atom=$(awk '/^ATOM/ {last_atom = $2} END {print last_atom}' "$output_file")  # Get last atom number
#     last_res=$(awk '/^ATOM/ {last_res = $6} END {print last_res}' "$output_file")    # Get last residue number
#     printf "TER   %5d      %3s %1s%4d\n" "$((last_atom + 1))" "PRO" "$chain" "$last_res" >> "$output_file"

#     # Add the "ENDMDL" and "END" lines (space-delimited)
#     printf "ENDMDL\n" >> "$output_file"
#     printf "END\n" >> "$output_file"

#     echo " "
#     echo "Target seq PDB: $output_file"
# }
