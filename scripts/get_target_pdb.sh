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
    local input_file="$1"
    local output_file="$2"
    local chain="$3"
    local start_pos="$4"
    local end_pos="$5"

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
            # Print each ATOM line with proper space-delimited formatting
            printf "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f  %5.2f %5.2f           %s\n",
                atom_counter, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12;
            atom_counter++;
        }
    }
    ' "$input_file" >> "$output_file"

    # Add the "TER" line (space-delimited)
    last_atom=$(awk '/^ATOM/ {last_atom = $2} END {print last_atom}' "$output_file")  # Get last atom number
    last_res=$(awk '/^ATOM/ {last_res = $6} END {print last_res}' "$output_file")    # Get last residue number
    printf "TER   %5d      %3s %1s%4d\n" "$((last_atom + 1))" "PRO" "$chain" "$last_res" >> "$output_file"

    # Add the "ENDMDL" and "END" lines (space-delimited)
    printf "ENDMDL\n" >> "$output_file"
    printf "END\n" >> "$output_file"

    # Print the result
    echo "Results written to $output_file"
}

# get_target_pdb() {
#     # Arguments: input file, output file, chain, start residue position, end residue position
#     local input_file="$1"
#     local output_file="$2"
#     local chain="$3"
#     local start_pos="$4"
#     local end_pos="$5"

#     # Create the output file and write the "MODEL 1" header (tab-delimited)
#     printf "MODEL\t\t1\n" > "$output_file"

#     # Process the input PDB file
#     awk -v chain="$chain" -v start="$start_pos" -v end="$end_pos" '
#     BEGIN {
#         atom_counter = 1;  # Counter for atom index
#         res_offset = 1 - start;  # Adjust residue numbers to start from 1
#     }
#     # Only process rows starting with "ATOM" and matching the specified chain
#     $1 == "ATOM" && $5 == chain {
#         res_id = $6 + res_offset;  # Adjust residue number
#         if ($6 >= start && $6 <= end) {
#             # Print updated line with renumbered atom and residue indices (tab-delimited)
#             printf "ATOM\t%6d\t%-4s\t%-3s\t%-1s\t%3d\t%8.3f\t%8.3f\t%8.3f\t%6.2f\t%6.2f\t%-2s\n",
#                 atom_counter, $3, $4, $5, res_id, $7, $8, $9, $10, $11, $12;
#             atom_counter++;
#         }
#     }
#     ' "$input_file" >> "$output_file"
#     echo " "
#     cat "Results: $output_file"

#     # Add the "TER" line (tab-delimited)
#     last_atom=$(awk '/^ATOM/ {last_atom = $2} END {print last_atom}' "$output_file")  # Get last atom number
#     last_res=$(awk '/^ATOM/ {last_res = $6} END {print last_res}' "$output_file")    # Get last residue number
#     printf "TER\t\t%6d\t\tPRO\tA\t%3d\n" "$((last_atom + 1))" "$last_res" >> "$output_file"

#     # Add the "ENDMDL" and "END" lines (tab-delimited)
#     printf "ENDMDL\n" >> "$output_file"
#     printf "END\n" >> "$output_file"
# }
