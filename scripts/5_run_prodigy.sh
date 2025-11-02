#!/bin/bash
######################################################
# Calculate binding free energy using PRODIGY
# Inputs: Two PDB files of the peptide binder and raw target sequence complex
# Output: A single combined PDB file with header notes; 
#       - including predicted peptide sequence,
#       - free energy (ΔG) and dissociation constant (Kd) etc..

# Usage: calc_prodigy.sh ${chain} ${start_pos} ${end_pos} ${diffusion} ${temp} ${i} ${num_seq} ${target_sequence} ${REPO_DIR} ${raw_pdb} $input_file
######################################################
chain=$1
start_pos=$2
end_pos=$3
diffusion=$4
temp=$5
i=$6
num_seq=$7
target_sequence=$8
REPO_DIR=$9
target_pdb=${10}
input_file=${11}
hotspot=${12}
remove_water=TRUE

# At least 11 arguments are provided
if [[ $# -lt 12 ]]; then
    echo "At least 12 arguments are required."
    echo "Usage: $0 {name} {params} {i} {num_seq} {target_sequence} {REPO_DIR} {target_pdb} {input_file} {hotspot}"
    continue
fi

echo "Chain: $chain"
echo "Start position: $start_pos"
echo "hotspot: $hotspot"
echo "End position: $end_pos"
echo "Diffusion steps: $diffusion"
echo "Temperature: $temp"
echo "Number of RFDiffusion iterations: $i"
echo "Number of ProteinMPNN sequences per target: $num_seq"
echo "Target sequence: $target_sequence"
echo "REPO_DIR: $REPO_DIR"
echo "Target PDB: $target_pdb"
echo "Input file: $input_file"
echo ""

# start_pos=6
# end_pos=157
# chain="A"
# diffusion=50
# temp=0.3
# i=4
# num_seq=2
# name="target'A6','A7','A8','A9','A34'" #name="target_${chain}${start_pos}_${end_pos}"
# params="${diffusion}diff_${temp}temp"
# iteration=1
# num=1

# for iteration in $(seq 1 $i); do
#     for num in $(seq 1 $num_seq); do

# REPO_DIR="/home/shadeform/protein-binder-design"
# aligned_pdb=${REPO_DIR}/${name}/5_${name}_${params}_i${iteration}_${num}_complex.pdb
# wc -l "$aligned_pdb"
# wc -l "${REPO_DIR}/${name}/3_${name}_${params}_i${iteration}.fasta"
# prodigy "$aligned_pdb" -np 4 --contact_list
# remove_water=TRUE
#         if [ "$remove_water" = "TRUE" ]; then
#             grep -v "^HETATM" "$aligned_pdb" > temp.pdb && mv temp.pdb "$aligned_pdb"
#         fi
#         grep -v "^#" "$aligned_pdb" > temp.pdb && mv temp.pdb "$aligned_pdb"
#         if grep -q "Final PDB complex of target protein and designed peptide binder" "$aligned_pdb"; then
#             echo "ERROR: The aligned PDB file contains headers already!!"
#             continue
#         fi
#             pmnn_file="${REPO_DIR}/${name}/3_${name}_${params}_i${iteration}.fasta"
#             if [ ! -f "$pmnn_file" ]; then
#                 echo "MISSING: $pmnn_file"
#                 echo ""
#                 continue
#             fi
#             start_line=$((num * 2 + 1))
#             end_line=$((num * 2 + 2))
#             pmnn_result=$(sed -n "${start_line},${end_line}p" "$pmnn_file")
#             prodigy_output=$(prodigy "$aligned_pdb" -np 4 --contact_list)
#             if echo "$prodigy_output" | grep -q -e "Invalid" -e "No contacts"; then
#                 echo "============================================================="
#                 echo ""
#                 echo "$prodigy_output"  # Print the error message from PRODIGY
#                 echo ""
#                 echo "[!] Stopping /5_run_prodigy.sh due to error"
#                 echo ""
#                 echo "============================================================="
#                 return 1  # Gracefully stop the script by returning a status code
#             fi
#             binding_affinity=$(echo "$prodigy_output" | grep "Predicted binding affinity" | awk '{print $6}')
#             diss_constant=$(echo "$prodigy_output" | grep "Predicted dissociation constant" | awk '{print $8}')
#             num_contacts=$(echo "$prodigy_output" | grep "No. of intermolecular contacts" | awk '{print $6}')
#             num_charged_charged=$(echo "$prodigy_output" | grep "No. of charged-charged contacts" | awk '{print $6}')
#             num_charged_polar=$(echo "$prodigy_output" | grep "No. of charged-polar contacts" | awk '{print $6}')
#             num_charged_apolar=$(echo "$prodigy_output" | grep "No. of charged-apolar contacts" | awk '{print $6}')
#             num_polar_polar=$(echo "$prodigy_output" | grep "No. of polar-polar contacts" | awk '{print $6}')
#             num_apolar_polar=$(echo "$prodigy_output" | grep "No. of apolar-polar contacts" | awk '{print $6}')
#             num_apolar_apolar=$(echo "$prodigy_output" | grep "No. of apolar-apolar contacts" | awk '{print $6}')
#             percent_apolar_nis=$(echo "$prodigy_output" | grep "Percentage of apolar NIS residues" | awk '{print $7}')
#             percent_charged_nis=$(echo "$prodigy_output" | grep "Percentage of charged NIS residues" | awk '{print $7}')
#             echo "Binding affinity ($binding_affinity kcal.mol-1)"
#             echo "Dissociation constant ($diss_constant Kb at 25.0˚C)"
#             {
#                 echo "# Final PDB complex of target protein and designed peptide binder";
#                 echo "# Date: $(date +%Y-%m-%d)";
#                 echo "# Target sequence: ${target_sequence}";
#                 echo "# RFDiffusion candidate ${iteration}, ProteinPMNN predicted sequence ${num}";
#                 echo "# ";
#                 echo "# ========= Input parameters ========= ";
#                 echo "# target_pdb: ${target_pdb}";
#                 echo "# input target coordinates: ${input_file}";
#                 echo "# chain=${chain}";
#                 echo "# start_pos=${start_pos}";
#                 echo "# end_pos=${end_pos}";
#                 echo "# hotspot: ${hotspot}";
#                 echo "# diffusion=${diffusion}";
#                 echo "# temp=${temp}";
#                 echo "# iteration=${i} total RFDiffusion candidates";
#                 echo "# num_seq=${num_seq} total ProteinMPNN sequences";
#                 echo "# ";
#                 echo "# ========= ProteinPMNN predicted peptide =========";
#                 echo "$pmnn_result" | sed 's/^/# /';   # Add a "#" in front of every line in pmnn_result
#                 echo "# ";
#                 echo "# ========= PRODIGY results ========= ";
#                 echo "# Binding affinity (kcal.mol-1): ${binding_affinity}";
#                 echo "# Dissociation constant (Kb) at 25.0˚C: ${diss_constant}";
#                 echo "# No. of intermolecular contacts: ${num_contacts}";
#                 echo "# No. of charged-charged contacts: ${num_charged_charged}";
#                 echo "# No. of charged-polar contacts: ${num_charged_polar}";
#                 echo "# No. of charged-apolar contacts: ${num_charged_apolar}";
#                 echo "# No. of polar-polar contacts: ${num_polar_polar}";
#                 echo "# No. of apolar-polar contacts: ${num_apolar_polar}";
#                 echo "# No. of apolar-apolar contacts: ${num_apolar_apolar}";
#                 echo "# Percentage of apolar NIS residues: ${percent_apolar_nis}";
#                 echo "# Percentage of charged NIS residues: ${percent_charged_nis}";
#                 echo "# ";
#                 cat "$aligned_pdb";
#             } > "$aligned_pdb".tmp && mv "$aligned_pdb".tmp "$aligned_pdb"
# done
# done

echo ""
echo "####################################################"
echo "5. Running PRODIGY to calculate dissociation constant"
echo "####################################################"
echo "REPO_DIR: $REPO_DIR"
if [[ "$hotspot_res" == *","* ]]; then
    name="target${hotspot}"
else
    name="target_${chain}${start_pos}_${end_pos}"
fi
params="${diffusion}diff_${temp}temp"
echo "name: $name"
echo ""

for iteration in $(seq 1 $i); do
    for num in $(seq 1 $num_seq); do

        aligned_pdb=${REPO_DIR}/${name}/5_${name}_${params}_i${iteration}_${num}_complex.pdb

        echo "Aligned PDB:  5_${name}_${params}_i${iteration}_${num}_complex.pdb"
        echo ""
        if [ ! -f "$aligned_pdb" ]; then
            echo "MISSING: $aligned_pdb"
            echo ""
            continue
        fi

        # Remove water molecules from PDB
        if [ "$remove_water" = "TRUE" ]; then
            grep -v "^HETATM" "$aligned_pdb" > temp.pdb && mv temp.pdb "$aligned_pdb"
        fi

        # Remove lines starting with '#' and overwrite the original file
        grep -v "^#" "$aligned_pdb" > temp.pdb && mv temp.pdb "$aligned_pdb"

        if grep -q "Final PDB complex of target protein and designed peptide binder" "$aligned_pdb"; then
            echo "ERROR: The aligned PDB file contains headers already!!"
            continue
        fi

        # Re-align pdb formatting (cols 5 and 6)
        if grep -q "A000" "$aligned_pdb"; then
            chmod +x "${REPO_DIR}/scripts/fix_pdb_format.sh"
            bash "${REPO_DIR}/scripts/fix_pdb_format.sh" "$aligned_pdb"
        else
            echo " "
        fi

        # Process the file and fix merged "HETATM" entries
        # awk '
        # # If the line starts with "HETATM" but column 2 is merged with column 1
        # /^HETATM[0-9]/ {
        #     # Separate "HETATM" and the numeric part into two columns
        #     $1 = substr($1, 1, 6) " " substr($1, 7)
        # }
        # { print }
        # ' "$aligned_pdb" > temp.pdb && mv temp.pdb "$aligned_pdb"

            # Get peptide binder sequence from proteinPMNN output
            pmnn_file="${REPO_DIR}/${name}/3_${name}_${params}_i${iteration}.fasta"
            if [ ! -f "$pmnn_file" ]; then
                echo "MISSING: $pmnn_file"
                echo ""
                continue
            fi
            start_line=$((num * 2 + 1))
            end_line=$((num * 2 + 2))
            pmnn_result=$(sed -n "${start_line},${end_line}p" "$pmnn_file")

            # 2. Run PRODIGY to predict binding affinity (kcal.mol-1)
            prodigy_output=$(prodigy "$aligned_pdb" -np 4 --contact_list)

            if echo "$prodigy_output" | grep -q -e "Invalid" -e "No contacts"; then
                echo "============================================================="
                echo ""
                echo "$prodigy_output"  # Print the error message from PRODIGY
                echo ""
                echo "[!] Stopping /5_run_prodigy.sh due to error"
                echo ""
                echo "============================================================="
                return 1  # Gracefully stop the script by returning a status code
            fi

            binding_affinity=$(echo "$prodigy_output" | grep "Predicted binding affinity" | awk '{print $6}')
            diss_constant=$(echo "$prodigy_output" | grep "Predicted dissociation constant" | awk '{print $8}')
            num_contacts=$(echo "$prodigy_output" | grep "No. of intermolecular contacts" | awk '{print $6}')
            num_charged_charged=$(echo "$prodigy_output" | grep "No. of charged-charged contacts" | awk '{print $6}')
            num_charged_polar=$(echo "$prodigy_output" | grep "No. of charged-polar contacts" | awk '{print $6}')
            num_charged_apolar=$(echo "$prodigy_output" | grep "No. of charged-apolar contacts" | awk '{print $6}')
            num_polar_polar=$(echo "$prodigy_output" | grep "No. of polar-polar contacts" | awk '{print $6}')
            num_apolar_polar=$(echo "$prodigy_output" | grep "No. of apolar-polar contacts" | awk '{print $6}')
            num_apolar_apolar=$(echo "$prodigy_output" | grep "No. of apolar-apolar contacts" | awk '{print $6}')
            percent_apolar_nis=$(echo "$prodigy_output" | grep "Percentage of apolar NIS residues" | awk '{print $7}')
            percent_charged_nis=$(echo "$prodigy_output" | grep "Percentage of charged NIS residues" | awk '{print $7}')
            echo "Binding affinity ($binding_affinity kcal.mol-1)"
            echo "Dissociation constant ($diss_constant Kb at 25.0˚C)"

            # 3. Update PDB with comments
            {
                echo "# Final PDB complex of target protein and designed peptide binder";
                echo "# Date: $(date +%Y-%m-%d)";
                echo "# Target sequence: ${target_sequence}";
                echo "# RFDiffusion candidate ${iteration}, ProteinPMNN predicted sequence ${num}";
                echo "# ";
                echo "# ========= Input parameters ========= ";
                echo "# target_pdb: ${target_pdb}";
                echo "# input target coordinates: ${input_file}";
                echo "# chain=${chain}";
                echo "# start_pos=${start_pos}";
                echo "# end_pos=${end_pos}";
                echo "# hotspot: ${hotspot}";
                echo "# diffusion=${diffusion}";
                echo "# temp=${temp}";
                echo "# iteration=${i} total RFDiffusion candidates";
                echo "# num_seq=${num_seq} total ProteinMPNN sequences";
                echo "# ";
                echo "# ========= ProteinPMNN predicted peptide =========";
                echo "$pmnn_result" | sed 's/^/# /';   # Add a "#" in front of every line in pmnn_result
                echo "# ";
                echo "# ========= PRODIGY results ========= ";
                echo "# Binding affinity (kcal.mol-1): ${binding_affinity}";
                echo "# Dissociation constant (Kb) at 25.0˚C: ${diss_constant}";
                echo "# No. of intermolecular contacts: ${num_contacts}";
                echo "# No. of charged-charged contacts: ${num_charged_charged}";
                echo "# No. of charged-polar contacts: ${num_charged_polar}";
                echo "# No. of charged-apolar contacts: ${num_charged_apolar}";
                echo "# No. of polar-polar contacts: ${num_polar_polar}";
                echo "# No. of apolar-polar contacts: ${num_apolar_polar}";
                echo "# No. of apolar-apolar contacts: ${num_apolar_apolar}";
                echo "# Percentage of apolar NIS residues: ${percent_apolar_nis}";
                echo "# Percentage of charged NIS residues: ${percent_charged_nis}";
                echo "# ";
                cat "$aligned_pdb";
            } > "$aligned_pdb".tmp && mv "$aligned_pdb".tmp "$aligned_pdb"
        done
done

# Clean up intermediate .ic files
rm ${REPO_DIR}/${name}/5_*.ic