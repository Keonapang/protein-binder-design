#!/bin/bash

# Calculate binding free energy using PRODIGY
# Inputs: PDB file of the peptide binder and target sequence complex
# Output: Binding free energy (ΔG) and dissociation constant (Kd)

# Usage: ./calc_prodigy.sh ${chain} ${start_pos} ${end_pos} ${diffusion} ${temp} ${num_seq} ${i}

chain=$1
start_pos=$2
end_pos=$3
diffusion=$4
temp=$5
num_seq=$6
i=$7
target_sequence=$8
raw_pdb=${9:-NA} # # Check if the 9th argument is provided; if not, default to "NA"

# Check if at least 8 arguments are provided
if [[ $# -lt 8 ]]; then
    echo "At least 8 arguments are required."
    echo "Usage: $0 <chain> <start_pos> <end_pos> <diffusion> <temp> <num_seq> <i> <target_sequence> [raw_pdb]"
    exit 1
fi

REPO_DIR="/home/ubuntu/protein-binder-design"

cycle="target_${chain}${start_pos}_${end_pos}"
params="${diffusion}diff_${temp}temp"

target_pdb="${REPO_DIR}/input/${cycle}.pdb" # Original PDB input file
if [ ! -f "$target_pdb" ]; then
    echo "Target PDB not found: $target_pdb"
    exit 1
fi

if [ "$raw_pdb" == "NA" ]; then
    for iteration in $(seq 1 $i); do
        for num in $(seq 1 $num_seq); do

            # Get peptide binder sequence from proteinPMNN output
            pmnn_file="${REPO_DIR}/${cycle}/3_${cycle}_${params}_i${iteration}.fasta"
            start_line=$((num * 2 + 1))
            end_line=$((num * 2 + 2))
            pmnn_result=$(sed -n "${start_line},${end_line}p" "$pmnn_file")

            binder_pdb="${REPO_DIR}/${cycle}/4_${cycle}_${params}_binder_i${iteration}_${num}.pdb"
            DIR_OUT="${REPO_DIR}/${cycle}"

            # Check if $binder_pdb exist
            if [ ! -f "$binder_pdb" ]; then
                echo "Binder PDB not found: $binder_pdb"
                continue 
            fi

            # 1. Build binder-target PDB complex, output PDB structure
            # python3.11 "${REPO_DIR}/scripts/conversion.py" "$target_pdb" "$binder_pdb" "$diffusion" "$temp" "$DIR_OUT"
            python3.11 "${REPO_DIR}/scripts/align_complex.py" "$target_pdb" "$binder_pdb" "$diffusion" "$temp" "$DIR_OUT"

            # 2. Run PRODIGY to predict binding affinity (kcal.mol-1)
            multi_model_file="${DIR_OUT}/5_${cycle}_${params}_binder_i${iteration}_${num}.pdb"

            output_file="${multi_model_file%.pdb}.txt"
            prodigy_output=$(prodigy "$multi_model_file" -np 4)
            binding_affinity=$(echo "$prodigy_output" | grep "Predicted binding affinity" | awk '{print $6}')
            echo "$binding_affinity" > "$output_file"
            echo "Predicted binding affinity ($binding_affinity kcal.mol-1)"

            # Extract all the necessary values from the output
            binding_affinity=$(echo "$prodigy_output" | grep "Predicted binding affinity" | awk '{print $6}')
            diss_constant=$(echo "$prodigy_output" | grep "Predicted dissociation constant" | awk '{print $8}')
            num_contacts=$(echo "$prodigy_output" | grep "No. of intermolecular contacts" | awk '{print $6}')
            num_charged_charged=$(echo "$prodigy_output" | grep "No. of charged-charged contacts" | awk '{print $6}')
            num_charged_polar=$(echo "$prodigy_output" | grep "No. of charged-polar contacts" | awk '{print $6}')
            num_charged_apolar=$(echo "$prodigy_output" | grep "No. of charged-apolar contacts" | awk '{print $6}')
            num_polar_polar=$(echo "$prodigy_output" | grep "No. of polar-polar contacts" | awk '{print $6}')
            num_apolar_polar=$(echo "$prodigy_output" | grep "No. of apolar-polar contacts" | awk '{print $6}')
            num_apolar_apolar=$(echo "$prodigy_output" | grep "No. of apolar-apolar contacts" | awk '{print $6}')
            percent_apolar_nis=$(echo "$prodigy_output" | grep "Percentage of apolar NIS residues" | awk '{print $6}')
            percent_charged_nis=$(echo "$prodigy_output" | grep "Percentage of charged NIS residues" | awk '{print $6}')

            # 3. Update PDB with comments
            {
                echo "# Final PDB complex of target protein and designed peptide binder";
                echo "# Date: $(date +%Y-%m-%d)";
                echo "# Target sequence: ${target_sequence}";
                echo "# Parameters: ";
                echo "#     target_pdb=${target_pdb}";
                echo "#     binder_pdb=${binder_pdb}";
                echo "#     diffusion=${diffusion}";
                echo "#     temp=${temp}";
                echo "#     iteration=${iteration} RFDiffusion candidates";
                echo "#     num=${num} seqs from ProteinMPNN per candidate";
                echo "# ";
                echo "# Predicted peptide binder from ProteinPMNN:";
                echo "# ${pmnn_result}";
                echo "# PRODIGY results: ";
                echo "#     Predicted binding affinity (kcal.mol-1): ${binding_affinity}";
                echo "#     Predicted dissociation constant (M) at 25.0˚C: ${diss_constant}";
                echo "#     No. of intermolecular contacts: ${num_contacts}";
                echo "#     No. of charged-charged contacts: ${num_charged_charged}";
                echo "#     No. of charged-polar contacts: ${num_charged_polar}";
                echo "#     No. of charged-apolar contacts: ${num_charged_apolar}";
                echo "#     No. of polar-polar contacts: ${num_polar_polar}";
                echo "#     No. of apolar-polar contacts: ${num_apolar_polar}";
                echo "#     No. of apolar-apolar contacts: ${num_apolar_apolar}";
                echo "#     Percentage of apolar NIS residues: ${percent_apolar_nis}";
                echo "#     Percentage of charged NIS residues: ${percent_charged_nis}";
                echo "# ";
                cat "$multi_model_file";
            } > "$multi_model_file".tmp && mv "$multi_model_file".tmp "$multi_model_file"
        done
    done

else

    for iteration in $(seq 1 $i); do
    for num in $(seq 1 $num_seq); do
        for pdb_file in "$target_pdb" "$raw_pdb"; do

            # Get peptide binder sequence from proteinPMNN output
            pmnn_file="${REPO_DIR}/${cycle}/3_${cycle}_${params}_i${iteration}.fasta"
            start_line=$((num * 2 + 1))
            end_line=$((num * 2 + 2))
            pmnn_result=$(sed -n "${start_line},${end_line}p" "$pmnn_file")

            binder_pdb="${REPO_DIR}/${cycle}/4_${cycle}_${params}_binder_i${iteration}_${num}.pdb"
            DIR_OUT="${REPO_DIR}/${cycle}"

            if [ ! -f "$pdb_file" ]; then
                    echo "Target PDB not found: $pdb_file"
                    continue
            fi
            if [ ! -f "$binder_pdb" ]; then
                    echo "Binder PDB not found: $binder_pdb"
                    continue 
            fi

            # 1. Build binder-target PDB complex
            python3.11 "${REPO_DIR}/scripts/conversion.py" "$pdb_file" "$binder_pdb" "$diffusion" "$temp" "$DIR_OUT"

            # 2. Run PRODIGY to predict binding affinity (kcal.mol-1)
            multi_model_file="${DIR_OUT}/5_${cycle}_${params}_binder_i${iteration}_${num}.pdb"
            if [ "$pdb_file" == "$raw_pdb" ]; then
                multi_model_file="${DIR_OUT}/5_${cycle}_${params}_binder_i${iteration}_${num}_raw.pdb"
            fi

            output_file="${multi_model_file%.pdb}.txt"
            prodigy_output=$(prodigy "$multi_model_file" -np 4)
            binding_affinity=$(echo "$prodigy_output" | grep "Predicted binding affinity" | awk '{print $6}')
            echo "$binding_affinity" > "$output_file"
            echo "Predicted binding affinity ($binding_affinity kcal.mol-1)"

            # Extract all the necessary values from the output
            binding_affinity=$(echo "$prodigy_output" | grep "Predicted binding affinity" | awk '{print $6}')
            diss_constant=$(echo "$prodigy_output" | grep "Predicted dissociation constant" | awk '{print $8}')
            num_contacts=$(echo "$prodigy_output" | grep "No. of intermolecular contacts" | awk '{print $6}')
            num_charged_charged=$(echo "$prodigy_output" | grep "No. of charged-charged contacts" | awk '{print $6}')
            num_charged_polar=$(echo "$prodigy_output" | grep "No. of charged-polar contacts" | awk '{print $6}')
            num_charged_apolar=$(echo "$prodigy_output" | grep "No. of charged-apolar contacts" | awk '{print $6}')
            num_polar_polar=$(echo "$prodigy_output" | grep "No. of polar-polar contacts" | awk '{print $6}')
            num_apolar_polar=$(echo "$prodigy_output" | grep "No. of apolar-polar contacts" | awk '{print $6}')
            num_apolar_apolar=$(echo "$prodigy_output" | grep "No. of apolar-apolar contacts" | awk '{print $6}')
            percent_apolar_nis=$(echo "$prodigy_output" | grep "Percentage of apolar NIS residues" | awk '{print $6}')
            percent_charged_nis=$(echo "$prodigy_output" | grep "Percentage of charged NIS residues" | awk '{print $6}')

            # 3. Update PDB with comments
            {
                echo "# Final PDB complex of target protein and designed peptide binder";
                echo "# Date: $(date +%Y-%m-%d)";
                echo "# Target sequence: ${target_sequence}";
                echo "# Parameters: ";
                echo "#     target_pdb=${pdb_file}";
                echo "#     binder_pdb=${binder_pdb}";
                echo "#     diffusion=${diffusion}";
                echo "#     temp=${temp}";
                echo "#     iteration=${iteration} RFDiffusion candidates";
                echo "#     num=${num} seqs from ProteinMPNN per candidate";
                echo "# ";
                echo "# Predicted peptide binder from ProteinPMNN:";
                echo "# ${pmnn_result}";
                echo "# PRODIGY results: ";
                echo "#     Predicted binding affinity (kcal.mol-1): ${binding_affinity}";
                echo "#     Predicted dissociation constant (M) at 25.0˚C: ${diss_constant}";
                echo "#     No. of intermolecular contacts: ${num_contacts}";
                echo "#     No. of charged-charged contacts: ${num_charged_charged}";
                echo "#     No. of charged-polar contacts: ${num_charged_polar}";
                echo "#     No. of charged-apolar contacts: ${num_charged_apolar}";
                echo "#     No. of polar-polar contacts: ${num_polar_polar}";
                echo "#     No. of apolar-polar contacts: ${num_apolar_polar}";
                echo "#     No. of apolar-apolar contacts: ${num_apolar_apolar}";
                echo "#     Percentage of apolar NIS residues: ${percent_apolar_nis}";
                echo "#     Percentage of charged NIS residues: ${percent_charged_nis}";
                echo "# ";
                cat "$multi_model_file";
            } > "$multi_model_file".tmp && mv "$multi_model_file".tmp "$multi_model_file"

            # Compare results from "$target_pdb" and "$raw_pdb"
            if [ "$pdb_file" == "$target_pdb" ]; then
                target_binding_affinity="$binding_affinity"
                target_diss_constant="$diss_constant"
            elif [ "$pdb_file" == "$raw_pdb" ]; then
                raw_binding_affinity="$binding_affinity"
                raw_diss_constant="$diss_constant"
            fi
            done

            # Print comparison results to terminal
            echo "Iteration: $iteration, Num: $num"
            echo "--------------------------------------------"
            echo "Target PDB Results:"
            echo "    Binding Affinity: ${target_binding_affinity} kcal/mol"
            echo "    Dissociation Constant: ${target_diss_constant} M"
            echo "Raw PDB Results:"
            echo "    Binding Affinity: ${raw_binding_affinity} kcal/mol"
            echo "    Dissociation Constant: ${raw_diss_constant} M"
            echo "--------------------------------------------"

            # OPTIONAL: Add conditional analysis (e.g., higher binding affinity)
            if (( $(echo "$target_binding_affinity < $raw_binding_affinity" | bc -l) )); then
                echo "Target PDB has a stronger binding affinity."
            else
                echo "Raw PDB has a stronger binding affinity."
            fi
        done
    done
fi

