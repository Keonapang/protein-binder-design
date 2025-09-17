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

REPO_DIR="/home/ubuntu/protein-binder-design"

cycle="target_${chain}${start_pos}_${end_pos}"
params="${diffusion}diff_${temp}temp"

for iteration in $(seq 1 $i); do
        for num in $(seq 1 $num_seq); do
            target_pdb="${REPO_DIR}/input/${cycle}.pdb" # Original PDB input file
            binder_pdb="${REPO_DIR}/${cycle}/4_${cycle}_${params}_binder_i${iteration}_${num}.pdb"
            DIR_OUT="${REPO_DIR}/${cycle}"

            # Check if $target_pdb and $binder_pdb exist
            if [ ! -f "$target_pdb" ]; then
                echo "Target PDB not found: $target_pdb"
                continue
            fi
            if [ ! -f "$binder_pdb" ]; then
                echo "Binder PDB not found: $binder_pdb"
                continue 
            fi

            # 1. Build binder-target PDB complex
            python3.11 "${REPO_DIR}/scripts/conversion.py" "$target_pdb" "$binder_pdb" "$diffusion" "$temp" "$DIR_OUT"

            # 2. Run PRODIGY to predict binding affinity (kcal.mol-1)
            multi_model_file="${DIR_OUT}/5_${cycle}_${params}_binder_i${iteration}_${num}.pdb"

            # Run PRODIGY and extract binding affinity
            output_file="${multi_model_file%.pdb}.txt"
            prodigy_output=$(prodigy "$multi_model_file" -np 4)
            binding_affinity=$(echo "$prodigy_output" | grep "Predicted binding affinity" | awk '{print $6}')
            echo "$binding_affinity" > "$output_file"
            echo "Predicted binding affinity ($binding_affinity kcal.mol-1)"
        done
done