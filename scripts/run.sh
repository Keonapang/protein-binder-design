#!/usr/bin/env bash
# Command-line to run the generative protein binder design pipeline. 
# Dates: Apr 24, 2025 - March 17, 2026

# start_pos=60
# end_pos=90
# chain="A"
# input_file="/home/ubuntu/protein-binder-design/input/pdb2e7a.pdb"
# output_file="/home/ubuntu/protein-binder-design/input/target_${chain}${start_pos}_${end_pos}.pdb"


# 1. Install Dependencies
sudo apt-get update # updated nvidia toolkit
sudo apt-get install -y docker-compose # docker compose version 2+
sudo apt install python3.11
mkdir -p ~/.cache/nim
chmod -R 777 ~/.cache/nim    
export HOST_NIM_CACHE=~/.cache/nim

export NGC_CLI_API_KEY=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
docker login nvcr.io --username='$oauthtoken' --password="${NGC_CLI_API_KEY}"
cd ~
cd protein-binder-design/deploy/
docker compose up

###########################################################################################################
# **Set up `openmm` and `pdbfier` to align backbones (requires python 3.13):**
###########################################################################################################
# cd ~
# git clone https://github.com/openmm/pdbfixer
# cd pdbfixer
# python setup.py install
# python3.13 -m pip install numpy prodigy-prot torch Bio biopython pdb-tools

# # Install conda
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash Miniconda3-latest-Linux-x86_64.sh # installed to $HOME/miniconda3
# echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
# source ~/.bashrc # conda --version
# conda create -n pdbfixer_env python=3.13 -y
# conda init
# # Open new temrinal window(s) and install these packages:
# conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
# conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
# conda activate pdbfixer_env
# conda install -c conda-forge pdbfixer
# conda install -c conda-forge openmm
# python3.13 -m pip install numpy prodigy-prot torch Bio biopython pdb-tools

# export NGC_CLI_API_KEY=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
# docker login nvcr.io --username='$oauthtoken' --password="${NGC_CLI_API_KEY}"

###########################################################################################################
curl localhost:8082/v1/health/ready # RFdiffusion
curl localhost:8083/v1/health/ready # Protein MPNN
    
# Define protein
protein="1TNF" # 1TNF, apob, tnf
REPO_DIR="$HOME/protein-binder-design"
raw_pdb="${REPO_DIR}/input/${protein}.pdb"             # target protein
input_file="${REPO_DIR}/input/target_file_${protein}_surface.txt"  # chain, hotspot residue, start/end pos 

# Clean up raw  format
REPO_DIR="$HOME/protein-binder-design"
sed -i 's/\r$//' "${REPO_DIR}/scripts/get_target_pdb.sh" # optional (to remove any hidden spaces from Windows)
sed -i 's/\r$//' "${REPO_DIR}/scripts/5_run_prodigy.sh" 
sed -i 's/^[ \t]*//;s/[ \t]*$//' "${REPO_DIR}/scripts/5_run_prodigy.sh"
bash "${REPO_DIR}/scripts/fix_pdb_format.sh" "$raw_pdb"

# Convert residue position file to tab-delimited
awk '{$1=$1; gsub(" ", "\t"); print}' "$input_file" > "$input_file.tmp" && mv "$input_file.tmp" "$input_file"
sed -i 's/\r$//' "$input_file"
cat -A "$input_file"
wc -l $input_file
wc -l $raw_pdb

sed -n '13p' "$input_file" | while IFS=$'\t' read -r chain hotspot_res_prefix start_pos end_pos; do
    # Define script input variables
    diffusion=50
    temp=0.3
    i=4
    num_seq=2
    peptide_length="15-25" # set a range (i.e."15-25") or a value ("25")
    chain="A" 

   # Variables (do not modify)
    hotspot_res="${chain}${hotspot_res_prefix}"
    contigs="${chain}${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
    name="target_${chain}${start_pos}_${end_pos}"
    params="${diffusion}diff_${temp}temp"
    echo "Processing chain=${chain}, hotspot_res=${hotspot_res}, start_pos=${start_pos}, end_pos=${end_pos}"
    echo "name=${name},    params=${params},    contigs=${contigs}"
    echo ""
    echo "=========================================================================="
    export CUDA_VISIBLE_DEVICES=0
    unset hotspot_res # NO HOTSPOTS SET

    # Step 1: Build target structure PDB and extract target seq amino acid
    target_pdb="${REPO_DIR}/input/${name}.pdb"
    bash ${REPO_DIR}/scripts/get_target_pdb.sh "${raw_pdb}" "${target_pdb}" "${chain}" "${start_pos}" "${end_pos}"
    if [ -f "$target_pdb" ]; then target_sequence=$(bash "${REPO_DIR}/scripts/get_target_seq.sh" "${target_pdb}"); fi

    # Step 2: Run the protein binder design script
    python3.11 "${REPO_DIR}/scripts/3_protein_binder_design.py" --root "${REPO_DIR}" \
    --num_seq "${num_seq}" --diffusion "${diffusion}" --temp "${temp}" --target_sequence "${target_sequence}" \
    --contigs "${contigs}" --i "${i}" --hotspot_res "${hotspot_res}" --target_pdb "${target_pdb}" --chain "${chain}"

    # Step 3: Generate merged binding alignment for peptide-target protein, and optimize alignment 
    python3.13 "${REPO_DIR}/scripts/4_merge_seq_to_backbone.py" "${REPO_DIR}" "${chain}" "${i}" "${num_seq}" "${name}" "${params}" --solvent
    bash "${REPO_DIR}/scripts/5_run_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${i}" "${num_seq}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}" "${input_file}" "${hotspot_res}"
done
# ---------------- Clean up if necessary 
# # remove all files from ${REPO_DIR}/${name} directory that begin with "5_target_" and end with ".pdb"
# ls ./5_target_*.pdb
# rm ./5_target_*.pdb
# # rename all files that begin with "5_target_" and end with "_complex.pdb_old" to remove the "_old" suffix
# for file in 5_target_*_complex.pdb_old; do
#     mv "$file" "${file%_old}"
# done

########### alternative code: Nov 1 ##########
# Define protein
protein="1TNF" # 1TNF, apob, tnf
REPO_DIR="/home/shadeform/protein-binder-design"
raw_pdb="${REPO_DIR}/input/${protein}.pdb"             # target protein
input_file="${REPO_DIR}/input/target_file_${protein}_surface.txt"  # chain, hotspot residue, start/end pos 
export NGC_CLI_API_KEY=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
docker login nvcr.io --username='$oauthtoken' --password="${NGC_CLI_API_KEY}"
export CUDA_VISIBLE_DEVICES=1

sed -n '9,21p' "$input_file" | while IFS=$'\t' read -r hotspot_res_prefix; do
    # Define script input variables
    diffusion=50
    temp=0.3
    i=4
    num_seq=2
    peptide_length="15-25" # set a range (i.e."15-25") or a value ("25")
    chain="A" 

   # Variables (do not modify)
    start_pos="6"
    end_pos="157"
    hotspot_res=${hotspot_res_prefix}
    name="target${hotspot_res}"

    contigs="${chain}${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
    params="${diffusion}diff_${temp}temp"
    echo "Processing chain=${chain}, hotspot_res=${hotspot_res}, start_pos=${start_pos}, end_pos=${end_pos}"
    echo "name=${name},    params=${params},    contigs=${contigs}"
    echo "=========================================================================="
    # unset hotspot_res # NO HOTSPOTS SET
    # Step 1: Build target structure PDB and extract target seq amino acid
    target_pdb="${REPO_DIR}/input/${name}.pdb"
    # bash ${REPO_DIR}/scripts/get_target_pdb.sh "${raw_pdb}" "${target_pdb}" "${chain}" "${start_pos}" "${end_pos}"
    if [ -f "$target_pdb" ]; then target_sequence=$(bash "${REPO_DIR}/scripts/get_target_seq.sh" "${target_pdb}"); fi
    echo "Target sequence: ${target_sequence}"
    # # Step 2: Run the protein binder design script
    # python3.11 "${REPO_DIR}/scripts/3_protein_binder_design.py" --root "${REPO_DIR}" \
    # --num_seq "${num_seq}" --diffusion "${diffusion}" --temp "${temp}" --target_sequence "${target_sequence}" \
    # --contigs "${contigs}" --i "${i}" --hotspot_res ${hotspot_res} --target_pdb "${target_pdb}" --chain "${chain}"

    # # Step 3: Generate merged binding alignment for peptide-target protein, and optimize alignment 
    # python3.13 "${REPO_DIR}/scripts/4_merge_seq_to_backbone.py" "${REPO_DIR}" "${chain}" "${i}" "${num_seq}" "${name}" "${params}" --solvent
    bash "${REPO_DIR}/scripts/5_run_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${i}" "${num_seq}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}" "${input_file}" "${hotspot_res}"
done

####################################################################
# CODE WORKS!
# Creates summary file .txt
####################################################################
# Initialize the summary file with headers
summary_file="${REPO_DIR}/summary_${protein}.txt"
echo -e "Target\titeration\tnum_seq\tbinding_affinity\tdiss_constant" > "$summary_file"

# Read input file and process for each chain, hotspot, start, and end position
sed -n '1,21p' "$input_file" | while IFS=$'\t' read -r hotspot_res_prefix; do
    diffusion=50
    temp=0.3
    i=4
    num_seq=2
    chain="A" 
    hotspot=${hotspot_res_prefix}
    for iteration in $(seq 1 $i); do
        for num in $(seq 1 $num_seq); do
            start_pos="6"
            end_pos="157"
            name="target${hotspot}"
            params="${diffusion}diff_${temp}temp"
            aligned_pdb="${REPO_DIR}/${name}/5_${name}_${params}_i${iteration}_${num}_complex.pdb"

            if [[ -f "$aligned_pdb" ]]; then
                binding_affinity=$(grep "# Binding affinity (kcal.mol-1): " "$aligned_pdb" | awk -F": " '{print $2}')
                diss_constant=$(grep "# Dissociation constant (Kb) at 25.0˚C: " "$aligned_pdb" | awk -F": " '{print $2}')
                echo -e "${name}\t${iteration}\t${num}\t${binding_affinity}\t${diss_constant}" >> "$summary_file"
                echo "5_${name}_${params}_i${iteration}_${num}_complex.pdb:  diss_constant=${diss_constant}"
            else
                echo -e "${name}\t${iteration}\t${num}\tNA\tNA" >> "$summary_file"
                echo "Warning: File '$aligned_pdb' not found. Skipping this entry."
            fi
            echo ""
            unset binding_affinity
            unset diss_constant
        done
    done
done
echo "Summary file created: $summary_file"

# Read input file and process for each chain, hotspot, start, and end position
# sed -n '1,21p' "$input_file" | while IFS=$'\t' read -r hotspot_res_prefix; do
#         diffusion=50
#         temp=0.3
#         i=4
#         num_seq=2
#         chain="A" 
#         hotspot=${hotspot_res_prefix}
#     for iteration in $(seq 1 $i); do
#         for num in $(seq 1 $num_seq); do
#         start_pos="6"
#         end_pos="157"
#         name="target${hotspot}"
#         params="${diffusion}diff_${temp}temp"
#         aligned_pdb="${REPO_DIR}/${name}/5_${name}_${params}_i${iteration}_${num}_complex.pdb"
#             if [[ -f "$aligned_pdb" ]]; then
#                 binding_affinity=$(grep "# Binding affinity (kcal.mol-1): " "$aligned_pdb" | awk -F": " '{print $2}')
#                 diss_constant=$(grep "# Dissociation constant (Kb) at 25.0˚C: " "$aligned_pdb" | awk -F": " '{print $2}')
#                 echo -e "${name}\t${iteration}\t${num}\t${binding_affinity}\t${diss_constant}" >> "$summary_file"
#             echo "5_${name}_${params}_i${iteration}_${num}_complex.pdb:  diss_constant=${diss_constant}"
#             else
#                 echo "Warning: File '$aligned_pdb' not found. Skipping this entry."
#             fi
#             echo ""
#             unset binding_affinity
#             unset diss_constant
#         done
#     done
# done
# echo "Summary file created: $summary_file"

####### If hotspots are a combination of different residues ##########

summary_file="${REPO_DIR}/summary_${protein}.txt"
echo -e "Target\titeration\tnum_seq\tbinding_affinity\tdiss_constant" > "$summary_file"
sed -n '1,21p' "$input_file" | while IFS=$'\t' read -r chain hotspot_res_prefix start_pos end_pos; do
        diffusion=50
        temp=0.3
        i=4
        num_seq=2
        chain="A" 
    for iteration in $(seq 1 $i); do
        for num in $(seq 1 $num_seq); do

            chain="A"
            hotspot_res="${chain}${hotspot_res_prefix}"
            contigs="${chain}${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
            name="target_${chain}${start_pos}_${end_pos}"
            params="${diffusion}diff_${temp}temp"
            aligned_pdb="${REPO_DIR}/${name}/5_${name}_${params}_i${iteration}_${num}_complex.pdb"

            if [[ -f "$aligned_pdb" ]]; then
                binding_affinity=$(grep "# Binding affinity (kcal.mol-1): " "$aligned_pdb" | awk -F": " '{print $2}')
                diss_constant=$(grep "# Dissociation constant (Kb) at 25.0˚C: " "$aligned_pdb" | awk -F": " '{print $2}')
                echo -e "${name}\t${iteration}\t${num}\t${binding_affinity}\t${diss_constant}" >> "$summary_file"
                echo "Name=${name}, iteration=${iteration}, num=${num}, binding_affinity=${binding_affinity}, diss_constant=${diss_constant}"
            else
                echo "Warning: File '$aligned_pdb' not found. Skipping this entry."
            fi
            echo ""
        done
    done
done
echo "Summary file created: $summary_file"


###################

# Loop through each line of $input_file
while IFS=$'\t' read -r chain hotspot_res_prefix start_pos end_pos; do # space-delimited

    # Variables (do not modify)
    hotspot_res="${chain}${hotspot_res_prefix}"
    contigs="A${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
    name="target_${chain}${start_pos}_${end_pos}"
    params="${diffusion}diff_${temp}temp"
    echo "Processing chain=${chain}, hotspot_res=${hotspot_res}, start_pos=${start_pos}, end_pos=${end_pos}"
    echo "name=${name},    params=${params},    contigs=${contigs}"
    echo ""

    # Step 1: Build target structure PDB and extract target seq amino acid
    target_pdb="${REPO_DIR}/input/${name}.pdb" # output
    bash ${REPO_DIR}/scripts/get_target_pdb.sh "${raw_pdb}" "${target_pdb}" "${chain}" "${start_pos}" "${end_pos}"
    if [ -f "$target_pdb" ]; then target_sequence=$(bash "${REPO_DIR}/scripts/get_target_seq.sh" "${target_pdb}"); fi

    # Step 2: Run the protein binder design script
    python3.11 "${REPO_DIR}/scripts/3_protein_binder_design.py" --root "${REPO_DIR}" \
    --num_seq "${num_seq}" --diffusion "${diffusion}" --temp "${temp}" --target_sequence "${target_sequence}" \
    --contigs "${contigs}" --i "${i}" --hotspot_res "${hotspot_res}" --target_pdb "${target_pdb}" --chain "${chain}"
done < "$input_file"

##################################################################################################################

while IFS=$'\t' read -r chain hotspot_res_prefix start_pos end_pos; do # space-delimited

    # Variables (do not modify)
    hotspot_res="${chain}${hotspot_res_prefix}"
    contigs="A${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
    name="target_${chain}${start_pos}_${end_pos}"
    params="${diffusion}diff_${temp}temp"
    echo "Processing chain=${chain}, hotspot_res=${hotspot_res}, start_pos=${start_pos}, end_pos=${end_pos}"
    echo "name=${name},    params=${params},    contigs=${contigs}"
    echo ""

    # Step 3: Generate merged binding alignment for peptide and target protein, and then optimize alignment 
    python3.13 ${REPO_DIR}/scripts/4_merge_seq_to_backbone.py "${REPO_DIR}" A ${i} ${num_seq} ${name} ${params} --solvent

    # Step 5 alternative: PRODIGY
    bash "${REPO_DIR}/scripts/calc_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${i}" "${num_seq}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}" "${input_file}"

done < "$input_file"

