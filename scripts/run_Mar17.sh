#!/usr/bin/env bash
# Command-line to run the generative protein binder design pipeline. 
# Dates: Apr 24, 2025 - March 17, 2026

# 1. Install Dependencies
git clone https://github.com/Keonapang/protein-binder-design

sudo apt-get update # updated nvidia toolkit
sudo apt-get install -y docker-compose # docker compose version 2+
sudo apt install python3.11

# Create cache directory for the NIMs
mkdir -p ~/.cache/nim
sudo chmod -R 777 ~/.cache/nim    
export HOST_NIM_CACHE=~/.cache/nim # HOST_NIM_CACHE is used in the docker-compose.yaml - DO NOT MODIFY!

# 3. Export NGC API key and login to NGC container registry
export NGC_CLI_API_KEY=<insert>
docker login nvcr.io --username='$oauthtoken' --password="${NGC_CLI_API_KEY}"
cd ~
cd protein-binder-design/deploy/
docker compose up

###########################################################################################################
# Alternative download of the models
###########################################################################################################

# docker run -it \
#     --runtime=nvidia \
#     --gpus='"device=0"' \
#     -p 8082:8000 \
#     -e NGC_CLI_API_KEY \
#     -v "$LOCAL_NIM_CACHE":/opt/nim/.cache \
#     nvcr.io/nim/ipd/rfdiffusion:2.2.0

# docker run -it \
#     --runtime=nvidia \
#     --gpus='"device=0"' \
#     -p 8083:8000 \
#     -e NGC_CLI_API_KEY \
#     -v "$LOCAL_NIM_CACHE":/home/nvs/.cache/nim \
#     nvcr.io/nim/ipd/proteinmpnn:latest

curl localhost:8082/v1/health/ready # RFdiffusion
curl localhost:8083/v1/health/ready # Protein MPNN
    

###########################################################################################################
# Clean up script format
###########################################################################################################

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

###########################################################################################################
# OPTION 1: No hotspots; generate i number of RFDiffusion backbones, and num_seq PorteinPMNN sequences per backbone
###########################################################################################################

# ------ Define protein ------
protein="1TNF" # 1TNF, apob, tnf
REPO_DIR="$HOME/protein-binder-design"
raw_pdb="${REPO_DIR}/input/${protein}.pdb"             # target protein

# ------ Input chain (region where the generated peptides will bind to) ------
#       in this case, it is the trimer's chain A (AA position 6 to 157)
chain="A" 
start_pos="6"
end_pos="157"

# ------ peptide design parameters ------
diffusion=50
temp=0.3
peptide_length="15-25" # set a range (i.e."15-25") or a value ("25")
i=10 # RFDiffusion iterations 
num_seq=4 #Protein_PMNN sequences per RFDiffusion structure (sequences in the same FASTA file)
export CUDA_VISIBLE_DEVICES=0

# ------ hotspots - yes or no? ------
unset hotspot_res # NO HOTSPOTS SET
# hotspot_res="${chain}${hotspot_res_prefix}" # SET A HOTSPOT

# Other variables 
contigs="${chain}${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
name="target_${chain}${start_pos}_${end_pos}"
params="${diffusion}diff_${temp}temp"

# Step 1: Build target structure PDB and extract target seq amino acid
target_pdb="${REPO_DIR}/input/${name}.pdb"
bash ${REPO_DIR}/scripts/get_target_pdb.sh "${raw_pdb}" "${target_pdb}" "${chain}" "${start_pos}" "${end_pos}"
if [ -f "$target_pdb" ]; then target_sequence=$(bash "${REPO_DIR}/scripts/get_target_seq.sh" "${target_pdb}"); fi

# Step 2: Run protein binder design
# for iteration in $(seq 1 $i); do
    # for num in $(seq 1 $num_seq); do

    python3.11 "${REPO_DIR}/scripts/3_protein_binder_design.py" \
        --root "${REPO_DIR}" --num_seq "${num_seq}" --diffusion "${diffusion}" --temp "${temp}" \
        --target_sequence "${target_sequence}" --contigs "${contigs}" --i "${i}" \
        --hotspot_res "${hotspot_res}" --target_pdb "${target_pdb}" --chain "${chain}"

    # Step 3: Generate merged binding alignment for peptide-target protein, and optimize alignment 
    python3.13 "${REPO_DIR}/scripts/4_merge_seq_to_backbone.py" "${REPO_DIR}" "${chain}" "${i}" "${num_seq}" "${name}" "${params}" --solvent
    
    # Step 4: Run PRODIGY
    bash "${REPO_DIR}/scripts/5_run_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${i}" "${num_seq}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}" "${input_file}" "${hotspot_res}"
    
    # done
# done

###########################################################################################################
# OPTION 2: Loop through each line within $input_file, defining hotspots and start/end positions
# Requires additionally:
#   -   Input file with each row being: Chain, hotspot residues, start/end positions

###########################################################################################################
input_file="${REPO_DIR}/input/target_file_${protein}_surface.txt"  # chain, hotspot residue, start/end pos 
head $input_file
# A 60 50 90
# A 156 50 157

sed -n '13p' "$input_file" | while IFS=$'\t' read -r chain hotspot_res_prefix start_pos end_pos; do
    
    # ------ Define variables WITHIN loop ------
    diffusion=50
    temp=0.3
    i=4
    num_seq=2
    peptide_length="15-25" # set a range (i.e."15-25") or a value ("25")

   # Variables (do not modify)
    hotspot_res="${chain}${hotspot_res_prefix}"
    contigs="${chain}${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
    name="target_${chain}${start_pos}_${end_pos}"
    params="${diffusion}diff_${temp}temp"
    
    echo "Processing chain=${chain}, hotspot_res=${hotspot_res}, start_pos=${start_pos}, end_pos=${end_pos}"
    echo "name=${name},    params=${params},    contigs=${contigs}"
    echo "=========================================================================="
    export CUDA_VISIBLE_DEVICES=0

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

# ---------------- Clean up NVIDIA workspace if necessary  ----------------
# # remove all files from ${REPO_DIR}/${name} directory that begin with "5_target_" and end with ".pdb"
# ls ./5_target_*.pdb
# rm ./5_target_*.pdb
# # rename all files that begin with "5_target_" and end with "_complex.pdb_old" to remove the "_old" suffix
# for file in 5_target_*_complex.pdb_old; do
#     mv "$file" "${file%_old}"
# done

###########################################################################################################
# OPTION 3: Loop through each line within $input_file, defining hotspots positions ONLY
# Requires additionally:
#   -   Input file with each row being: a list of hotspot residues
#   -   Must specify which lines (in this case, lines 1 - 21)
###########################################################################################################
# Define protein
protein="1TNF" # 1TNF, apob, tnf
REPO_DIR="$HOME/protein-binder-design"
raw_pdb="${REPO_DIR}/input/${protein}.pdb"             # target protein

input_file="${REPO_DIR}/input/target_file_${protein}_surface.txt"
head $input_file
# 'A6','A7','A8','A9','A34'
# 'A10','A11','A38','A39','A156'

sed -n '1,21p' "$input_file" | while IFS=$'\t' read -r hotspot_res_prefix; do
    # Define script input variables
    diffusion=50
    temp=0.3
    i=10 # RFDiffusions structures
    num_seq=4 #Protein_PMNN sequences per structure
    peptide_length="20-50" # set a range (i.e."15-25") or a value ("25")
    chain="A" 
    export CUDA_VISIBLE_DEVICES=1

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
    bash ${REPO_DIR}/scripts/get_target_pdb.sh "${raw_pdb}" "${target_pdb}" "${chain}" "${start_pos}" "${end_pos}"
    if [ -f "$target_pdb" ]; then target_sequence=$(bash "${REPO_DIR}/scripts/get_target_seq.sh" "${target_pdb}"); fi
    echo "Target sequence: ${target_sequence}"

    # # Step 2: Run the protein binder design script
    python3.11 "${REPO_DIR}/scripts/3_protein_binder_design.py" --root "${REPO_DIR}" \
    --num_seq "${num_seq}" --diffusion "${diffusion}" --temp "${temp}" --target_sequence "${target_sequence}" \
    --contigs "${contigs}" --i "${i}" --hotspot_res ${hotspot_res} --target_pdb "${target_pdb}" --chain "${chain}"

    # # Step 3: Generate merged binding alignment for peptide-target protein, and optimize alignment 
    python3.13 "${REPO_DIR}/scripts/4_merge_seq_to_backbone.py" "${REPO_DIR}" "${chain}" "${i}" "${num_seq}" "${name}" "${params}" --solvent
    bash "${REPO_DIR}/scripts/5_run_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${i}" "${num_seq}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}" "${input_file}" "${hotspot_res}"
done



#################################################################################
# OPTION 4: Loop through each line within $input_file, defining hotspots positions ONLY
# CODE WORKS!
# Creates summary file .txt
#################################################################################

# Initialize the summary file with headers
summary_file="${REPO_DIR}/summary_${protein}.txt"
echo -e "Target\titeration\tnum_seq\tbinding_affinity\tdiss_constant" > "$summary_file"

# Read input file and process for each chain, hotspot, start, and end position
sed -n '1,21p' "$input_file" | while IFS=$'\t' read -r hotspot_res_prefix; do
    diffusion=50
    temp=0.3
    i=10 # RFDiffusions structures
    num_seq=4 #Protein_PMNN sequences per structure
    peptide_length="20-50" # set a range (i.e."15-25") or a value ("25")
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
                echo "Warning: '$aligned_pdb' not found."
            fi
            echo ""
            unset binding_affinity
            unset diss_constant
        done
    done
done
echo "Summary: $summary_file"


####### If hotspots are a combination of different residues ##########

summary_file="${REPO_DIR}/summary_${protein}.txt"
echo -e "Target\titeration\tnum_seq\tbinding_affinity\tdiss_constant" > "$summary_file"

sed -n '1,21p' "$input_file" | while IFS=$'\t' read -r chain hotspot_res_prefix start_pos end_pos; do
    diffusion=50
    temp=0.3
    i=10 # RFDiffusions structures
    num_seq=4 #Protein_PMNN sequences per structure
    peptide_length="20-50" # set a range (i.e."15-25") or a value ("25")
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
                echo "Warning: '$aligned_pdb' not found."
            fi
            echo ""
#             unset binding_affinity
#             unset diss_constant
        done
    done
done
echo "Summary file: $summary_file"


###########################################################################
# OPTION 5: Loop through each line of $input_file
###########################################################################

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

    # Step 3: Generate merged binding alignment for peptide and target protein, and then optimize alignment 
    python3.13 ${REPO_DIR}/scripts/4_merge_seq_to_backbone.py "${REPO_DIR}" A ${i} ${num_seq} ${name} ${params} --solvent

    # Step 5 alternative: PRODIGY
    bash "${REPO_DIR}/scripts/calc_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${i}" "${num_seq}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}" "${input_file}"

done < "$input_file"
