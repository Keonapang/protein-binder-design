#!/usr/bin/env bash
# Command-line tool for PRODIGY

# chmod +x /home/ubuntu/protein-binder-design/scripts/get_target_pdb.sh
# get_target_pdb ${input_file} ${output_file} ${chain} ${start_pos} ${end_pos}

# start_pos=60
# end_pos=90
# chain="A"
# input_file="/home/ubuntu/protein-binder-design/input/pdb2e7a.pdb"
# output_file="/home/ubuntu/protein-binder-design/input/target_${chain}${start_pos}_${end_pos}.pdb"

sudo apt install python3.11 
sudo apt update

# Install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh # installed to /home/shadeform/miniconda3

# Assign path to conda
echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc # conda --version

# Allow override permissions
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# **Set up `openmm` and `pdbfier` to align backbones (requires python 3.13):**
cd ~
git clone https://github.com/openmm/pdbfixer
cd pdbfixer
python setup.py install
python3.13 -m pip install numpy prodigy-prot torch Bio biopython pdb-tools
python3.13 -m pip prodigy-prot

# use Conda to install pdbfixer and openmm (writes to python3.13)
conda create -n pdbfixer_env python=3.13 -y
conda init
conda activate pdbfixer_env
conda install -c conda-forge pdbfixer
conda install -c conda-forge openmm
python3.13 -m pip install numpy prodigy-prot torch Bio biopython pdb-tools
python3.13 -m pip prodigy-prot
########################################################################################################### 

# Optional 
sed -i 's/\r$//' "${REPO_DIR}/scripts/get_target_pdb.sh" # optional (to remove any hidden spaces from Windows)
sed -i 's/\r$//' "${REPO_DIR}/scripts/5_run_prodigy.sh" 
sed -i 's/^[ \t]*//;s/[ \t]*$//' "${REPO_DIR}/scripts/5_run_prodigy.sh"
sed -i 's/\r$//' "${REPO_DIR}/input/target_hotspots.txt" 
chmod +x "${REPO_DIR}/scripts/get_target_pdb.sh"

# Ensure any space delimited regions are all converted to tab-delimited 
awk '{$1=$1; gsub(" ", "\t"); print}' "$input_file" > "$input_file.tmp" && mv "$input_file.tmp" "$input_file"
sed -i 's/\r$//' "$input_file"
###########################################################################################################

# start_pos=91 
# end_pos=130
# name="target_${chain}${start_pos}_${end_pos}"
# params="${diffusion}diff_${temp}temp"
# aligned_pdb=${REPO_DIR}/${name}/5_${name}_${params}_i${i}_${num_seq}_complex.pdb
# aligned_pdb=${REPO_DIR}/target_A91_130/5_target_A91_130_50diff_0.3temp_i1_1_complex.pdb
wc -l $aligned_pdb

# Define protein
protein="1TNF" # 1TNF, apob, tnf

# Define repo directory  
REPO_DIR="/home/shadeform/protein-binder-design"

# Two input files
raw_pdb="${REPO_DIR}/input/${protein}.pdb"             # target protein
input_file="${REPO_DIR}/input/target_file_${protein}.txt"  # chain, hotspot residue, start/end pos 

# Define script input variables
diffusion=50
temp=0.4
i=2
num_seq=2
peptide_length="15-25" # set a range (i.e."15-25") or a value ("25")
chain="C"

# check before running!
wc -l $raw_pdb
wc -l $input_file

# Clean up raw PDB format
chmod +x ${REPO_DIR}/scripts/fix_pdb_format.sh
bash "${REPO_DIR}/scripts/fix_pdb_format.sh" "$raw_pdb"

# Convert residue position file to tab-delimited
awk '{$1=$1; gsub(" ", "\t"); print}' "$input_file" > "$input_file.tmp" && mv "$input_file.tmp" "$input_file"
sed -i 's/\r$//' "$input_file"
head $input_file

export NGC_CLI_API_KEY=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
docker login nvcr.io --username='$oauthtoken' --password="${NGC_CLI_API_KEY}"

sed -n '6p' "$input_file" | while IFS=$'\t' read -r chain hotspot_res_prefix start_pos end_pos; do
    # Variables (do not modify)
    hotspot_res="${chain}${hotspot_res_prefix}"
    contigs="${chain}${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
    name="target_${chain}${start_pos}_${end_pos}"
    params="${diffusion}diff_${temp}temp"
    echo "Processing chain=${chain}, hotspot_res=${hotspot_res}, start_pos=${start_pos}, end_pos=${end_pos}"
    echo "name=${name},    params=${params},    contigs=${contigs}"
    echo ""
    echo "=========================================================================="
    export CUDA_VISIBLE_DEVICES=1
    # unset hotspot_res # NO HOTSPOTS SET

    # Step 1: Build target structure PDB and extract target seq amino acid
    target_pdb="${REPO_DIR}/input/${name}.pdb" # Output PDB path
    bash ${REPO_DIR}/scripts/get_target_pdb.sh "${raw_pdb}" "${target_pdb}" "${chain}" "${start_pos}" "${end_pos}"
    if [ -f "$target_pdb" ]; then target_sequence=$(bash "${REPO_DIR}/scripts/get_target_seq.sh" "${target_pdb}"); fi
    head $target_pdb
    echo $target_sequence
    echo ""
    # # # Step 2: Run the protein binder design script
    # python3.11 "${REPO_DIR}/scripts/3_protein_binder_design.py" --root "${REPO_DIR}" \
    # --num_seq "${num_seq}" --diffusion "${diffusion}" --temp "${temp}" --target_sequence "${target_sequence}" \
    # --contigs "${contigs}" --i "${i}" --hotspot_res "${hotspot_res}" --target_pdb "${target_pdb}" --chain "${chain}"

    # # Step 3: Generate merged binding alignment for peptide and target protein, and then optimize alignment 
    # python3.13 ${REPO_DIR}/scripts/4_merge_seq_to_backbone.py "${REPO_DIR}" ${chain} ${i} ${num_seq} ${name} ${params} --solvent
    bash "${REPO_DIR}/scripts/5_run_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${i}" "${num_seq}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}" "${input_file}"
done


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

