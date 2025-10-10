#!/usr/bin/env bash
# Command-line tool for PRODIGY

# chmod +x /home/ubuntu/protein-binder-design/scripts/get_target_pdb.sh
# get_target_pdb ${input_file} ${output_file} ${chain} ${start_pos} ${end_pos}

# start_pos=60
# end_pos=90
# chain="A"
# input_file="/home/ubuntu/protein-binder-design/input/pdb2e7a.pdb"
# output_file="/home/ubuntu/protein-binder-design/input/target_${chain}${start_pos}_${end_pos}.pdb"


# Optional 
sed -i 's/\r$//' "${REPO_DIR}/scripts/get_target_pdb.sh" # optional (to remove any hidden spaces from Windows)
sed -i 's/\r$//' "${REPO_DIR}/scripts/calc_prodigy.sh" 
sed -i 's/\r$//' "${REPO_DIR}/input/target_hotspots.txt" 
sed -i 's/^[ \t]*//;s/[ \t]*$//' "${REPO_DIR}/scripts/calc_prodigy.sh"
chmod +x "${REPO_DIR}/scripts/get_target_pdb.sh"

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

# use Conda to install pdbfixer and openmm (writes to python3.13)
conda create -n pdbfixer_env python=3.13 -y
conda init
conda activate pdbfixer_env
conda install -c conda-forge pdbfixer
conda install -c conda-forge openmm
python3.13 -m pip install numpy prodigy-prot torch Bio biopython pdb-tools

########################################################################################################### 

# Define protein
protein="apob"

# Define repo directory  
REPO_DIR="/home/shadeform/protein-binder-design"

# Two input files
raw_pdb="${REPO_DIR}/input/${protein}.pdb"             # target protein
input_file="${REPO_DIR}/input/target_hotspots_${protein}.txt"  # chain, hotspot residue, start/end pos 

# Define script input variables
diffusion=50
temp=0.3 
i=5
num_seq=2
peptide_length="15-25" # set a range (i.e."15-25") or a value ("25")
chain="A"

# check before running!
wc -l $raw_pdb
awk '$1 == "ATOM" && $5 == "A"' /home/shadeform/protein-binder-design/input/apob.pdb | head

# Convert space-delimited file to tab-delimited
awk '{$1=$1; gsub(" ", "\t"); print}' "$input_file" > "$input_file.tmp" && mv "$input_file.tmp" "$input_file"
sed -i 's/\r$//' "$input_file"

while read -r line; do
    chain=$(echo "$line" | cut -f1)              # First column
    hotspot_res_prefix=$(echo "$line" | cut -f2) # Second column
    start_pos=$(echo "$line" | cut -f3)          # Third column
    end_pos=$(echo "$line" | cut -f4)            # Fourth column
    hotspot_res="${chain}${hotspot_res_prefix}"
    contigs="A${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
    name="target_${chain}${start_pos}_${end_pos}"
    params="${diffusion}diff_${temp}temp"
    echo "Processing chain=${chain}, hotspot_res=${hotspot_res}, start_pos=${start_pos}, end_pos=${end_pos}"
    echo "name=${name},    params=${params},    contigs=${contigs}"
    echo ""

    # Step 1: Build target structure PDB and extract target seq amino acid
    # source "${REPO_DIR}/scripts/get_target_pdb.sh" # Load the get_target_pdb function
    # get_target_pdb "${raw_pdb}" "${target_pdb}" "${chain}" "${start_pos}" "${end_pos}"
    target_pdb="/home/shadeform/protein-binder-design/input/${name}.pdb" # output
    # bash ${REPO_DIR}/scripts/get_target_pdb.sh "${raw_pdb}" "${target_pdb}" "${chain}" "${start_pos}" "${end_pos}"
    echo $target_pdb
    echo ""
done < "$input_file"


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
    target_pdb="/home/shadeform/protein-binder-design/input/${name}.pdb" # output
    bash ${REPO_DIR}/scripts/get_target_pdb.sh "${raw_pdb}" "${target_pdb}" "${chain}" "${start_pos}" "${end_pos}"
    if [ -f "$target_pdb" ]; then target_sequence=$(bash "${REPO_DIR}/scripts/get_target_seq.sh" "${target_pdb}"); fi

    # Step 2: Run the protein binder design script
    python3.11 "${REPO_DIR}/scripts/3_protein_binder_design.py" --root "${REPO_DIR}" \
    --num_seq "${num_seq}" --diffusion "${diffusion}" --temp "${temp}" --target_sequence "${target_sequence}" \
    --contigs "${contigs}" --i "${i}" --hotspot_res "${hotspot_res}" --target_pdb "${target_pdb}"
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
    # time: ~2mins per structure
    python3.13 ${REPO_DIR}/scripts/4_merge_seq_to_backbone.py "${REPO_DIR}" A ${i} ${num_seq} ${name} ${params} --solvent

    # Step 5 alternative: PRODIGY
    bash "${REPO_DIR}/scripts/calc_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${i}" "${num_seq}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}" "${input_file}"

done < "$input_file"


# Process Lines 3-4
sed -n '3,4p' "$input_file" | while IFS=$'\t' read -r chain hotspot_res_prefix start_pos end_pos; do
    # Variables (do not modify)
    hotspot_res="${chain}${hotspot_res_prefix}"
    contigs="A${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
    name="target_${chain}${start_pos}_${end_pos}"
    params="${diffusion}diff_${temp}temp"
    echo "Processing chain=${chain}, hotspot_res=${hotspot_res}, start_pos=${start_pos}, end_pos=${end_pos}"
    echo "name=${name},    params=${params},    contigs=${contigs}"
    echo ""

    # Step 3: Generate merged binding alignment for peptide and target protein, and then optimize alignment 
    # time: ~2mins per structure
    python3.13 ${REPO_DIR}/scripts/4_merge_seq_to_backbone.py "${REPO_DIR}" A ${i} ${num_seq} ${name} ${params} --solvent
    chmod +x "${REPO_DIR}/scripts/6_run_mmpbsa.sh"
    bash "${REPO_DIR}/scripts/calc_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${i}" "${num_seq}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}" "${input_file}"
done


sed -n '5,6p' "$input_file" | while IFS=$'\t' read -r chain hotspot_res_prefix start_pos end_pos; do
    # Variables (do not modify)
    hotspot_res="${chain}${hotspot_res_prefix}"
    contigs="A${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
    name="target_${chain}${start_pos}_${end_pos}"
    params="${diffusion}diff_${temp}temp"
    echo "Processing chain=${chain}, hotspot_res=${hotspot_res}, start_pos=${start_pos}, end_pos=${end_pos}"
    echo "name=${name},    params=${params},    contigs=${contigs}"
    echo ""

    # Step 3: Generate merged binding alignment for peptide and target protein, and then optimize alignment 
    # time: ~2mins per structure
    python3.13 ${REPO_DIR}/scripts/4_merge_seq_to_backbone.py "${REPO_DIR}" A ${i} ${num_seq} ${name} ${params} --solvent
    chmod +x "${REPO_DIR}/scripts/6_run_mmpbsa.sh"
    bash "${REPO_DIR}/scripts/calc_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${i}" "${num_seq}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}" "${input_file}"
done