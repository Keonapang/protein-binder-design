# Adaption of NVIDIA BioNeMo: Protein Binder Design

This workflow is designed for _in silico_ protein binder design by generating binder sequences and predicted structures for the binder and target.

Two input files are required: (1) a `.pdb` structure of the target protein, and (2) a `.txt` containing selected target binding hotspots and start/stop positions.

**RFDiffusion** will first generate the protein backbones for binder design. Subsequently, **ProteinMPNN** predicts the amino acid sequence of this peptide binder, and then **AlphaFold2** outputs the peptide PDB structure. Finally, the aligned binder-target complex is validated (step #4) on [PRODIGY](https://rascar.science.uu.nl/prodigy/). Optionally, you may visualize the complex using Rosetta [FlexPepDoc](https://r2.graylab.jhu.edu/auth/login?next=%2Fapps%2Fsubmit%2Fflexpepdock) or PyMOL.

### System Requirements

- at least **1500 GB** of fast NVMe SSD space
- at least **24 CPU** cores
- at least **128 GB RAM**

### Hardware requirements

3 or more NVIDIA L40s, A100, or H100 GPUs
- **RFdiffusion** runs on 1 x GPU, ≥12 GiB GPU memory, 15GB free SSD drive space
- **ProteinMPNN** runs on 1 x GPU, ≥3 GiB GPU memory, 10GB free SSD drive space
- **AlphaFold2** runs on 1 x GPU, ≥32 GiB GPU memory, 1300GB free SSD drive space

### Software Pre-requisites

- Python 3.11+


## 1. Ensure you have prepared two input files

This script takes in two arguments:

**1. Raw `.PDB` file** of the target protein, most likely downloaded through a protein database. The script will ignore all irrelevant rows in the PDB file except the 'ATOM' rows. *Note: the file doesn't have to begin with amino acid residue 1*. See [/input/pdb2e7a.pdb](https://github.com/Keonapang/protein-binder-design/blob/main/input/pdb2e7a.pdb).

```bash
HEADER    CYTOKINE                                09-JAN-07   2E7A              
TITLE     TNF RECEPTOR SUBTYPE ONE-SELECTIVE TNF MUTANT WITH ANTAGONISTIC       
JRNL       DOI    10.1074/JBC.M707933200                                       
REMARK   2 RESOLUTION.    1.80 ANGSTROMS.                                       
REMARK   3     
ATOM      1  N   PRO A   8      18.727  24.301  31.792  1.00 56.82           N  
ATOM      2  CA  PRO A   8      17.276  24.324  31.476  1.00 57.03           C  
```

**2. Space-delimited `.txt` file (no header)** containing four columns that store: chain identifier (e.g. A), hotspot residue (e.g. 80), start and end residue position (e.g. 60 and 90). See [/input/target_hotspots.txt](https://github.com/Keonapang/protein-binder-design/blob/main/input/target_hotspots.txt)

```bash
    A 80 60 90
    A 90 80 110
    A 103 83 113
    A 113 103 133
```

## 2. Launch a NVIDIA cloud virtual machine (VM)

1. Sign into [brev.nvidia](https://brev.nvidia.com/).
2. Click on this [Launchable](https://brev.nvidia.com/launchable/deploy?launchableID=env-32hOwzKwjmWbUwZHWYtSoW9vrxO). If you would like to create your own, click on **Launchables** in the menu bar, then `Create Launchables` and follow these settings:
    - Select "I have codes in a git repository" and enter the URL of this repo
    - Select 'VM-mode'
    - Click "Next" until you finally reach **select compute**. We recommended **A100 (80GiB) 2 GPUs x 24 CPUs | 240GiB  (80GiB GPU memory) ($3.96/hr)**

![NVIDIA VM settings](docs/NVIDIA_VM.png)

3. Click **"Deploy Launchable"** and **"Go to Instance Page"**. Wait ~10 minutes for VM to start
4. Enter the VM, then drag and drop the two input files from **step (2)** into your workspace under `protein-binder-design/input/`
5. Start a new terminal session by clicking the top-right **"+"** button.

6. An **NGC Personal API Key** is required to download and run any NVIDIA NIMs. If this is your first time, start by creating an account on [NGC](https://catalog.ngc.nvidia.com/). Then [generate the key](https://org.ngc.nvidia.com/setup/api-key) and note it down somewhere secure for future use.

```bash
    export NGC_CLI_API_KEY=<enter-key> 
    # export NGC_CLI_API_KEY=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
    docker login nvcr.io --username='$oauthtoken' --password="${NGC_CLI_API_KEY}"
```

7. Build inference models via **docker** container. Wait 20-30mins for the docker to build.

```bash
    # Install Dependencies
    sudo apt-get update # updated nvidia toolkit
    sudo apt-get install -y docker-compose # docker compose version 2+
    sudo apt install python3.11

    #  NIM cache allows you to download models and store previously-downloaded models on your local/server disk
    mkdir -p ~/.cache/nim
    chmod -R 777 ~/.cache/nim    
    export HOST_NIM_CACHE=~/.cache/nim

    # From the root of the cloned protein-binder-design repository:
    cd ~
    cd protein-binder-design/deploy/
    docker compose up # runs /deploy/docker-compose.yaml
```

8. Open a new terminal, leaving the current terminal open with the launched service.

9. In the new terminal, view running containers with the following codes. Wait until the health check returns `{"status":"ready"}` before proceeding.

```bash
    # 1. Check storage space
    df -h 

    # 2. Check which dockers are currently active/on stand-by
    docker container ls
    docker ps 

    # 3. Health check  # shoud be {"status":"ready"}
    # curl localhost:8081/v1/health/ready # AlphaFold2
    curl localhost:8082/v1/health/ready # RFdiffusion
    curl localhost:8083/v1/health/ready # Protein MPNN
    
    # Example: check download log
    docker logs -f protein-binder-design-alphafold-1
```

## 2. Install conda

Download and activate `conda`:

```bash
# Install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh # installed to /home/shadeform/miniconda3

# Assign path to conda
echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc # conda --version

# Allow override permissions
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
```

## 3. Install pre-requisites

```bash
# Install python, if container doesn't come built with it
sudo apt install python3.11 
pip install prodigy-prot torch Bio
# Install essentials
sudo apt update
```

**Set up `openmm` and `pdbfier` to align backbones (requires python 3.13):**

```bash
cd ~
git clone https://github.com/openmm/pdbfixer
cd pdbfixer
python setup.py install

# use Conda to install pdbfixer and openmm (writes to python3.13)
conda create -n pdbfixer_env python=3.13 -y
conda activate pdbfixer_env
conda install -c conda-forge pdbfixer
conda install -c conda-forge openmm
python3.13 -m pip install numpy prodigy-prot torch Bio biopython pdb-tools

```

**Set up environment to run `gromacs` for free energy calculation**

**Build via conda:**
The Conda package typically defaults to OpenC for GPU compatibility, not CUDA.
OpenCL is not as well-optimized for NVIDIA GPUs.
Use if the system uses an AMD or Intel GPU, or you don’t need GPU acceleration.

```bash
# Download of GROMACS via Conda (5-10mins), however this doesn't enable GPU support
cd ~
conda create -n mmpbsa python=3.11 -y
conda activate mmpbsa
conda install -c conda-forge libgcc-ng mpi mpi4py compilers gcc_linux-64 gxx_linux-64
conda install -c conda-forge libmpich-dev libopenmpi-dev libopenmpi-dev openmpi-bin -y 
conda install -c conda-forge ocl-icd-system
conda install -c conda-forge clinfo
conda install -c conda-forge intel-compute-runtime
conda install -c conda-forge pocl oclgrind
conda install -c conda-forge gromacs -y  # not ideal!!
```

**Build GROMACS from source with cmake:**
Download and build the latest stable or development version of GROMACS directly from the source.
Can compile GROMACS with options that are optimized for your specific CPU (e.g., AVX2, AVX512) and GPU (e.g., CUDA).

```bash
#######################################################
# Ensure GROMACS is compiled with CUDA support and 
# links correctly to the NVIDIA drivers and libraries. 
#######################################################
conda create -n gromacs_gpu python=3.11 -y
conda activate gromacs_gpu
conda install -c conda-forge cudatoolkit=12.8

# Install latest cmake (importnat)
sudo apt update
sudo apt install snapd -y
sudo snap install cmake --classic
export PATH=/snap/bin:$PATH
cmake --version # cmake version 4.1.2
echo 'export PATH=/snap/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

conda install -c conda-forge cudatoolkit # version 11.8 installed
export PATH=$CONDA_PREFIX/bin:$PATH
export CUDACXX=$CONDA_PREFIX/bin/nvcc
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

conda install -c conda-forge cmake gcc_linux-64 gxx_linux-64 fftw
conda install -c conda-forge openmpi
sudo apt install openmpi-bin openmpi-common libopenmpi-dev
mpicc --version
mpicxx --version

export PATH=$CONDA_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export PATH=/usr/bin:$PATH

# Download and Extract GROMACS from source
cd ~
wget ftp://ftp.gromacs.org/gromacs/gromacs-2025.3.tar.gz
tar xfz gromacs-2025.3.tar.gz
cd 

# Build GROMACS with CUDA Support
mkdir build
cd build
sudo apt install libopenblas-dev liblapack-dev
sudo apt install libblas-dev liblapack-dev
sudo apt install nvidia-cuda-toolkit
sudo apt install openmpi-bin libopenmpi-dev
export CUDACXX=/usr/local/cuda/bin/nvcc
cmake .. \
    -DGMX_BUILD_OWN_FFTW=OFF \
    -DGMX_MPI=ON \
    -DGMX_GPU=CUDA \
    -DCUDA_TOOLKIT_ROOT_DIR=$CONDA_PREFIX \
    -DGMX_CUDA_TARGET_COMPUTE="80" \
    -DCMAKE_INSTALL_PREFIX=$HOME/gromacs \
    -DCMAKE_C_COMPILER=/usr/bin/mpicc \
    -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx \
    -DCMAKE_CUDA_COMPILER=$CUDACXX \
    -DFFTWF_LIBRARY=$HOME/fftw/lib/libfftw3f.so \
    -DFFTWF_INCLUDE_DIR=$HOME/fftw/include


# mpicxx --version
# gcc (Ubuntu 11.4.0-1ubuntu1~22.04.2) 11.4.0
# Copyright (C) 2021 Free Software Foundation, Inc.
# This is free software; see the source for copying conditions.  There is NO
# warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# g++ (Ubuntu 11.4.0-1ubuntu1~22.04.2) 11.4.0
# Copyright (C) 2021 Free Software Foundation, Inc.
# This is free software; see the source for copying conditions.  There is NO
# warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


####################################################################
# Download GROMACS and build manually (https://manual.gromacs.org/)
# Everything is installed system-wide under /usr/local
# requires cmake 4.1 and also CUDA 12.1
####################################################################
sudo apt update
sudo apt install -y build-essential cmake gcc g++ libfftw3-dev libopenmpi-dev openmpi-bin

# # If CUDA is not pre-built, install it (v12.1 or higher)
# nvcc --version 

# wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-ubuntu2204.pin
# sudo mv cuda-ubuntu2204.pin /etc/apt/preferences.d/cuda-repository-pin-600
# wget https://developer.download.nvidia.com/compute/cuda/12.1.0/local_installers/cuda-repo-ubuntu2204-12-1-local_12.1.0-530.30.02-1_amd64.deb
# sudo dpkg -i cuda-repo-ubuntu2204-12-1-local_12.1.0-530.30.02-1_amd64.deb
# sudo cp /var/cuda-repo-ubuntu2204-12-1-local/cuda-*-keyring.gpg /usr/share/keyrings/
# sudo apt-get update
# sudo apt-get -y install cuda
# sudo apt update
# sudo apt install -y cuda
# echo 'export PATH=/usr/local/cuda/bin:$PATH' >> ~/.bashrc
# echo 'export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
# source ~/.bashrc

# Install latest cmake (importnat)
sudo apt update
sudo apt install snapd -y
sudo snap install cmake --classic
export PATH=/snap/bin:$PATH
cmake --version

# Download and Extract GROMACS from source
cd ~
wget ftp://ftp.gromacs.org/gromacs/gromacs-2025.3.tar.gz
tar xfz gromacs-2025.3.tar.gz
cd gromacs-2025.3

nano ~/.bashrc
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
export PATH=/snap/bin:$PATH
source ~/.bashrc

conda install -c conda-forge cudatoolkit=11.8y
conda install -c conda-forge openmpi
sudo apt install openmpi-bin openmpi-common libopenmpi-dev
mpicc --version
mpicxx --version

# Build GROMACS with CUDA Support
mkdir build
cd build
sudo apt install libopenblas-dev liblapack-dev
cmake .. \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DGMX_GPU=CUDA \
    -DGMX_OPENMP=ON \
    -DREGRESSIONTEST_DOWNLOAD=ON \
    -DCUDA_ARCHITECTURES=ALL
make -j$(nproc)
make check
sudo make install
source /usr/local/gromacs/bin/GMXRC
gmx --version
gmx mdrun -deviceinfo # check compatability 

# Download nvpl
wget https://developer.download.nvidia.com/compute/nvpl/25.5/local_installers/nvpl-local-repo-ubuntu2404-25.5_1.0-1_arm64.deb
sudo dpkg -i nvpl-local-repo-ubuntu2404-25.5_1.0-1_arm64.deb
sudo cp /var/nvpl-local-repo-ubuntu2404-25.5/nvpl-*-keyring.gpg /usr/share/keyrings/
sudo apt-get update
sudo apt-get -y install nvpl

# Update GCC and G++ to match versions
sudo apt install g++-12
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-12 12
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-12 12

python -m pip install gmx-MMPBSA # (doesn't work for python3.11)

# check GPU acceleration
nvidia-smi
gmx mdrun -deviceinfo
gmx mdrun -gpu_id 0 -nsteps 1000 -v
```

## 4. Run clean up and installation

```bash
    # clean up scripts before running
    sed -i 's/\r$//' "${REPO_DIR}/scripts/get_target_pdb.sh" # optional (to remove any hidden spaces from Windows)
    sed -i 's/\r$//' "${REPO_DIR}/scripts/calc_prodigy.sh" 
    sed -i 's/\r$//' "${REPO_DIR}/input/target_hotspots.txt" 
    sed -i 's/^[ \t]*//;s/[ \t]*$//' "${REPO_DIR}/scripts/calc_prodigy.sh"
```

Export API key again

```bash
    # export NGC_CLI_API_KEY=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
    export NGC_CLI_API_KEY=<enter-key>
```

## 5. Run script

The command `get_target_pdb` checks that `start_pos` and `end_pos` are valid given the input PDB file. this script runs 3 prediction models sequentially, followed by aligning the designed peptide PDB to the original target PDB to create a combined PDB using BioPython's `Superimposer` module. Lastly, it calculates the **dissociation constant (Kd)** using [PRODIGY](https://github.com/haddocking/prodigy):

- **RFDiffusion**: takes in target protein PDB, outputs a peptide binder PDB. Note: every output is a glycine, see [Github](https://github.com/RosettaCommons/RFdiffusion?tab=readme-ov-file#understanding-the-output-files) for more info.
    - `contigs`: range of amino acid positions, and expected length of peptide (i.e."A1-30/0 15-25")
    - `diffusion`: number of diffusion_steps (default: 50)
    - `i`: number of RFDiffusion iterations for each target sequence

- **ProteinPMNN**: takes in the peptide PDB, and outputs a .fasta file storing a list of their predicted amino acid sequences
    - `num_seq`: how many seqs to generate for a given structure from RFDiffusion
    - `hotspot_res`: A20 # array (i.e. "A20")
    - `temp`: sample temperature controls the diversity of designed peptides. Higher values will lead to more diversity (range:0-1)

**Define variables**

```bash
# Define protein
protein="apob"

# Define repo directory  
REPO_DIR="/home/shadeform/protein-binder-design"

# Two input files
raw_pdb="${REPO_DIR}/input/${protein}.pdb"                     # target protein
input_file="${REPO_DIR}/input/target_hotspots_${protein}.txt"  # chain, hotspot residue, start/end pos 

```

```bash
# Define script input variables
chain="A"
diffusion=50
temp=0.3 
i=5
num_seq=2
peptide_length="15-25" # set a range (i.e."15-25") or a value ("25")

```

**Clean up the format of raw PDB file**
Just in case the chain ID and residue number are not properly spaced (larger protein files are more susceptible)


```bash
${REPO_DIR}/scripts/fix_pdb_format.sh $raw_pdb
```

**Clean up the format of the input file**


```bash
# Ensure any space delimited regions are all converted to tab-delimited 
awk '{$1=$1; gsub(" ", "\t"); print}' "$input_file" > "$input_file.tmp" && mv "$input_file.tmp" "$input_file"
sed -i 's/\r$//' "$input_file"

```

**Run through each line of the input file**

```bash

# Loop through each line of $input_file
while IFS=$'\t' read -r chain hotspot_res_prefix start_pos end_pos; do # space-delimited

    # Variables (do not modify)
    hotspot_res="${chain}${hotspot_res_prefix}"
    contigs="A${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
    name="target_${chain}${start_pos}_${end_pos}"
    params="${diffusion}diff_${temp}temp"
    echo "Processing chain=${chain}, hotspot_res=${hotspot_res}, start_pos=${start_pos}, end_pos=${end_pos}"
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

while IFS=$'\t' read -r chain hotspot_res_prefix start_pos end_pos; do # space-delimited

    # Variables (do not modify)
    hotspot_res="${chain}${hotspot_res_prefix}"
    contigs="A${start_pos}-${end_pos}/0 ${peptide_length}" # e.g. "A60-90/0 15-25"
    name="target_${chain}${start_pos}_${end_pos}"
    params="${diffusion}diff_${temp}temp"

    # Step 3: Generate merged binding alignment for peptide and target protein, and then optimize alignment 
    # time: ~2-5mins per structure
    python3.13 ${REPO_DIR}/scripts/4_merge_seq_to_backbone.py "${REPO_DIR}" A ${i} ${num_seq} ${name} ${params} --solvent

    # Step 4: Main binding free energy calculation  
    bash "${REPO_DIR}/scripts/5_run_mmpbsa.sh" ${REPO_DIR} ${name} ${params} ${i} ${num_seq}

    # Step 5 alternative: PRODIGY
    bash "${REPO_DIR}/scripts/5_run_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${i}" "${num_seq}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}" "${input_file}"

done < "$input_file"

```

### Parallel processing option (ensure its on CUDA GPU)

```bash
# Get the total number of lines in the input file
total_lines=$(wc -l < "$input_file")

sed -n "1,${total_lines}p" "$input_file" | while IFS=$'\t' read -r chain hotspot_res_prefix start_pos end_pos; do
    {
        # Variables (do not modify)
        hotspot_res="${chain}${hotspot_res_prefix}"
        contigs="A${start_pos}-${end_pos}/0 ${peptide_length}"
        name="target_${chain}${start_pos}_${end_pos}"
        params="${diffusion}diff_${temp}temp"

        # Assign GPUs dynamically (e.g., alternate between GPU 0 and GPU 1)
        line_index=$((++index % 2))  # Alternate between 0 and 1
        export CUDA_VISIBLE_DEVICES=$line_index

        echo "Running on GPU ${CUDA_VISIBLE_DEVICES}"
        echo "Processing chain=${chain}, hotspot_res=${hotspot_res}, start_pos=${start_pos}, end_pos=${end_pos}"
        echo "name=${name},    params=${params},    contigs=${contigs}"
        echo ""

        # Step 3: Generate merged binding alignment for peptide and target protein, and then optimize alignment 
        python3.13 "${REPO_DIR}/scripts/4_merge_seq_to_backbone.py" "${REPO_DIR}" A ${i} ${num_seq} ${name} ${params} --solvent
        bash "${REPO_DIR}/scripts/5_run_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${i}" "${num_seq}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}" "${input_file}"
    } &
done
# Wait for all background jobs to finish
wait
```

### Workflow outputs

Example of **final aligned binder-target PDB**:

```bash
# Final PDB complex of target protein and designed peptide binder
# Date: 2025-10-08
# Target sequence: SSQLSRALLSLLLAL
# RFDiffusion candidate 1, ProteinPMNN predicted sequence 1
#
# ========= Input parameters ========= 
# target_pdb= ./input/pdb2e7a.pdb
# input coordinates file: ./input/target_hotspots.txt
# input line num=1
# chain=A
# start_pos=60
# end_pos=90
# diffusion=50
# temp=30
# iteration=5 total RFDiffusion candidates
# num_seq=2 total ProteinMPNN sequences
# 
# ========= ProteinPMNN predicted peptide =========
# >T=0.3, sample=1, score=1.9195, global_score=2.2284, seq_recovery=0.0000
# SSQLSRALLSLLLAL
# 
# ========= PRODIGY results ========= 
# Binding affinity (kcal.mol-1): -7.6 
# Dissociation constant (Kb) at 25.0˚C: 2.8e-06
# No. of intermolecular contacts: 11
# No. of charged-charged contacts: 0.0
# No. of charged-polar contacts: 0.0
# No. of charged-apolar contacts: 3.0
# No. of polar-polar contacts: 0.0
# No. of apolar-polar contacts: 2.0
# No. of apolar-apolar contacts: 6.0
# Percentage of apolar NIS residues: 45.65
# Percentage of charged NIS residues: 4.35
#
# ATOM      1  N   PRO A   8      18.727  24.301  31.792  1.00 56.82           N  
# ATOM      2  CA  PRO A   8      17.276  24.324  31.476  1.00 57.03           C  
# ATOM      3  C   PRO A   8      16.970  25.033  30.160  1.00 59.51           C  
# ATOM      4  O   PRO A   8      17.158  26.246  30.040  1.00 60.32           O  

```

Example of **ProteinPMNN** .fasta output if `num_seq=2`:

```bash
    # binder_target_pairs
    [
        ['RIAELLAQLLKELLE','SQVLFSGQGCPSTHVLLTHTISRISTTHNQP'], # binder design 1
        ['AIEEALARLLLEQLL', 'SQVLFSGQGCPSTHVLLTHTISRISTTHNQP'] # binder design 2
    ]
```


Example of **PRODIGY** output in the terminal:

![PRODIGY](docs/prodigy.png)

> [!NOTE]
> Once you're done running the scrips above, you may shutdown the NVIDIA VM so it doesn't keep charging money.

## 7. Visualization and validation of binder-target (local)

Before you proceed, ensure that all your predicted peptide binder structures are now stored in local directory `$DIR_WORK`.

### (a) PRODIGY Gibbs Free Energy

The [PRODIGY web server](https://rascar.science.uu.nl/prodigy/) tool can provide a prediction of binder-target binding affinity. The outputs include the Gibbs free energy change (ΔG) of the binding interaction, dissociation constant (Kd), number of interfacial contacts (ICs) between residues, and the non-interacting surface (NIS) percentage, which reflects the proportion of charged interface surface area not directly involved in binding.

1. Generate a **combined PDB of the entire binder-target complex**. Run script below on your local device, which takes in two inputs: (a) PDB of the target sequence and (b) PDB of the RFDiffusion-generated protein binder. The output of this script will be found in `$DIR_WORK/prodigy_input_pdb`.

- An example of the resulting PDB can be found in `example/prodigy_input_pdb/cycle1A_50diff_0.5temp.pdb`. You will find that the **peptide binder** (Chain B) has been merged to the target sequence on **ApoB-100** (Chain A).

2. On the web page, select the **PRODIGY (protein-protein)** setting and upload the PDB file. Set **"Interactor 1"** as A, and **"Interactor 2"** as B. 

3. Finally, click **Submit Prodigy**.

### (b) FlexPepDoc

You could also use FlexPepDoc to visualize the 3D structures in `$DIR_WORK/prodigy_input_pdb`.

1. Log into [FlexPepDoc ROSIE web server](https://r2.graylab.jhu.edu/auth/login?next=%2Fapps%2Fsubmit%2Fflexpepdock) via Github.
2. Upload a combined pdb file.
3. Specify **Docking partner** as "A_B".

![flexpepdoc](docs/flexpepdoc1.png)
**Fig 1**. Example of FlexPepDoc visualization.

### (c) PyMol

Optionally, you could visualize the entire binding complex (multiple peptides binding to target protein) on PyMOL (download latest version [HERE](https://www.pymol.org/)).

1. Generate the entire target protein structure on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2) (or download from a PDB repository).

    - Note: if the protein is super big (i.e. ApoB-100 is > 4000aa), then you don't need to generate the entire structure. Just generate a portion of the protein that is **sufficiently large enough to cover all the binding sites** of the peptides that you plan to visualize. See `example/pep1.pdb`.

2. Use this script to generate a **multi-PDB file**, which takes in the large PDB from step (1) and PDB files of peptide binders. An example of the output can be found in `example/pymol_pdb/cycle1_50diff_0.5temp_complex.pdb`.

```bash
    cycle="1" 
    diffusion="50"
    temp="0.5"
    DIR_WORK="/Users/keonapang/Desktop/NVIDIA/Sept12" 
    Rscript "./conversion_all.R" $cycle $diffusion $temp $DIR_WORK
```

3. Open PyMoL and enter the following code:

```bash
    # Cycle 1 - COMPLETE PRODUCT (APR 30, 2025)
    load /Users/keonapang/Desktop/NVIDIA/1_AlphaFold/pep1.pdb, scaffold
    load /Users/keonapang/Desktop/NVIDIA/5_pdb_original/cycle1A_1seqs_50diff_05temp.pdb, peptide1 # 
    load /Users/keonapang/Desktop/NVIDIA/5_pdb_original/cycle1B_1seqs_50diff_05temp.pdb, peptide2 # 
    load /Users/keonapang/Desktop/NVIDIA/5_pdb_original/cycle1DD_1seqs_50diff_05temp.pdb, peptide3 # note: this is actually 1C--
    load /Users/keonapang/Desktop/NVIDIA/5_pdb_original/cycle1D_1seqs_50diff_05temp.pdb, peptide4 #
    
    # Step 2: Align Chain A of each peptide to the scaffold
    align peptide1 and chain A, scaffold
    align peptide2 and chain A, scaffold
    align peptide3 and chain A, scaffold
    align peptide4 and chain A, scaffold
    
    # Step 3: Extract only the peptide chains (Chain B) and rename
    extract peptide1_chainB, peptide1 and chain B
    alter peptide1_chainB, chain="C"
    extract peptide2_chainB, peptide2 and chain B
    alter peptide2_chainB, chain="D"
    extract peptide3_chainB, peptide3 and chain B
    alter peptide3_chainB, chain="E"
    extract peptide4_chainB, peptide4 and chain B
    alter peptide4_chainB, chain="F"
    
    # Step 4: Merge scaffold and peptide chains into a single object
    create combined, scaffold or peptide1_chainB or peptide2_chainB or peptide3_chainB or peptide4_chainB
    save /Users/keonapang/Desktop/NVIDIA/scaffold_peptides_1ABCD_new.pdb, combined
```

![Fig 2.](docs/pymol2.png) **Fig 2**. Example of the 3D peptide-target binding complex rendering on PyMOL.

## The end

In summary, we started off with manually selecting binding sites on a target region (**cycle 1**) of the protein via SWISS-MODEL. Finally, we output 4 non-overlapping, unique peptide binders for this region.

![Fig 6](docs/Fig1_visualization.png) **Fig 3**. First 4 peptides (1A, 1B, 1C, 1D) binding to ApoB, visualized on PyMOL (left) and SWISS-MODEL (right).


## Version updates

| Date         | Update          |
|---------------|---------------------------------------------|
| Sept 13, 2025      | New parameter 'i' for iterations per target seq  |
| Sept 14, 2025      | Added AlphaMissense-multimer code     |

## Jupyter notebook version

A similar example of the workflow is available in jupyter notebook format [/scripts/protein-binder-design_Sept2025.ipynb](/scripts/protein-binder-design_Sept2025.ipynbscript/protein-binder-design_Sept2025.ipynb)

## Useful resources

- [dl_binder_design](https://github.com/nrbennet/dl_binder_design?tab=readme-ov-file#inf3)
- [alphafold](https://github.com/google-deepmind/alphafold)
- HADDOCK [web server](https://rascar.science.uu.nl/haddock2.4/)
- PRODIGY [web server](https://rascar.science.uu.nl/prodigy/)
- RosettaCommons [RFDiffusion](https://github.com/RosettaCommons/RFdiffusion)
- [ProteinPMNN](https://github.com/dauparas/ProteinMPNN)