# Adaption of NVIDIA BioNeMo: Protein Binder Design

This workflow is designed for _in silico_ protein binder design by generating binder sequences and predicted structures for the binder and target. Unlike the original NVIDIA pipeline, this approach **does not** require running Alphafold2 on cloud GPU. 

Instead, you will first pre-compute the target protein structure on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2) (model #1), and save it as a `.PDB` file. This becomes the input for **RFDiffusion** (model #2), which generates the protein backbones for binder design. Subsequently, **ProteinMPNN** (model #3) will back-generate the amino acid sequence of this protein binder. Finally, this generated peptide structure is validated (step #4) on [PRODIGY](https://rascar.science.uu.nl/prodigy/) (Gibbs Free Energy) and Rosetta [FlexPepDoc](https://r2.graylab.jhu.edu/auth/login?next=%2Fapps%2Fsubmit%2Fflexpepdock). 

The example below generates 8 peptide binders for target protein ApoB-100.

### System Requirements

- at least **240 GB** of fast NVMe SSD space *(assuming you're not running AlphaFold2. If so, you will need at least 2TB)*
- at least **24 CPU** cores
- at least **64 GB RAM**
- 2 or more NVIDIA L40s, A100, or H100 GPUs

### Hardware requirements

- Total: 2 x GPU, 47 GiB GPU memory, 1.3TB GB SSD drive space, 60GiB RAM,24 CPU

    - **RFdiffusion** runs on 1 x GPU, ≥12 GiB GPU memory, 15GB free SSD drive space
    - **ProteinMPNN** runs on 1 x GPU, ≥3 GiB GPU memory, 10GB free SSD drive space

### Software Pre-requisites

- Python 3.11+


## 1. Ensure you have prepared two input files

This script takes in two arguments:

**1. Raw `.PDB` file** of the target protein, most likely downloaded through a protein database. The script will ignore all irrelevant rows in the PDB file except the 'ATOM' rows. *Note: the file doesn't have to begin with amino acid residue 1*.

```bash
HEADER    CYTOKINE                                09-JAN-07   2E7A              
TITLE     TNF RECEPTOR SUBTYPE ONE-SELECTIVE TNF MUTANT WITH ANTAGONISTIC       
JRNL       DOI    10.1074/JBC.M707933200                                       
REMARK   2 RESOLUTION.    1.80 ANGSTROMS.                                       
REMARK   3     
ATOM      1  N   PRO A   8      18.727  24.301  31.792  1.00 56.82           N  
ATOM      2  CA  PRO A   8      17.276  24.324  31.476  1.00 57.03           C  
```

**2. Space-delimited `.txt` file (no header)** containing four columns that store: chain identifier (e.g. A), hotspot residue (e.g. 80), start and end residue position (e.g. 60 and 90)

```bash
    A 80 60 90
    A 90 80 110
    A 103 83 113
    A 113 103 133
```

## 2. Launch a NVIDIA cloud virtual machine (VM)

1. Go to [brev.nvidia](https://brev.nvidia.com/) and click on this [Launchable](https://brev.nvidia.com/launchable/deploy/now?launchableID=env-32aLABBLqme9fNaaSdVL94Bollg). If you would like to create your own, click on Launchables in the menu bar, then `Create Launchables` and follow these settings:
    - Select "I have codes in a git repository" and enter the URL of this repo
    - Select 'VM-mode'
    - Click "Next" until you finally reach **select compute**. We recommended **A100 (80GiB) 2 GPUs x 24 CPUs | 240GiB  (80GiB GPU memory) ($3.96/hr)**

![NVIDIA VM settings](docs/NVIDIA_VM.png)


2. Click **"Deploy Launchable"** and **"Go to Instance Page"**. Wait ~10 minutes for VM to start

3. Enter the VM, then drag and drop the AlphaFold2 .pdb files from **step (2)** into the VM workspace, under `protein-binder-design/input/`

4. Start a new terminal session by clicking the "+" button at the top right.

5. An **NGC Personal API Key** is required to download and run any NVIDIA NIMs. If this is your first time, start by creating an account on [NGC](https://catalog.ngc.nvidia.com/). Then [generate the key](https://org.ngc.nvidia.com/setup/api-key) and note it down somewhere secure for future use.

```bash
    export NGC_CLI_API_KEY=<enter-key> 
    docker login nvcr.io --username='$oauthtoken' --password="${NGC_CLI_API_KEY}"
```

6. Build inference models via docker container. Wait 20-30mins for the docker to build.

```bash
    # Install Dependencies
    sudo apt-get update # updated nvidia toolkit
    sudo apt-get install -y docker-compose # docker compose version 2+
    sudo apt install python3.11

    # The NIM cache allows you to download models and store previously-downloaded models on your local/server disk,
    mkdir -p ~/.cache/nim
    chmod -R 777 ~/.cache/nim    
    export HOST_NIM_CACHE=~/.cache/nim

    # From the root of the cloned protein-binder-design repository:
    cd protein-binder-design/deploy/
    docker compose up # runs /deploy/docker-compose.yaml
```

7. Open a new terminal, leaving the current terminal open with the launched service.

8. In the new terminal, view running containers with the following codes. Wait until the health check returns `{"status":"ready"}` before proceeding.

```bash
    # 1. Check storage space
    df -h 

    # 2. Check which dockers are currently active/on stand-by
    docker container ls
    docker ps 

    # 3. Health check  # shoud be {"status":"ready"}
    curl localhost:8081/v1/health/ready # AlphaFold2
    curl localhost:8082/v1/health/ready # RFdiffusion
    curl localhost:8083/v1/health/ready # Protein MPNN

    docker logs -f protein-binder-design-alphafold-1
```

## 3. Download command-line tools

```bash
    pip install prodigy-prot torch Bio
    python3.11 -m pip install biopython pdb-tools
```

## 4. Run clean up and installation

```bash
    # clean up scripts before running
    sed -i 's/\r$//' "${REPO_DIR}/scripts/get_target_pdb.sh" # optional (to remove any hidden spaces from Windows)
    sed -i 's/\r$//' "${REPO_DIR}/scripts/calc_prodigy.sh" 
    sed -i 's/\r$//' "${REPO_DIR}/input/target_hotspots.txt" 
    sed -i 's/^[ \t]*//;s/[ \t]*$//' "${REPO_DIR}/scripts/calc_prodigy.sh"
```

Export API key again (assume it has already been generated)

```bash
    # export NGC_CLI_API_KEY=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
    export NGC_CLI_API_KEY=<enter-key>
```

## 5. Run script 

The command `get_target_pdb` checks that `start_pos` and `end_pos` are valid given the input PDB file. Then the script runs 3 prediction models sequentially, followed by aligning the designed peptide PDB to the original target sequence PDB to create a combined PDB using **BioPython's** `Superimposer` module. Lastly, it calculates the **dissociation constant (Kd)** using [PRODIGY](https://github.com/haddocking/prodigy):

- **RFDiffusion**: takes in target protein PDB, outputs a peptide binder PDB. Note: every output is a glycine, see [Github](https://github.com/RosettaCommons/RFdiffusion?tab=readme-ov-file#understanding-the-output-files) for more info.
    - `contigs`: range of amino acid positions, and expected length of peptide (i.e."A1-30/0 15-25")
    - `diffusion`: number of diffusion_steps (default: 50)
    - `i`: number of RFDiffusion iterations for each target sequence

- **ProteinPMNN**: takes in the output from RFDiffusion, and outputs a .fasta file storing a list of their predicted amino acid sequences
    - `num_seq`: how many seqs to generate for a given structure from RFDiffusion
    - `hotspot_res`: A20 # array (i.e. "A20")
    - `temp`: sample temperature controls the diversity of designed peptides. Higher values will lead to more diversity (range:0-1)

- **AlphaFold2**: takes in a list of peptide amino acid sequence(s) from ProteinMPNN, and predicts each of their PDB structures

```bash
    # Define repo directory  
    REPO_DIR="/home/shadeform/protein-binder-design"

    # Two input files
    raw_pdb="${REPO_DIR}/input/pdb2e7a.pdb" # input <- modify!
    input_file="${REPO_DIR}/input/target_hotspots.txt"
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

```bash
# Loop through each line of $input_file
while IFS=$' ' read -r chain hotspot_res_prefix start_pos end_pos; do # space-delimited

    # No need to edit these variables
    hotspot_res="${chain}${hotspot_res_prefix}"
    contigs="A${start_pos}-${end_pos}/0 ${peptide_length}$" # e.g. "A60-90/0 15-25"
    echo "Processing chain=$chain, hotspot_res=$hotspot_res, start_pos=$start_pos, end_pos=$end_pos"

    # Step 1: Build target structure PDB and extract target seq amino acid
    source "${REPO_DIR}/scripts/get_target_pdb.sh" # Load the get_target_pdb function
    get_target_pdb "${raw_pdb}" "${target_pdb}" "${chain}" "${start_pos}" "${end_pos}"
    target_pdb="${REPO_DIR}/input/target_${chain}${start_pos}_${end_pos}.pdb" # output
    target_sequence=$(bash "${REPO_DIR}/scripts/get_target_seq.sh" "${target_pdb}")
    
    # Step 2: Run the protein binder design script
    python3.11 "${REPO_DIR}/scripts/4_protein_binder_design.py" --num_seq "${num_seq}" --diffusion ${diffusion} --temp ${temp} --target_sequence "${target_sequence}" --contigs "${contigs}" --i "${i}" --hotspot_res ${hotspot_res} --target_pdb "${target_pdb}"

    # Step 3: Calculate binding free energy - adding the $raw_pdb argument is optional
    bash "${REPO_DIR}/scripts/calc_prodigy.sh" "${chain}" "${start_pos}" "${end_pos}" "${diffusion}" "${temp}" "${num_seq}" "${i}" "${target_sequence}" "${REPO_DIR}" "${raw_pdb}"

done < "$input_file"

```

### Workflow outputs

Example of **ProteinPMNN** .fasta output if `num_seq=2`:

```bash
    # binder_target_pairs
    [
        ['RIAELLAQLLKELLE','SQVLFSGQGCPSTHVLLTHTISRISTTHNQP'], # binder design 1
        ['AIEEALARLLLEQLL', 'SQVLFSGQGCPSTHVLLTHTISRISTTHNQP'] # binder design 2
    ]
```

Example of **final aligned binder-target PDB**:

```bash
# Final PDB complex of target protein and designed peptide binder
# Date: 2025-09-18
# Target sequence: ISRISTTHNQPVNLLSAIRSPCQRETPEGAE
# Parameters: 
#     target_pdb=/home/shadeform/protein-binder-design/input/pdb2e7a.pdb
#     binder_pdb=/home/shadeform/protein-binder-design/target_A60_90/4_target_A60_90_50diff_0.3temp_binder_i5_2.pdb
#     diffusion=50
#     temp=0.3
#     iteration=5 RFDiffusion candidates
#     num=2 seqs from ProteinMPNN per candidate
# 
# Predicted peptide binder from ProteinPMNN:
# >T=0.3, sample=2, score=1.6251, global_score=2.0374, seq_recovery=0.0000
# LTAEELAKLLAAAAAALALLALL
# PRODIGY results: 
#     Predicted binding affinity (kcal.mol-1): -24.3
#     Predicted dissociation constant (M) at 25.0˚C: 1.6e-18
#     No. of intermolecular contacts: 283
#     No. of charged-charged contacts: 8.0
#     No. of charged-polar contacts: 10.0
#     No. of charged-apolar contacts: 40.0
#     No. of polar-polar contacts: 21.0
#     No. of apolar-polar contacts: 81.0
#     No. of apolar-apolar contacts: 123.0
#     Percentage of apolar NIS residues: residues:
#     Percentage of charged NIS residues: residues:
# ATOM      1  N   PRO A   8      18.727  24.301  31.792  1.00 56.82           N  
# ATOM      2  CA  PRO A   8      17.276  24.324  31.476  1.00 57.03           C  
# ATOM      3  C   PRO A   8      16.970  25.033  30.160  1.00 59.51           C  
# ATOM      4  O   PRO A   8      17.158  26.246  30.040  1.00 60.32           O  

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
    # In this example, I wish to visualize peptides 1A and 1B binding to ApoB-100
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