# NVIDIA BioNeMo: Protein Binder Design for Drug Discovery

This workflow is designed for _in silico_ protein binder design by generating binder sequences and predicted structures for the binder and target. Unlike the original NVIDIA pipeline, this approach **does not** require running Alphafold2 on cloud GPU. 

Instead, you will first pre-compute the target protein structure on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2) (step #1), and save it as a `.PDB` file. This becomes the input for **RFDiffusion** (step #2), which generates the protein backbones for binder design. Subsequently, **ProteinMPNN** (step #3) will back-generate the amino acid sequence of this protein binder. Finally, this generated peptide structure is validated (step #4) on [PRODIGY](https://rascar.science.uu.nl/prodigy/) (Gibbs Free Energy) and Rosetta [FlexPepDoc](https://r2.graylab.jhu.edu/auth/login?next=%2Fapps%2Fsubmit%2Fflexpepdock). 

The example workflow illustrated below generates 8 peptide binders for target protein ApoB-100. We will also generate their multimer structures (in PDB format).

### System Requirements

- at least 1300 GB (1.3 TB) of fast NVMe SSD space
- modern CPU with at least 24 CPU cores
- at least 64 GB RAM
- 2 or more NVIDIA L40s, A100, or H100 GPUs

### Software Pre-requisites

- Python 3.11+

### Hardware requirements

- Total: 2 x GPU, 47 GiB GPU memory, 1.3TB GB SSD drive space, 60GiB RAM,24 CPU

    - **RFdiffusion** runs on 1 x GPU, ≥12 GiB GPU memory, 15GB free SSD drive space
    - **ProteinMPNN** runs on 1 x GPU, ≥3 GiB GPU memory, 10GB free SSD drive space


## 1. Manually select binding sites on target protein

![Fig 1. 3D model of ApoB-100 on SWISS-MODEL](docs/ApoB_3D.png) **Fig 1**. 3D model of ApoB-100 (target) on SWISS-MODEL

1. Get on the [SWISS-MODEL](https://swissmodel.expasy.org/) repository and in the search bar, type in your target protein of interest to view the interactive 3D model. In this example, we used [ApoB-100](https://swissmodel.expasy.org/repository/uniprot/P04114?template=9eag.1.A&range=38-4563) (**Fig 1**).
2. Define your 'binding window' size, which affects how long your generated peptide binder will be. In this example, we designed 8 binding windows of 40 amino acids in length.
3. Manually select these 8 binding sites on the 3D interactive model. In this example, we chose 4 peptides within the first half of the ApoB protein sequence (residues A91–357) and 4 peptides within the second half (residues A390–642). The full ApoB-100 protein backbone consists of 4,563 amino acids (**Fig 2**).

![Fig 2. Manually identify the binding sites](docs/ApoB_3D_seq.png) **Fig 2**. Manually identify the 8 target binding sites

**Table 1**. 8 target sequences of length 40 amino acids.

| Cycle         | ApoB-100 target sequence          | Amino Acid Position |
|---------------|---------------------------------------------|---------------------|
| 1A      | LKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYEL    | A91-130            |
| 1B      | EEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLG    | A170-209           |
| 1C      | VAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKIN    | A255-294           |
| 1D      | PKQAEAVLKTLQELKKLTISEQNIQRANLFNKLVTELRGL    | A318-357           |
| 2A      | CSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREI    | A390-429           |
| 2B      | GTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQT    | A459-498           |
| 2C      | LRKMEPKDKDQEVLLQTFLDDASPGDKRLAAYLMLMRSPS    | A531-570           |
| 2D      | EQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPT    | A587-626           |

## 2. AlphaFold2 Colab (free)

For each of the 8 target sequences (**Table 1**), compute their 3D structure (.PDB) on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2). On this colab notebook, perform the following steps for each sequence:

1. Input sequence under `query_sequence`.
2. Hit `Runtime` > `Run all` and wait ~ 5mins.
3. Download .zip results, unpack it and identify the model with this suffix `...rank_001_alphafold2_ptm_model_1_seed_000.pdb` (This is because AlphaFold2 automatically generates 5 possible structures, with the first-ranked structure being the one with highest confidence, based on [pLDDT](https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/)).
4. Rename only the first-ranked model file to `pep${cycle}.pdb` (i.e. if `cycle="1A"`, then this would be pep1A.pdb).
5. Store this model in a new working directory and assign it's path to `$DIR_WORK` (which will also be used later on).

```bash
    DIR_WORK="/your/local/working/directory" # DIR_WORK="/Users/keonapang/Desktop/NVIDIA/Sept12"
    mkdir -p $DIR_WORK
    cd $DIR_WORK # contains RFDiffusion and ProteinPMNN results
```

## 3. Launch a NVIDIA cloud virtual machine (VM)

1. Go to [brev.nvidia](https://brev.nvidia.com/) and click on this [Launchable](https://brev.nvidia.com/launchable/deploy/now?launchableID=env-32aLABBLqme9fNaaSdVL94Bollg). If you would like to create your own, click on Launchables in the menu bar, then `Create Launchables` and follow these settings:
    - Select "I have codes in a git repository" and enter the URL of this repo
    - Select 'VM-mode'
    - Click "Next" until you finally reach **select compute**. We recommended **A100 (80GiB) 2 GPUs x 24 CPUs | 240GiB  (80GiB GPU memory) ($3.96/hr)**

![NVIDIA VM settings](docs/NVIDIA_VM.png)


2. Click **"Deploy Launchable"** and **"Go to Instance Page"**. Wait ~10 minutes for VM to start

3. Enter the VM, then drag and drop the AlphaFold2 .pdb files from **step (2)** into your workspace

4. Start a new terminal session by clicking the "+" button at the top right.

5. An **NGC Personal API Key** is required to download and run any NVIDIA NIMs. If this is your first time, start by creating an account on [NGC](https://catalog.ngc.nvidia.com/). Then [generate the key](https://org.ngc.nvidia.com/setup/api-key) and note it down somewhere secure for future use.

```bash
    export NGC_CLI_API_KEY=<enter-key> 
    # Example: 
    #   export NGC_CLI_API_KEY=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
    docker login nvcr.io --username='$oauthtoken' --password="${NGC_CLI_API_KEY}"
```

6. Build inference models via docker container. Wait 20-30mins for the docker to build.

```bash
    # Install Dependencies
    sudo apt-get update # updated nvidia toolkit
    sudo apt-get install -y docker-compose # docker compose version 2+
    sudo apt install python3.11

    # The NIM cache allows you to download models and store previously-downloaded models on your local/server disk,
    # so that you don’t need to download them again later when you run the NIM again. 
    mkdir -p ~/.cache/nim
    chmod -R 777 ~/.cache/nim    
    export HOST_NIM_CACHE=~/.cache/nim

    # From the root of the cloned protein-binder-design repository:
    cd deploy/
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
    # Example:
        # CONTAINER ID   IMAGE                             COMMAND                   CREATED         STATUS         PORTS                                                             NAMES
        # 8ac6ad7cbb27   nvcr.io/nim/ipd/rfdiffusion:2.0   "/bin/sh -c 'exec \"$…"   2 minutes ago   Up 2 minutes   0.0.0.0:8082->8000/tcp, [::]:8082->8000/tcp                       protein-binder-design-rfdiffusion-1
        # 87cf41c4e7c6   nvcr.io/nim/ipd/proteinmpnn:1.0   "/bin/sh -c 'exec \"$…"   2 minutes ago   Up 2 minutes   6006/tcp, 8888/tcp, 0.0.0.0:8083->8000/tcp, [::]:8083->8000/tcp   protein-binder-design-proteinmpnn-1

    # 3. Health check 
    curl localhost:8082/v1/health/ready # RFdiffusion
    curl localhost:8083/v1/health/ready # Protein MPNN
    curl localhost:8084/v1/health/ready # AlphaFold multimer
    # Example:
        # {"status":"ready"}
```

## 4. Run RFDiffusion and ProteinMPNN

```bash
    # Define the cycle-to-target_sequence mapping. For simplicity, we use the first two target sequences as example.
    declare -A cycle_to_sequence=(
        ["1A"]="LKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYEL"  
        ["1B"]="EEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLG" 
    )

    contigs="15-20" # RFDiffusion input; sets the expected length of peptide to be between 15-20aa
    diffusion=50 # recommended range: 20-50
    temp=0.5 # recommended range range: 0.2 to 0.8
    num_seq=1 # one peptide binder per target sequence (do not modify)

    # Export API key
    export NGC_CLI_API_KEY=<enter-key> # Example: export NGC_CLI_API_KEY=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

    cycles=("1A" "1B")    # cycles=("1A" "1B" "1C" "1D" "2A" "2B" "2C" "2D")
    for cycle in "${cycles[@]}"; do
        target_sequence=${cycle_to_sequence[$cycle]}

        echo "Running script for $cycle..."
        dir="/home/ubuntu/protein-binder-design/scripts"
        python3.11 ${dir}/4_protein_binder_design.py  --cycle ${cycle} --num_seq ${num_seq} --diffusion ${diffusion} --temp ${temp} --target_sequence ${target_sequence} --contigs ${contigs}
    done
```

While results are being generated on the cloud (**Fig 3**), ensure you download them into `$DIR_WORK` directory on your local computer.

![results](docs/results.png)

**Fig 3**. Example of the output directory on the cloud VM.


> [!NOTE]
> Now, you may exit and shutdown the NVIDIA VM so it doesn't keep charging money.


## 5. Visualization and validation of binder-target (locally)

### (a) PRODIGY Gibbs Free Energy

The [PRODIGY web server](https://rascar.science.uu.nl/prodigy/) tool can provide a prediction of binder-target binding affinity. The outputs include the Gibbs free energy change (ΔG) of the binding interaction, dissociation constant (Kd), number of interfacial contacts (ICs) between residues, and the non-interacting surface (NIS) percentage, which reflects the proportion of charged interface surface area not directly involved in binding.

1. Generate a **combined PDB of the entire binder-target complex**. Run script below on your local device, which takes in two inputs: (a) PDB of the target sequence and (b) PDB of the RFDiffusion-generated protein binder. The output of this script will be found in `$DIR_WORK/prodigy_input_pdb`.

```bash
    cycle="1A"
    diffusion="50"
    temp="0.5" 
    DIR_WORK="/Users/keonapang/Desktop/NVIDIA/Sept12" # Example (please modify)
    
    root="/location/of/this/script" # # root="/Users/keonapang/Desktop/NVIDIA/Sept12"
    for cycle in "1A" "1B"; do
        Rscript "${root}/conversion.R" $cycle $diffusion $temp $DIR_WORK
    done
```

An example of the resulting PDB of the entire complex can be found in `example/prodigy_input_pdb/cycle1A_50diff_0.5temp.pdb`. You will find that the **peptide binder** (Chain B) has been merged to the target sequence on **ApoB-100** (Chain A).

2. On the web page, select the **PRODIGY (protein-protein)** setting and upload the PDB file. Set **"Interactor 1"** as A, and **"Interactor 2"** as B. 

3. Finally, click **Submit Prodigy**.

### (b) FlexPepDoc

You could also use FlexPepDoc to visualize the 3D structures in `$DIR_WORK/prodigy_input_pdb`.

1. Log into [FlexPepDoc ROSIE web server](https://r2.graylab.jhu.edu/auth/login?next=%2Fapps%2Fsubmit%2Fflexpepdock) via Github.
2. Upload a combined pdb file.
3. Specify **Docking partner** as "A_B".

![flexpepdoc](docs/flexpepdoc1.png)
**Fig 4**. Example of FlexPepDoc visualization.

### (c) PyMol

Optionally, you could visualize the entire binding complex (multiple peptides binding to target protein) on PyMOL (download latest version [HERE](https://www.pymol.org/)).

1. Generate the entire target protein structure on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2) (or download from a PDB repository). 

- Note that if the protein is super big (i.e. ApoB-100 is > 4000aa in length), then you don't need to generate the entire structure. Just generate a portion of the protein that is **sufficiently large enough to cover all the binding sites** of the peptides that you plan to visualize. See `example/pep1.pdb`.

2. Use this script to generate a **multi-PDB file**, which takes in the large PDB from step (1) and PDB files of peptide binders. An example of the output can be found in `example/pymol_pdb/cycle1_50diff_0.5temp_complex.pdb`.

```bash
    # In this example, I wish to visualize peptides 1A and 1B binding to ApoB-100
    cycle="1" 
    diffusion="50"
    temp="0.5"
    DIR_WORK="/Users/keonapang/Desktop/NVIDIA/Sept12" 
    root="/location/of/this/script" # root="/Users/keonapang/Desktop/NVIDIA/Sept12"

    Rscript "${root}/conversion_all.R" $cycle $diffusion $temp $DIR_WORK
```


3. Open PyMoL and visualize `cycle1_50diff_0.5temp_complex.pdb`.

![Fig 5.](docs/Fig1_visualization.png) **Fig 5**. First 4 peptides (cycle 1A, 1B, 1C, 1D) binding to ApoB, visualized on PyMOL (left) and SWISS-MODEL (right).

![Fig 6.](docs/Fig2_visualization.png)  **Fig 6**. Last 4 peptides (cycle 2A, 2B, 2C, 2D) binding to ApoB, visualized on PyMOL (left) and SWISS-MODEL (right).

## Notebook

A similar example of the workflow is available in jupyter notebook format [/scripts/protein-binder-design_Sept2025.ipynb](/scripts/protein-binder-design_Sept2025.ipynbscript/protein-binder-design_Sept2025.ipynb)

## Full NVIDIA pipeline

Downloading the [AlphaFold2](https://docs.nvidia.com/nim/bionemo/alphafold2/latest/prerequisites.html) and [AlphaFold2-Multimer](https://docs.nvidia.com/nim/bionemo/alphafold2-multimer/latest/quickstart-guide.html) model requires an additional 1250GB and 512GB of free SSD drive space respectively. Their download time (running the `docker compose pull` step) is also very long, up to 4-10 hours on 100+ Mbps internet connection for both models.

These models can:

- AlphaFold: Predict protein structure given a protein sequence
- AlphaFold2-Multimer: Predict protein structure given multiple protein sequences

