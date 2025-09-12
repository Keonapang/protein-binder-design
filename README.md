# NVIDIA BioNeMo: Protein Binder Design

The NVIDIA BioNeMo blueprint for _in silico_ protein binder design by automatically generating binder sequences and predicted structures for the binder and target.

This Blueprint takes as input an amino acid protein sequence. It utilizes the following models:

- **AlphaFold2**: for predicting protein structure from amino acid sequence.
- **ProteinMPNN**: for predicting amino acid sequences for protein backbones.
- **RFDiffusion**: of protein backbones for protein binder design.
- **AlphaFold2-Multimer**: for predicting protein structure of multimers from a list of amino acid sequences

Once completed, this Blueprint outputs predicted multimer structures (in PDB format) for the target protein sequence and any generated peptide binders. These binder-target multimeric structures can then be assessed to find binders that effectively bind the target protein.

## System Requirements

- At least 1300 GB (1.3 TB) of fast NVMe SSD space
- A modern CPU with at least 24 CPU cores
- At least 64 GB of RAM
- Two or more NVIDIA L40s, A100, or H100 GPUs

## Software Pre-requisites

- Python 3.11+

## Hardware requirements

- Total: 2 x GPU, 47 GiB GPU memory, 1.3TB GB SSD drive space, 60GiB RAM,24 CPU
- **RFdiffusion** runs on 1 x GPU, ≥12 GiB GPU memory, 15GB free SSD drive space
- **ProteinMPNN** runs on 1 x GPU, ≥3 GiB GPU memory, 10GB free SSD drive space

## Ensure that you have these files

-`protein-binder-design_v3.ipynb` from [this repo](scripts/protein-binder-design_v3.ipynb)

-`docker-compose.yaml` (3MB) from [this repo](deploy/docker-compose.yaml)

-`cycle1_alphafold2_output.pdb` (~80KB) pre-computed on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2)

-`cycle2_alphafold2_output.pdb` (~80KB) pre-computed on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2)


## Notebook

A detailed example of the workflow is located in [/scripts/protein-binder-design_Sept2025.ipynb](/scripts/protein-binder-design_Sept2025.ipynbscript/protein-binder-design_Sept2025.ipynb)

![Fig 1. 3D rendering of ApoB protein, with two distinct regions for binding](docs/ApoB_3D.png)

![Fig 2. Finding sequences](docs/ApoB_3D_seq.png)

![Fig 3. Visualization](docs/Fig1_visualization.png)

![Fig 4. Visualization 2](docs/Fig2_visualization.png)
