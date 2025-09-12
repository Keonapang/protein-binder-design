# NVIDIA BioNeMo: Protein Binder Design for Drug Discovery

The NVIDIA BioNeMo blueprint for _in silico_ protein binder design by generating binder sequences and predicted structures for the binder and target. Unlike the original NVIDIA pipeline, this approach **does not** require running Alphafold2 on cloud GPU. 

Instead, you will first pre-compute the target protein structure on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2) (step #1), and save it as a `.PDB` file.

This is the input for **RFDiffusion** (step #2), which then generates the protein backbones for binder design. Subsequently, **ProteinMPNN** (step #3) will back-generate the amino acid sequence. Finally, the generated peptide structure is validated (step #4) on [PRODIGY](https://rascar.science.uu.nl/prodigy/) (Gibbs Free Energy) and Rosetta [FlexPepDoc](https://r2.graylab.jhu.edu/auth/login?next=%2Fapps%2Fsubmit%2Fflexpepdock). 

This workflow outputs peptide binders for each of the target protein sequences of interest. Optionally, you could also generate their multimer structures (in PDB format).

## System Requirements

- at least 1300 GB (1.3 TB) of fast NVMe SSD space
- modern CPU with at least 24 CPU cores
- at least 64 GB RAM
- 2 or more NVIDIA L40s, A100, or H100 GPUs

## Software Pre-requisites

- Python 3.11+

## Hardware requirements

- Total: 2 x GPU, 47 GiB GPU memory, 1.3TB GB SSD drive space, 60GiB RAM,24 CPU

    - **RFdiffusion** runs on 1 x GPU, ≥12 GiB GPU memory, 15GB free SSD drive space
    - **ProteinMPNN** runs on 1 x GPU, ≥3 GiB GPU memory, 10GB free SSD drive space

## Ensure that you have these files

- `protein-binder-design_v3.ipynb` from [this repo](scripts/protein-binder-design_v3.ipynb)

- `docker-compose.yaml` (3MB) from [this repo](deploy/docker-compose.yaml)

## 1. AlphaFold2

library(knitr)
library(dplyr)

df <- tribble(
  ~Sequence, ~Amino_acid_sequence, ~Position,
  "1A", "LKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYEL", "A91-130",
  "1B", "EEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLG", "A170-209",
  "1C", "VAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKIN", "A255-294",
  "2A", "CSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREI", "A390-429",
  "2B", "GTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQT", "A459-498",
  "2D", "EQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPT", "A587-626"
)

kable(df, col.names = c("Sequence", "Amino acid sequence", "Position"), align = "l")
* Sequence "1A": "LKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYEL" # A91-130
* Sequence "1B": "EEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLG"  #A170-209
* Sequence "1C": "VAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKIN" # A255-294
* Sequence "2A": "CSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREI", # A390-429
* Sequence "2B": "GTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQT", # A459-498
* Sequence "2D": "EQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPT" # A587-626


- `cycle1_alphafold2_output.pdb` (~80KB) pre-computed on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2)

- `cycle2_alphafold2_output.pdb` (~80KB) pre-computed on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2)


## Notebook

A detailed example of the workflow is located in [/scripts/protein-binder-design_Sept2025.ipynb](/scripts/protein-binder-design_Sept2025.ipynbscript/protein-binder-design_Sept2025.ipynb)

![Fig 1. 3D rendering of ApoB protein, with two distinct regions for binding](docs/ApoB_3D.png)

![Fig 2. Finding sequences](docs/ApoB_3D_seq.png)

![Fig 3. Visualization](docs/Fig1_visualization.png)

![Fig 4. Visualization 2](docs/Fig2_visualization.png)
