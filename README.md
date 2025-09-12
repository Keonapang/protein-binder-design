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

For each of the 8 target sequences (**Table 1**), compute their 3D structure (.PDB) on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2). Perform these steps for each sequence:

1. On the AlphaFold2 colab notebook, input sequence under `query_sequence`
2. Hit `Runtime` > `Run all` and wait ~5-10mins
3. Download .zip results, decompress it and save the model with this suffix `...rank_001_alphafold2_ptm_model_1_seed_000.pdb` (this is because AlphaFold2 automatically generates 5 possible structures, with the first model being the one with the highest confidence, based on [pLDDT](https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/)).


## 3. RFDiffusion

- `protein-binder-design_v3.ipynb` from [this repo](scripts/protein-binder-design_v3.ipynb)

- `docker-compose.yaml` (3MB) from [this repo](deploy/docker-compose.yaml)




## 4.



## 5. Visualization and validation


## Notebook

A detailed example of the workflow is located in [/scripts/protein-binder-design_Sept2025.ipynb](/scripts/protein-binder-design_Sept2025.ipynbscript/protein-binder-design_Sept2025.ipynb)

![Fig 3. Visualization](docs/Fig1_visualization.png)

![Fig 4. Visualization 2](docs/Fig2_visualization.png)
