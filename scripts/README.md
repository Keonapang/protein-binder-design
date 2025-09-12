# NVIDIA BioNeMo Blueprint for Protein Binder Design in Drug Discovery

This directory contains a Jupyter notebook showing a complete example. Check out the
[notebook](./protein-binder-design.ipynb) itself for more explanation of
how the protein binder design pipeline works.

Deploy the NIMs using the provided docker compose yaml found in [deploy](../deploy), you should be able to start up a Jupyter notebook instance and get going immediately!

```bash
jupyter notebook
```


## Software Pre-requisites

  -Python 3.11+

## Hardware requirements

Total: 2 x GPU, 47 GiB GPU memory, 1.3TB GB SSD drive space, 60GiB RAM,24 CPU
  -**RFdiffusion** runs on 1 x GPU, ≥12 GiB GPU memory, 15GB free SSD drive space
  -**ProteinMPNN** runs on 1 x GPU, ≥3 GiB GPU memory, 10GB free SSD drive space

## Ensure that you have these files

  -`protein-binder-design_v3.ipynb` uploaded to a [public Github repo](https://github.com/Keonapang/generative-protein-binder-design/blob/main/src/protein-binder-design.ipynb)
  -`docker-compose.yaml` (3MB) from the original [BioNeMo repo](https://github.com/NVIDIA-BioNeMo-blueprints/generative-protein-binder-design/blob/main/deploy/docker-compose.yaml)
  -`cycle1_alphafold2_output.pdb` (~80KB) pre-computed on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2)
  -`cycle2_alphafold2_output.pdb` (~80KB) pre-computed on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2)