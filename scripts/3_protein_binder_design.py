# # De Novo Protein Design Workflow using NIMs (skipped AlphaFold2 steps)
#   modified to accept a precomputed AlphaFold2 PDB structure
# Keona Pang
# Apr 29, 2025

# Hardware requirements: 
#   Run on 2 x L40S GPUs x 16 CPUs | 294GiB; 128 GiB storage CRUSOE
#   doesn't require updating docker, just reuqires python3.11 
##########################################################################

import argparse
import json
import os
import requests
from enum import StrEnum, Enum  # must be Python 3.11+
from typing import Tuple, Dict, Any, List
from pathlib import Path

# Load arguments
parser = argparse.ArgumentParser(description="De Novo Protein Design Workflow")
parser.add_argument("--cycle", type=str, required=True, help="Cycle number (e.g., '1', '1A', '1B', or '2')")
parser.add_argument("--num_seq", type=int, default=1, help="Number of sequences to generate per target")
parser.add_argument("--diffusion", type=int, default=20, help="Number of diffusion steps (15-30 recommended)")
parser.add_argument("--temp", type=float, default=0.2, help="Sampling temperature (range: 0-1)")
args = parser.parse_args()

# Assign input arguments to variables
cycle = args.cycle
num_seq = args.num_seq
diffusion = args.diffusion
temp = args.temp

# cycle = "1A"
# num_seq = 1
# diffusion = 40
# temp = 0.4

# Set up all paths and variables
root = "/home/ubuntu/nvidia-workbench"
os.makedirs(root, exist_ok=True)
print(f"Generating {num_seq} sequences per target for cycle {cycle}...")

if num_seq > 0: # if num_seq > 1, then run code below
    contigs="15-25"
    precomputed_pdb_path = f"/home/ubuntu/pep{cycle}.pdb" 
    if "1" in cycle:
        if cycle == "1A": 
            target_sequence="LKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYEL" # A91-130
        if cycle == "1B":
            target_sequence="EEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLG"  #A170-209
        if cycle == "1C":
            target_sequence="VAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKIN" # A255-294
        if cycle == "1D":
            target_sequence="PKQAEAVLKTLQELKKLTISEQNIQRANLFNKLVTELRGL" # A318-357
    elif "2" in cycle:
        if cycle == "2A": 
            target_sequence="CSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREI", # A390-429
        if cycle == "2B":
            target_sequence="GTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQT", # A459-498
        if cycle == "2C":
            target_sequence="LRKMEPKDKDQEVLLQTFLDDASPGDKRLAAYLMLMRSPS", # A531-570
        if cycle == "2D":
            target_sequence = "EQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPT" # A587-626
    else:
        raise ValueError("Invalid cycle number.")

# Set up variables part 2
name = f"cycle{cycle}_{num_seq}seqs_{diffusion}diff_{temp}temp"
outdir = f"{root}/{diffusion}diff_{temp}temp_{num_seq}seq"
os.makedirs(outdir, exist_ok=True)
print(f"Output dir : {outdir}")

# check if this outdir exists
if os.path.exists(outdir):
    print(f"Output dir {outdir} already exists. Please check the path.")

# list all files in outdir
print("Current directory contents:")
for file in os.listdir(outdir):
    print(file)
    
# Check if AlphaFold2 PDB exists
if not os.path.exists(precomputed_pdb_path):
    raise FileNotFoundError(f"Precomputed PDB file {precomputed_pdb_path} does not exist.")

##############################################################
# SET UP 
##############################################################
# Authenticator for NGC API

NVIDIA_API_KEY = os.getenv("NGC_CLI_API_KEY") # NVIDIA_API_KEY = os.getenv("NVIDIA_API_KEY") or input("Paste Run Key: ")
if not NVIDIA_API_KEY:
    raise ValueError("NGC_CLI_API_KEY environment variable is not set. Please export it before running the script.")

# Connect to NIMs
HEADERS = {
    "Content-Type": "application/json",
    "Authorization": f"Bearer {NVIDIA_API_KEY}",
    "poll-seconds": "900"
}
NIM_HOST_URL_BASE = "http://localhost"

# 3 different endpoints for the models
class NIM_PORTS(Enum):
    RFDIFFUSION_PORT = 8082
    PROTEINMPNN_PORT = 8083
class NIM_ENDPOINTS(StrEnum):
    RFDIFFUSION =  "biology/ipd/rfdiffusion/generate"
    PROTEINMPNN =  "biology/ipd/proteinmpnn/predict"

def query_nim(
            payload: Dict[str, Any],
            nim_endpoint: str,
            headers: Dict[str, str] = HEADERS,
            base_url: str = "http://localhost",
            nim_port: int = 8080,
            echo: bool = False) -> Tuple[int, Dict]:
    function_url = f"{base_url}:{nim_port}/{nim_endpoint}"
    if echo:
        print("*"*80)
        print(f"\tURL: {function_url}")
        print(f"\tPayload: {payload}")
        print("*"*80)
    response = requests.post(function_url,
                            json=payload,
                            headers=headers)
    if response.status_code == 200:
        return response.status_code, response.json()
    else:
        raise Exception(f"Error: response returned code [{response.status_code}], with text: {response.text}")

def check_nim_readiness(nim_port: NIM_PORTS,
                        base_url: str = NIM_HOST_URL_BASE,
                        endpoint: str = "v1/health/ready") -> bool:
    """
    Return true if a NIM is ready.
    """
    try:
        response = requests.get(f"{base_url}:{nim_port}/{endpoint}")
        d = response.json()
        if "status" in d:
            if d["status"] == "ready":
                return True
        return False
    except Exception as e:
        print(e)
        return False
    
def get_reduced_pdb(pdb_id: str, rcsb_path: str = None) -> str:
    pdb = Path(pdb_id)
    if not pdb.exists() and rcsb_path is not None:
        pdb.write_text(requests.get(rcsb_path).text)
    lines = filter(lambda line: line.startswith("ATOM"), pdb.read_text().split("\n"))
    return "\n".join(list(lines))

class ExampleRequestParams:
    def __init__(self,
                target_sequence: str,
                contigs: str, 
                hotspot_res: List[str],
                input_pdb_chains: List[str],
                ca_only: bool,
                use_soluble_model: bool,
                sampling_temp: List[float],
                diffusion_steps: int = 15,
                num_seq_per_target: int = 20):
        self.target_sequence = target_sequence
        self.contigs = contigs
        self.hotspot_res = hotspot_res
        self.input_pdb_chains = input_pdb_chains
        self.ca_only = ca_only
        self.use_soluble_model = use_soluble_model
        self.sampling_temp = sampling_temp
        self.diffusion_steps = diffusion_steps
        self.num_seq_per_target = num_seq_per_target
status = check_nim_readiness(NIM_PORTS.RFDIFFUSION_PORT.value)
print(f"RFDiffusion ready: {status}")
status = check_nim_readiness(NIM_PORTS.PROTEINMPNN_PORT.value)
print(f"ProteinMPNN ready: {status}")
print()
print(f"------------- Cycle {cycle} ------------------")

##############################################################
# Query code 
##############################################################
cycle = ExampleRequestParams(
    target_sequence= target_sequence,
    contigs=contigs, 
    hotspot_res=[], 
    input_pdb_chains=["A"], # [Optional] default is to design for all chains in the protein
    ca_only=False, # [Optional]  CA-only model helps to address specific needs in protein design where focusing on the alpha carbon (CA)
    use_soluble_model=True, 
    sampling_temp=[temp], # (range: 0 - 1) adjust the probability values for the 20 amino acids at each position, controls the diversity of the design outcomes
    diffusion_steps=diffusion, # 15-30 steps are recommended for protein design tasks
    num_seq_per_target=num_seq  # Generate 4 binders
)
example=cycle

##############################################################
# 2. RFdiffusion
##############################################################
precomputed_pdb=get_reduced_pdb(precomputed_pdb_path, rcsb_path=None)

print(f"Running RFdiffusion....")
rfdiffusion_query = {
    "input_pdb": precomputed_pdb,  # Now using the precomputed PDB structure
    "contigs": example.contigs,
    "diffusion_steps": example.diffusion_steps
}
rc, rfdiffusion_response = query_nim(
    payload=rfdiffusion_query,
    nim_endpoint=NIM_ENDPOINTS.RFDIFFUSION.value,
    nim_port=NIM_PORTS.RFDIFFUSION_PORT.value
)

print(rfdiffusion_response["output_pdb"][0:160])
with open(f"{outdir}/2_{name}_rfdiffusion.pdb", "w") as pdb_file:
    pdb_file.write(rfdiffusion_response["output_pdb"])

##############################################################
# 3. ProteinMPNN
##############################################################
print()
print(f"Running ProteinMPNN....")
proteinmpnn_query = {
    "input_pdb" : rfdiffusion_response["output_pdb"],
    "input_pdb_chains" : example.input_pdb_chains,
    "ca_only" : example.ca_only,
    "use_soluble_model" : example.use_soluble_model,
    "num_seq_per_target" : example.num_seq_per_target,
    "sampling_temp" : example.sampling_temp
}
rc, proteinmpnn_response = query_nim(
    payload=proteinmpnn_query,
    nim_endpoint=NIM_ENDPOINTS.PROTEINMPNN.value,
    nim_port=NIM_PORTS.PROTEINMPNN_PORT.value
)

# binder sequences are stored in fasta_sequences
fasta_sequences = [x.strip() for x in proteinmpnn_response["mfasta"].split("\n") if '>' not in x][2:]
binder_target_pairs = [[binder, example.target_sequence] for binder in fasta_sequences]
print()
print(fasta_sequences)
print()
print(proteinmpnn_response["mfasta"])
print()
print(proteinmpnn_response)
print()

# Save binder_target_pairs as .json file
fasta_sequences = []
lines = proteinmpnn_response["mfasta"].split("\n")
for i in range(len(lines)):
    if lines[i].startswith(">T="):  # Identify lines with binder headers
        if i + 1 < len(lines):  # Ensure the next line exists
            fasta_sequences.append(lines[i + 1].strip())  # Collect the sequence

binder_target_pairs = [[binder, example.target_sequence] for binder in fasta_sequences]
with open(f"{outdir}/3_{name}_proteinmpnn_pairs.json", "w") as json_file:
    json.dump(binder_target_pairs, json_file, indent=4)
print(binder_target_pairs)
print()

# Save proteinmpnn_response["mfasta"] to a .fasta file
with open(f"{outdir}/3_{name}_proteinmpnn.fasta", "w") as fasta_file:
    fasta_file.write(proteinmpnn_response["mfasta"])

# probs = proteinmpnn_response["probs"]
# with open(f"{outdir}/3_{name}_proteinmpnn_probs.txt", "w") as probs_file:
#     for i, prob_matrix in enumerate(probs):
#         probs_file.write(f"Sequence {i+1}:\n")
#         for position_probs in prob_matrix:
#             probs_file.write(",".join(map(str, position_probs)) + "\n")
#         probs_file.write("\n")

# Save scores and probs to files
# scores = proteinmpnn_response["scores"]
# with open(f"{outdir}/3_{name}_proteinmpnn_scores.txt", "w") as scores_file:
#     for i, score in enumerate(scores):
#         scores_file.write(f"Sequence {i+1}: Score = {score}\n")
