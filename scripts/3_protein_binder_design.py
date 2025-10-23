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
import time
import requests
from enum import StrEnum, Enum  # must be Python 3.11+
from typing import Tuple, Dict, Any, List
from pathlib import Path
import gc

# Load arguments
parser = argparse.ArgumentParser(description="De Novo Protein Design Workflow")
parser.add_argument("--num_seq", type=int, default=1, required=True, help="Number of sequences to generate per target")
parser.add_argument("--diffusion", type=int, default=50,help="Number of diffusion steps (15-30 recommended)")
parser.add_argument("--temp", type=float, default=0.2,help="Sampling temperature (range: 0-1)")
parser.add_argument("--target_sequence", type=str, required=True, help="Input sequence")
parser.add_argument("--contigs",type=str,required=True,help="Contigs string, e.g., 'A60-90/0 15-25'")
parser.add_argument("--i", type=int, required=True, default=1,help="iterations")
parser.add_argument("--hotspot_res", nargs='+',help="Hotspot residues (e.g., A67 A80)")
parser.add_argument("--target_pdb", type=str, required=True, help="Path to precomputed PDB file of target protein")
parser.add_argument("--cycle", type=str,help="Cycle/peptide name (e.g., '1', '1A', '1B', or '2')")
parser.add_argument("--root", type=str)
parser.add_argument("--chain", type=str)
args = parser.parse_args()

# Assign input arguments to variables
cycle = args.cycle
num_seq = args.num_seq
diffusion = args.diffusion
temp = args.temp
target_sequence = args.target_sequence
contigs = args.contigs
i = args.i
hotspot_res = args.hotspot_res
target_pdb = args.target_pdb
root = args.root
chain = args.chain

# target_pdb = "/home/shadeform/protein-binder-design/input/tnf1.pdb" ==
# num_seq = 1
# diffusion = 50
# temp = 0.4
# contigs="C88-128/0 15-25"
# target_sequence="QTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEK"
# i=2
# hotspot_res=['C108']
# root = "/home/shadeform/protein-binder-design"
# chain = "C"

print(f"Number of seqs to generate per target: {num_seq}"
      f"\nDiffusion: {diffusion}"
      f"\nSampling temp: {temp}"
      f"\nContigs: {contigs}"
      f"\nIterations: {i}"
      f"\nHotspot residues: {hotspot_res}\n"
)
# Set input directory 
# root = "/home/ubuntu/protein-binder-design"
os.makedirs(root, exist_ok=True)

# Set input variables and path
print(f"\nGenerating {num_seq} peptide binder(s) per input sequence...")
# Extract the filename without the extension
cycle = os.path.splitext(os.path.basename(target_pdb))[0]
name = f"{cycle}_{diffusion}diff_{temp}temp"

# Set output directory
# outdir = f"{root}/{diffusion}diff_{temp}temp"
outdir = f"{root}/{cycle}"
os.makedirs(outdir, exist_ok=True)
print(f"Output dir: {outdir}")
    
pdb_path = target_pdb
if not os.path.exists(pdb_path):
    raise FileNotFoundError(f"PDB file not found.\nMissing: {target_pdb}\n")

start_time = time.time()

##############################################################
# SET UP 
##############################################################
# Authenticator for NGC API

NVIDIA_API_KEY = os.getenv("NGC_CLI_API_KEY") # NVIDIA_API_KEY = os.getenv("NVIDIA_API_KEY") or input("Paste Run Key: ")
if not NVIDIA_API_KEY:
    raise ValueError("NGC_CLI_API_KEY environment variable is not set. Please export it before running script.")

# Connect to NIMs
HEADERS = {
    "Content-Type": "application/json",
    "Authorization": f"Bearer {NVIDIA_API_KEY}",
    "poll-seconds": "900"
}
NIM_HOST_URL_BASE = "http://localhost"

# 3 different endpoints for the models
NIM_HOST_URL_BASE = "http://localhost"

class NIM_PORTS(Enum):
    # ALPHAFOLD2_PORT = 8081
    RFDIFFUSION_PORT = 8082
    PROTEINMPNN_PORT = 8083
    # AF2_MULTIMER_PORT = 8084

class NIM_ENDPOINTS(StrEnum):
    # ALPHAFOLD2 = "protein-structure/alphafold2/predict-structure-from-sequence"
    RFDIFFUSION =  "biology/ipd/rfdiffusion/generate"
    PROTEINMPNN =  "biology/ipd/proteinmpnn/predict"
    # AF2_MULTIMER = "protein-structure/alphafold2/multimer/predict-structure-from-sequences"

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
    response = requests.post(function_url,json=payload,headers=headers)
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

if hotspot_res:  # Add hotspot_res only if it's provided
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
else:
    class ExampleRequestParams:
        def __init__(self,
                    target_sequence: str,
                    contigs: str, 
                    input_pdb_chains: List[str],
                    ca_only: bool,
                    use_soluble_model: bool,
                    sampling_temp: List[float],
                    diffusion_steps: int = 15,
                    num_seq_per_target: int = 20):
            self.target_sequence = target_sequence
            self.contigs = contigs
            self.hotspot_res = hotspot_res
            self.input_pdb_chains = input_pdb_chains or []  # Default to empty list
            self.ca_only = ca_only
            self.use_soluble_model = use_soluble_model
            self.sampling_temp = sampling_temp
            self.diffusion_steps = diffusion_steps
            self.num_seq_per_target = num_seq_per_target

status = check_nim_readiness(NIM_PORTS.RFDIFFUSION_PORT.value)
print(f"RFDiffusion ready: {status}")
status = check_nim_readiness(NIM_PORTS.PROTEINMPNN_PORT.value)
print(f"ProteinMPNN ready: {status}")

##############################################################
# Query code 
##############################################################
print(f"\n ================== {cycle} ==================\n")
chain="A"
example = ExampleRequestParams(
    target_sequence= target_sequence,
    contigs=contigs, 
    hotspot_res=hotspot_res, # hotspot_res=["A67", "A80"], optional 
    input_pdb_chains=[chain], # [Optional] default is to design for all chains in the protein
    ca_only=False, # [Optional]  CA-only model helps to address specific needs in protein design where focusing on the alpha carbon (CA)
    use_soluble_model=True, 
    sampling_temp=[temp], # (range: 0 - 1) adjust the probability values for the 20 amino acids at each position, controls the diversity of the design outcomes
    diffusion_steps=diffusion, # 15-30 steps are recommended for protein design tasks
    num_seq_per_target=num_seq  # Generate 4 binders
)

##############################################################
# 2. RFdiffusion
##############################################################

precomputed_pdb = get_reduced_pdb(pdb_path, rcsb_path=None)
print("\nPreview of input PDB file:\n")
print(precomputed_pdb[0:130])
print("")

# iterate through i iterations
for iteration in range(i):
    print(f"\n-----------[Iteration {iteration + 1} of {i}]-----------")
    outfile=f"{outdir}/2_{name}_i{iteration + 1}.pdb"
    
    if hotspot_res:  # Add hotspot_res only if it's provided
        print(f"1. Running RFdiffusion...")
        rfdiffusion_query = {
            "input_pdb": precomputed_pdb,  # Now using the precomputed PDB structure
            "contigs": contigs,
            "hotspot_res": hotspot_res,
            "diffusion_steps": diffusion
        }
    else: 
        print(f"1. Running RFdiffusion without hotspot_res ...")
        rfdiffusion_query = {
        "input_pdb": precomputed_pdb,  # Now using the precomputed PDB structure
        "contigs": contigs,
        "diffusion_steps": diffusion
        }
    rc, rfdiffusion_response = query_nim(
        payload=rfdiffusion_query,
        nim_endpoint=NIM_ENDPOINTS.RFDIFFUSION.value,
        nim_port=NIM_PORTS.RFDIFFUSION_PORT.value
    )
    print(rfdiffusion_response["output_pdb"][0:130])
    with open(outfile, "w") as pdb_file:
        pdb_file.write(rfdiffusion_response["output_pdb"])
    
    # clear python memory 
    del rfdiffusion_query
    gc.collect()

##############################################################
# 3. ProteinMPNN
##############################################################    
    print(f"\n2.Running ProteinMPNN to generate {num_seq} seq per target....")
    proteinmpnn_query = {
        "input_pdb" : rfdiffusion_response["output_pdb"],
        "input_pdb_chains" : example.input_pdb_chains,
        "ca_only" : example.ca_only,
        "use_soluble_model" : example.use_soluble_model,
        "num_seq_per_target" : num_seq,
        "sampling_temp" : example.sampling_temp
    }
    rc, proteinmpnn_response = query_nim(
        payload=proteinmpnn_query,
        nim_endpoint=NIM_ENDPOINTS.PROTEINMPNN.value,
        nim_port=NIM_PORTS.PROTEINMPNN_PORT.value
    )
    # Binder sequences are stored in fasta_sequences
    fasta_sequences = [x.strip() for x in proteinmpnn_response["mfasta"].split("\n") if '>' not in x][2:]
    binder_target_pairs = [[binder, example.target_sequence] for binder in fasta_sequences]
    # print(proteinmpnn_response["mfasta"])
    # Save binder_target_pairs as .json file
    fasta_sequences = []
    lines = proteinmpnn_response["mfasta"].split("\n")
    for i in range(len(lines)):
        if lines[i].startswith(">T="):  # Identify lines with binder headers
            if i + 1 < len(lines):  # Ensure the next line exists
                fasta_sequences.append(lines[i + 1].strip())  # Collect the sequence
    # Save proteinmpnn_response["mfasta"] to a .fasta file
    with open(f"{outdir}/3_{name}_i{iteration + 1}.fasta", "w") as fasta_file:
        fasta_file.write(proteinmpnn_response["mfasta"])
    binder_target_pairs = [[binder, example.target_sequence] for binder in fasta_sequences]
    # with open(f"{outdir}/3_{name}_proteinmpnn_pairs_i{iteration + 1}.json", "w") as json_file:
    #     json.dump(binder_target_pairs, json_file, indent=4)
    print(binder_target_pairs)
    print()
    # [['RIAELLAQLLKELLE', 'SQVLFSGQGCPSTHVLLTHTISRISTTHNQP'], ['AIEEALARLLLEQLL', 'SQVLFSGQGCPSTHVLLTHTISRISTTHNQP']]

    # clear python memory 
    del fasta_sequences
    del proteinmpnn_response
    del proteinmpnn_query
    gc.collect()
    ##############################################################
    # 3. AlphaFold to predict the structure of the binder alone
    ##############################################################
#     print(f"3. Running AlphaFold2 for {num_seq} seqs...")

#     # Preview binder_target_pairs
#     counter = 0
#     response_codes = [0 for i in binder_target_pairs]
#     results = [None for i in binder_target_pairs]

#     for binder_target_pair in binder_target_pairs:
#         output_file = os.path.join(outdir, f"4_{name}_binder_i{iteration + 1}_{counter+1}.pdb")
        
#         # Check if the output file already exists
#         if os.path.exists(output_file):
#             print(f"Output file already exists: {output_file}. Skipping this binder.")
#             counter += 1
#             continue  # Skip this iteration if the file exists
        
#         current_time = time.time()
#         predicted_binder = binder_target_pair[0]
#         print(f"\nDesigned binder ({counter+1} of {len(binder_target_pairs)}): {predicted_binder}")
#         alphafold2_query = {
#             "sequence" : predicted_binder,
#             "algorithm" : "mmseqs2",
#         }
#         rc, alphafold2_response = query_nim(
#             payload=alphafold2_query,
#             nim_endpoint=NIM_ENDPOINTS.ALPHAFOLD2.value,
#             nim_port=NIM_PORTS.ALPHAFOLD2_PORT.value
#         )
#         alphafold2_response[0][0:160]
#         if alphafold2_response and len(alphafold2_response) > 0:
#             first_structure = alphafold2_response[0]
            
#             with open(output_file, "w") as f:
#                 f.write(first_structure)
#             # print(f"Saved best structure to {output_file}")
#         else:
#             print("WARNING: no structure predictions found!")
#         # response_codes[counter] = rc
#         # results[counter] = alphafold2_response
#         end_time = time.time()
#         elapsed_minutes = (end_time - current_time) / 60
#         print(f"Time: {elapsed_minutes:.2f} mins")
#         counter += 1
#         if counter >= num_seq:
#             break

#         # clear python memory 
#         del first_structure
#         del alphafold2_response
#         del output_file
#         del alphafold2_query
#         gc.collect()

# print(f"Results saved in : {outdir}")
# end_time = time.time()
# elapsed_minutes = (end_time - start_time) / 60
# print(f"Total time: {elapsed_minutes:.2f} mins")

#     # Print probabilities and sequence scores
    # probs = proteinmpnn_response["probs"]
    # with open(f"{outdir}/3_{name}_proteinmpnn_probs_i{iteration + 1}.txt", "w") as probs_file:
    #     for i, prob_matrix in enumerate(probs):
    #         probs_file.write(f"Sequence {i+1}:\n")
    #         for position_probs in prob_matrix:
    #             probs_file.write(",".join(map(str, position_probs)) + "\n")
    #         probs_file.write("\n")
    # Save scores and probs to files
    # scores = proteinmpnn_response["scores"]
    # with open(f"{outdir}/3_{name}_proteinmpnn_scores_i{iteration + 1}.txt", "w") as scores_file:
    #     for i, score in enumerate(scores):
    #         scores_file.write(f"Sequence {i+1}: Score = {score}\n")
    
    
    ##############################################################
    # 4. AlphaFold-Multimer
    # inputs a binder-target pair (peptide from ProteinMPNN plus the original protein sequence used as input to this workflow).
    ##############################################################
    # print(f"Loading AlphaFold-Multimer...")

    # # Preview binder_target_pairs
    # n_processed = 0
    # multimer_response_codes = [0 for i in binder_target_pairs]
    # multimer_results = [None for i in binder_target_pairs]
    # pairs_to_process = 1

    # for binder_target_pair in binder_target_pairs:
    #     multimer_query = {
    #         "sequences" : binder_target_pair,
    #         "selected_models" : [1]
    #     }
    #     print(f"Processing binder-target pair {n_processed+1} of {len(binder_target_pairs)}")
    #     rc, multimer_response = query_nim(
    #         payload=multimer_query,
    #         nim_endpoint=NIM_ENDPOINTS.AF2_MULTIMER.value,
    #         nim_port=NIM_PORTS.AF2_MULTIMER_PORT.value
    #     )
    #     multimer_response_codes[n_processed] = rc
    #     multimer_results[n_processed] = multimer_response
    #     n_processed += 1
    #     if n_processed >= pairs_to_process:
    #         break

    # ## Print just the first 160 characters of the first multimer response
    # result_idx = 0
    # prediction_idx = 0
    # print(multimer_results[result_idx][prediction_idx][0:160])
    # print()

    # print(type(multimer_results))  # <class 'list'>
    # print(len(multimer_results))  # 1; Number of binder-target pairs processed
    # print(type(multimer_results[0])) # <class 'list'> (5 predictions)
    # print(len(multimer_results[0]))  # 5 (5 PDB predictions per pair)
    # # Check the first prediction for the first pair
    # print(type(multimer_results[0][0]))  # <class 'str'> (first PDB structure prediction)
    # print(multimer_results[0][0][0:160])

    # Save all PDB structures for each binder-target pair
    # for idx, predictions in enumerate(multimer_results):
    #     for pred_idx, pdb_structure in enumerate(predictions):
    #         output_file = os.path.join(outdir, f"4_AF2_complex_{cycle}_i{iteration+1}_{pred_idx+1}.pdb")
    #         with open(output_file, "w") as f:
    #             f.write(pdb_structure)
    #         print(f"Saved PDB structure for pair {idx+1}, prediction {pred_idx+1} to {output_file}")

    ##############################################################
    # 5. Analyze AlphaFold-Multimer results
    ##############################################################

    # # Function to calculate average pLDDT over all residues 
    # def calculate_average_pLDDT(pdb_string):
    #     total_pLDDT = 0.0
    #     atom_count = 0
    #     pdb_lines = pdb_string.splitlines()
    #     for line in pdb_lines:
    #         # PDB atom records start with "ATOM"
    #         if line.startswith("ATOM"):
    #             atom_name = line[12:16].strip() # Extract atom name
    #             if atom_name == "CA":  # consider atoms with name "CA"
    #                 try:
    #                     # Extract the B-factor value from columns 61-66 (following PDB format specifications)
    #                     pLDDT = float(line[60:66].strip())
    #                     total_pLDDT += pLDDT
    #                     atom_count += 1
    #                 except ValueError:
    #                     pass 
    #     if atom_count == 0:
    #         return 0.0  # Return 0 if no N atoms were found
    #     average_pLDDT = total_pLDDT / atom_count
    #     return average_pLDDT
    # plddts = []
    # for idx in range(0, len(multimer_results)):
    #     if multimer_results[idx] is not None:
    #         plddts.append(calculate_average_pLDDT(multimer_results[idx][0]))
    ##############################################################
    ## Combine the results with their pLDDTs, and print top 5 results
    # binder_target_results is a list of tuples, where each tuple contains:
        # - binder-target pair (binder_target_pairs[idx]).
        # - list of 5 PDB predictions (multimer_results[idx]).
        # - average pLDDT score (plddts[idx]).
    ##############################################################
    # ## Sort the results by plddt
    # binder_target_results = list(zip(binder_target_pairs, multimer_results, plddts))
    # sorted_binder_target_results = sorted(binder_target_results, key=lambda x : x[2])
    # for i in range(0, len(sorted_binder_target_results)):
    #     print("-"*80)
    #     print(f"rank: {i}")
    #     print(f"binder: {sorted_binder_target_results[i][0][0]}")
    #     print(f"target: {sorted_binder_target_results[i][0][1]}")
    #     print(f"pLDDT: {sorted_binder_target_results[i][2]}")
    #     print("-"*80)

    # Save PDB structure prediction for this binder-target pair
    # for idx, predictions in enumerate(sorted_binder_target_results):
    #     if predictions is not None and len(predictions) > 0:  # Ensure predictions exist
    #         # Get the first PDB structure (Prediction 1)
    #         first_pdb = predictions[0]  # list
    #         first_pdb_str = "\n".join(first_pdb)  # join list elements into a single string
    #         output_file = os.path.join(outdir, f"4_AF2_complex_{cycle}_i{iteration + 1}.pdb")
    #         with open(output_file, "w") as f:
    #             f.write(first_pdb_str)
    #         print(f"Saved first PDB structure for pair {idx+1} to {output_file}")
    # # Save pLDDT values for all 5 structures to .txt
    # for idx, (pair, predictions) in enumerate(zip(binder_target_pairs, multimer_results)):
    #     if predictions is not None:  
    #         structure_pLDDTs = [calculate_average_pLDDT(pdb) for pdb in predictions]
    #         output_file = os.path.join(outdir, f"4_AF2_complex_pLDDT_{cycle}_i{idx + 1}.txt")
    #         with open(output_file, "w") as f:
    #             f.write(f"AlphaFold2 Multimer\n")
    #             f.write(f"Binder: {pair[0]}\n")
    #             f.write(f"Target: {pair[1]}\n\n")
    #             f.write("pLDDTs for 5 predictions:\n")
    #             for i, plddt in enumerate(structure_pLDDTs):
    #                 f.write(f"Prediction {i + 1}: {plddt:.2f}\n")
    #         print(f"Saved pLDDT values for pair {idx + 1} to {output_file}")

            