# # De Novo Protein Design Workflow using NIMs (skipped AlphaFold2 prediction step)
#   modified to accept a precomputed AlphaFold2 PDB structure
# Keona Pang
# Apr 18, 2025

# Hardware requirements: 4 x A100 CRUSOE instance
#   4 x GPU, 47 GiB GPU memory
#   128 GB SSD drive space, 60GiB RAM,24 CPU

# INPUT FILES TO UPLOAD TO JUPYTER NOTEBOOK
#    - `protein-binder-design_v3.ipynb` uploaded to a [public Github repo](https://github.com/Keonapang/generative-protein-binder-design/blob/main/src/protein-binder-design.ipynb)
#    - `docker-compose.yaml` (3MB) from the original [BioNeMo repo](https://github.com/NVIDIA-BioNeMo-blueprints/generative-protein-binder-design/blob/main/deploy/docker-compose.yaml)
#    - `cycle1_alphafold2_output.pdb` (80KB) pre-computed on [AlphaFold2 colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=R_AH6JSXaeb2)

# EXECUTION
# cycle = 1 # cycle 1 or cycle 2
# num_seq = 4 # number of sequences to generate per target
# diffusion = 30 # number of diffusion steps (15-30 recommended)
# temp = 0.2 # sampling temperature (range: 0-1) to adjust the probability values for the 20 amino acids at each position, controls the diversity of the design outcomes

# brew install brevdev/homebrew-/brev && brev login --token <****> # Install the CLI
# brev shell <instance-name< # find instance-name on brev.nvidias # Open a terminal locally

# update
# sudo apt-get install -y docker-compose
# sudo apt install python3.11
# export NGC_CLI_API_KEY=<enter-key> # enter personal API key
# docker login nvcr.io --username='$oauthtoken' --password="${NGC_CLI_API_KEY}"

# ## Create the nim cache directory to download any model data to your local/server disk 
# mkdir -p ~/.cache/nim
# chmod -R 777 ~/.cache/nim    ## Make it writable by the NIM
# export HOST_NIM_CACHE=~/.cache/nim

# # Run Docker compose
# docker compose 

# pip install requests
# python SCRIPT_protein_binder_design.py --cycle 1 --num_seq 4 --diffusion 30 --temp 0.2

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

# Set up all paths and variables
root = "/home/ubuntu/nvidia-workbench"
os.makedirs(root, exist_ok=True)
print(f"Generating {num_seq} sequences per target for cycle {cycle}...")

if num_seq > 1: # if num_seq > 1, then run code below
    if cycle == "1":
        contigs="A400-600/0 15-25"
        precomputed_pdb_path = "cycle1_alphafold2_output.pdb" # converts pdb_id to a Path Object;
        target_sequence="MDPPRPALLALLALPALLLLLLAGARAEEEMLENVSLVCPKDATRFKHLRKYTYNYEAESSSGVPGTADSRSATRINCKVELEVPQLCSFILKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYELKLAIPEGKQVFLYPEKDEPTYILNIKRGIISALLVPPETEEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLGQCDRFKPIRTGISPLALIKGMTRPLSTLISSSQSCQYTLDAKRKHVAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKINSRFFGEGTKKMGLAFESTKSTSPPKQAEAVLKTLQELKKLTISEQNIQRANLFNKLVTELRGLSDEAVTSLLPQLIEVSSPITLQALVQCGQPQCSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREIFNMARDQRSRATLYALSHAVNNYHKTNPTGTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQTMEQLTPELKSSILKCVQSTKPSLMIQKAAIQALRKMEPKDKDQEVLLQTFLDDASPGDKRLAAYLMLMRSPSQADINKIVQILPWEQNEQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPTVMDFRKFSRNYQLYKSVSLPSLDPASAKIEGNLIFDPNNYLPKESMLKTTLTAFGFASADLIEIGLEGKGFEPTLEALFGKQGFFPDSVNKALYWVNGQVPDGVSKVLVDHFGYTKDDKHEQDMVNGIMLSVEKLIKDLKSKEVPEARAYLRILGEELGFASLHDLQLLGKLLLMGARTLQGIPQMIGEVIRKGSKNDFFLHYIFMENAFELPTGAGLQLQISSSGVIAPGAKAGVKLEVANMQAELVAKPSVSVEFVTNMGIIIPDFARSGVQMNTNFFHESGLEAHVALKAGKLKFIIPSPKRPVKLLSGGNTLHLVSTTKTEVIPPLIENRQSWSVCKQVFPGLNYCTSGAYSNASSTDSASYYPLTGDTRLELELRPTGEIEQYSVSATYELQREDRALVDTLKFVTQAEGAKQTEATMTFKYNRQSMTLSSEVQIPDFDVDLGTILRVNDESTEGKTSYRLTLDIQNKKITEVALMGHLSCDTKEERKIKGVISIPRLQAEARSEILAHWSPAKLLLQMDSSATAYGSTVSKRVAWHYDEEKIEFEWNTGTNVDTKKMTSNFPVDLSDYPKSLHMYANRLLDHRVPQTDMTFRHVGSKLIVAMSSWLQKASGSLPYTQTLQDHLNSLKEFNLQNMGLPDFHIPENLFLKSDGRVKYTLNKN"
    elif cycle == "2":
        contigs="A100-300/0 15-25"
        target_sequence="GKIDFLNNYALFLSPSAQQASWQVSARFNQYKYNQNFSAGNNENIMEAHVGINGEANLDFLNIPLTIPEMRLPYTIITTPPLKDFSLWEKTGLKEFLKTTKQSFDLSVKAQYKKNKHRHSITNPLAVLCEFISQSIKSFDRHFEKNRNNALDFVTKSYNETKIKFDKYKAEKSHDELPRTFQIPGYTVPVVNVEVSPFTIEMSAFGYVFPKAVSMPSFSILGSDVRVPSYTLILPSLELPVLHVPRNLKLSLPDFKELCTISHIFIPAMGNITYDFSFKSSVITLNTNAELFNQSDIVAHLLSSSSSVIDALQYKLEGTTRLTRKRGLKLATALSLSNKFVEGSHNSTVSLTTKNMEVSVATTTKAQIPILRMNFKQELNGNTKSKPTVSSSMEFKYDFNSSMLYSTAKGAVDHKLSLESLTSYFSIESSTKGDVKGSVLSREYSGTIASEANTYLNSKSTRSSVKLQGTSKIDDIWNLEVKENFAGEATLQRIYSLWEHSTKNHLQLEGLFFTNGEHTSKATLELSPWQMSALVQVHASQPSSFHDFPDLGQEVALNANTKNQKIRWKNEVRIHSGSFQSQVELSNDQEKAHLDIAGSLEGHLRFLKNIILPVYDKSLWDFLKLDVTTSIGRRQHLRVSTAFVYTKNPNGYSFSIPVKVLADKFIIPGLKLNDLNSVLVMPTFHVPFTDLQVPSCKLDFREIQIYKKLRTSSFALNLPTLPEVKFPEVDVLTKYSQPEDSLIPFFEITVPESQLTVSQFTLPKSVSDGIAALDLNAVANKIADFELPTIIVPEQTIEIPSIKFSVPAGIVIPSFQALTARFEVDSPVYNATWSASLKNKADYVETVLDSTCSSTVQFLEYELNVLGTHKIEDGTLASKTKGTFAHRDFSAEYEEDGKYEGLQEWEGKAHLNIKSPAFTDLHLRYQKDKKGISTSAASPAVGTVGMDMDEDDDFSKWNFYYSPQSSPDKKLTIFKTELRVRESDEETQIKVNWEEEAASGLLTSLKDNVPKATGVLYDYVNKYHWEHTGLTLREVSSKLRRNLQNNAEWVYQGAIRQIDDIDVRFQKAASGTTGTYQEWKDKAQNLYQELLTQEGQASFQGLKDNVFDGLVRVTQEFHMKVKHLIDSLIDFLNFPRFQFPGKPGIYTREELCTMFIREVGTVLSQVYSKVHNGSEILFSYFQDLVITLPFELRKHKLIDV"
        precomputed_pdb_path = "cycle2_alphafold2_output.pdb" # 
    else:
        raise ValueError("Invalid cycle number. Please set cycle to 1 or 2.")
    
if num_seq == 1: # if num_seq > 1, then run code below
    
    if "1" in cycle:
        target_sequence="MDPPRPALLALLALPALLLLLLAGARAEEEMLENVSLVCPKDATRFKHLRKYTYNYEAESSSGVPGTADSRSATRINCKVELEVPQLCSFILKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYELKLAIPEGKQVFLYPEKDEPTYILNIKRGIISALLVPPETEEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLGQCDRFKPIRTGISPLALIKGMTRPLSTLISSSQSCQYTLDAKRKHVAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKINSRFFGEGTKKMGLAFESTKSTSPPKQAEAVLKTLQELKKLTISEQNIQRANLFNKLVTELRGLSDEAVTSLLPQLIEVSSPITLQALVQCGQPQCSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREIFNMARDQRSRATLYALSHAVNNYHKTNPTGTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQTMEQLTPELKSSILKCVQSTKPSLMIQKAAIQALRKMEPKDKDQEVLLQTFLDDASPGDKRLAAYLMLMRSPSQADINKIVQILPWEQNEQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPTVMDFRKFSRNYQLYKSVSLPSLDPASAKIEGNLIFDPNNYLPKESMLKTTLTAFGFASADLIEIGLEGKGFEPTLEALFGKQGFFPDSVNKALYWVNGQVPDGVSKVLVDHFGYTKDDKHEQDMVNGIMLSVEKLIKDLKSKEVPEARAYLRILGEELGFASLHDLQLLGKLLLMGARTLQGIPQMIGEVIRKGSKNDFFLHYIFMENAFELPTGAGLQLQISSSGVIAPGAKAGVKLEVANMQAELVAKPSVSVEFVTNMGIIIPDFARSGVQMNTNFFHESGLEAHVALKAGKLKFIIPSPKRPVKLLSGGNTLHLVSTTKTEVIPPLIENRQSWSVCKQVFPGLNYCTSGAYSNASSTDSASYYPLTGDTRLELELRPTGEIEQYSVSATYELQREDRALVDTLKFVTQAEGAKQTEATMTFKYNRQSMTLSSEVQIPDFDVDLGTILRVNDESTEGKTSYRLTLDIQNKKITEVALMGHLSCDTKEERKIKGVISIPRLQAEARSEILAHWSPAKLLLQMDSSATAYGSTVSKRVAWHYDEEKIEFEWNTGTNVDTKKMTSNFPVDLSDYPKSLHMYANRLLDHRVPQTDMTFRHVGSKLIVAMSSWLQKASGSLPYTQTLQDHLNSLKEFNLQNMGLPDFHIPENLFLKSDGRVKYTLNKN",
        precomputed_pdb_path = "cycle1_alphafold2_output.pdb" # 

        if cycle == "1A": 
            contigs="A400-440/0 15-25"
        if cycle == "1B":
            contigs="A450-490/0 15-25"
        if cycle == "1C":
            contigs="A500-540/0 15-25"
        if cycle == "1D":
            contigs="A550-590/0 15-25"
    
    elif "2" in cycle:
        target_sequence="GKIDFLNNYALFLSPSAQQASWQVSARFNQYKYNQNFSAGNNENIMEAHVGINGEANLDFLNIPLTIPEMRLPYTIITTPPLKDFSLWEKTGLKEFLKTTKQSFDLSVKAQYKKNKHRHSITNPLAVLCEFISQSIKSFDRHFEKNRNNALDFVTKSYNETKIKFDKYKAEKSHDELPRTFQIPGYTVPVVNVEVSPFTIEMSAFGYVFPKAVSMPSFSILGSDVRVPSYTLILPSLELPVLHVPRNLKLSLPDFKELCTISHIFIPAMGNITYDFSFKSSVITLNTNAELFNQSDIVAHLLSSSSSVIDALQYKLEGTTRLTRKRGLKLATALSLSNKFVEGSHNSTVSLTTKNMEVSVATTTKAQIPILRMNFKQELNGNTKSKPTVSSSMEFKYDFNSSMLYSTAKGAVDHKLSLESLTSYFSIESSTKGDVKGSVLSREYSGTIASEANTYLNSKSTRSSVKLQGTSKIDDIWNLEVKENFAGEATLQRIYSLWEHSTKNHLQLEGLFFTNGEHTSKATLELSPWQMSALVQVHASQPSSFHDFPDLGQEVALNANTKNQKIRWKNEVRIHSGSFQSQVELSNDQEKAHLDIAGSLEGHLRFLKNIILPVYDKSLWDFLKLDVTTSIGRRQHLRVSTAFVYTKNPNGYSFSIPVKVLADKFIIPGLKLNDLNSVLVMPTFHVPFTDLQVPSCKLDFREIQIYKKLRTSSFALNLPTLPEVKFPEVDVLTKYSQPEDSLIPFFEITVPESQLTVSQFTLPKSVSDGIAALDLNAVANKIADFELPTIIVPEQTIEIPSIKFSVPAGIVIPSFQALTARFEVDSPVYNATWSASLKNKADYVETVLDSTCSSTVQFLEYELNVLGTHKIEDGTLASKTKGTFAHRDFSAEYEEDGKYEGLQEWEGKAHLNIKSPAFTDLHLRYQKDKKGISTSAASPAVGTVGMDMDEDDDFSKWNFYYSPQSSPDKKLTIFKTELRVRESDEETQIKVNWEEEAASGLLTSLKDNVPKATGVLYDYVNKYHWEHTGLTLREVSSKLRRNLQNNAEWVYQGAIRQIDDIDVRFQKAASGTTGTYQEWKDKAQNLYQELLTQEGQASFQGLKDNVFDGLVRVTQEFHMKVKHLIDSLIDFLNFPRFQFPGKPGIYTREELCTMFIREVGTVLSQVYSKVHNGSEILFSYFQDLVITLPFELRKHKLIDV",
        precomputed_pdb_path = "cycle2_alphafold2_output.pdb" # 

        if cycle == "2A":
            contigs="A100-140/0 15-25"
        if cycle == "2B":
            contigs="A150-190/0 15-25"
        if cycle == "2C":
            contigs="A200-240/0 15-25"
        if cycle == "2D":
            contigs="A250-290/0 15-25"
    else:
        raise ValueError("Invalid cycle number.")

# Set up variables part 2
name = f"cycle{cycle}_{num_seq}seqs_{diffusion}diff_{temp}temp"
outdir = f"{root}/cycle{cycle}_{diffusion}diff_{temp}temp"
os.makedirs(outdir, exist_ok=True)
print(f"Output dir : {outdir}")

# Check if AlphaFold2 PDB exists
if not os.path.exists(precomputed_pdb_path):
    raise FileNotFoundError(f"Precomputed PDB file {precomputed_pdb_path} does not exist. Please check the path.")
print()

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
    AF2_MULTIMER_PORT = 8084
class NIM_ENDPOINTS(StrEnum):
    RFDIFFUSION =  "biology/ipd/rfdiffusion/generate"
    PROTEINMPNN =  "biology/ipd/proteinmpnn/predict"
    AF2_MULTIMER = "protein-structure/alphafold2/multimer/predict-structure-from-sequences"

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
    contigs=contigs,  # Region A400-A600, peptides 15-25 residues long
    hotspot_res=[], 
    input_pdb_chains=[], # [Optional] default is to design for all chains in the protein
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
    # "hotspot_res": example.hotspot_res
}
rc, rfdiffusion_response = query_nim(
    payload=rfdiffusion_query,
    nim_endpoint=NIM_ENDPOINTS.RFDIFFUSION.value,
    nim_port=NIM_PORTS.RFDIFFUSION_PORT.value
)

# save
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

# Save binder_target_pairs to a .txt file
with open(f"{outdir}/3_{name}_proteinmpnn.txt", "w") as txt_file:
    for i, pair in enumerate(binder_target_pairs):
        binder, target = pair
        txt_file.write(f"Binder {i+1}: {binder}\n")
        txt_file.write(f"Target {i+1}: {target}\n")
        txt_file.write("\n")  # Add a blank line between pairs for readability

with open(f"{outdir}/3_{name}_proteinmpnn_pairs.json", "w") as json_file:
    json.dump(binder_target_pairs, json_file, indent=4)

# Save proteinmpnn_response["mfasta"] to a .fasta file
with open(f"{outdir}/3_{name}_proteinmpnn.fasta", "w") as fasta_file:
    fasta_file.write(proteinmpnn_response["mfasta"])

# Save scores and probs to files
scores = proteinmpnn_response["scores"]
with open(f"{outdir}/3_{name}_proteinmpnn_scores.txt", "w") as scores_file:
    for i, score in enumerate(scores):
        scores_file.write(f"Sequence {i+1}: Score = {score}\n")

# probs = proteinmpnn_response["probs"]
# with open(f"{outdir}/3_{name}_proteinmpnn_probs.txt", "w") as probs_file:
#     for i, prob_matrix in enumerate(probs):
#         probs_file.write(f"Sequence {i+1}:\n")
#         for position_probs in prob_matrix:
#             probs_file.write(",".join(map(str, position_probs)) + "\n")
#         probs_file.write("\n")

##############################################################
# 4. AlphaFold2-Multimer
##############################################################

# print(f"Loading AlphaFold-Multimer...")
# print()

# # Load binder_target_pairs from the JSON file
# binder_target_pairs_path = f"{outdir}/3_{name}_binder_target_pairs.json"
# with open(binder_target_pairs_path, "r") as json_file:
#     binder_target_pairs = json.load(json_file)

# # print preview of binder_target_pairs
# print(binder_target_pairs[:2])  # Print the first 2 pairs for preview


# n_processed = 0
# multimer_response_codes = [0 for i in binder_target_pairs]
# multimer_results = [None for i in binder_target_pairs]

# # change this value to process more or fewer target-binder pairs.
# pairs_to_process = 1

# for binder_target_pair in binder_target_pairs:
#     multimer_query = {
#         "sequences" : binder_target_pair,
#         "selected_models" : [1]
#     }
#     print(f"Processing pair number {n_processed+1} of {len(binder_target_pairs)}")
#     rc, multimer_response = query_nim(
#         payload=multimer_query,
#         nim_endpoint=NIM_ENDPOINTS.AF2_MULTIMER.value,
#         nim_port=NIM_PORTS.AF2_MULTIMER_PORT.value
#     )
#     multimer_response_codes[n_processed] = rc
#     multimer_results[n_processed] = multimer_response
#     print(f"Finished binder-target pair number {n_processed+1} of {len(binder_target_pairs)}")
#     n_processed += 1
#     if n_processed >= pairs_to_process:
#         break

# ## Print just the first 160 characters of the first multimer response
# result_idx = 0
# prediction_idx = 0
# print(multimer_results[result_idx][prediction_idx][0:160])
# print()

# # Save all AlphaFold-Multimer results to a .txt file
# with open(f"{root}/4_multimer_{name}.txt", "w") as results_file:
#     for i, result in enumerate(multimer_results):
#         results_file.write(f"Result {i+1}:\n")  # Add a header for each result
#         if result is not None:  # Check if the result exists
#             for prediction_idx, prediction in enumerate(result):
#                 results_file.write(f"Prediction {prediction_idx+1}:\n")
#                 results_file.write(prediction)
#                 results_file.write("\n\n")  # Add spacing between predictions
#         else:
#             results_file.write("No result available for this pair.\n")
#         results_file.write("-" * 50 + "\n")  # Separator between results

