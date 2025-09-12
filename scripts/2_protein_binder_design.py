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
# python SCRIPT_protein_binder_design.py --cycle 1 --num_seq 5 --diffusion 30 --temp 0.2

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
# num_seq = 5
# diffusion = 25
# temp = 0.2

# Set up all paths and variables
root = "/home/ubuntu/nvidia-workbench"
os.makedirs(root, exist_ok=True)
print(f"Generating {num_seq} sequences per target for cycle {cycle}...")

if num_seq > 0: # if num_seq > 1, then run code below
    if cycle == "1":
        contigs="15-25" # contigs = "A400-700/15-25"???
        precomputed_pdb_path ="/home/ubuntu/cycle1_alphafold2_output.pdb" # converts pdb_id to a Path Object;
        target_sequence="MDPPRPALLALLALPALLLLLLAGARAEEEMLENVSLVCPKDATRFKHLRKYTYNYEAESSSGVPGTADSRSATRINCKVELEVPQLCSFILKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYELKLAIPEGKQVFLYPEKDEPTYILNIKRGIISALLVPPETEEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLGQCDRFKPIRTGISPLALIKGMTRPLSTLISSSQSCQYTLDAKRKHVAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKINSRFFGEGTKKMGLAFESTKSTSPPKQAEAVLKTLQELKKLTISEQNIQRANLFNKLVTELRGLSDEAVTSLLPQLIEVSSPITLQALVQCGQPQCSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREIFNMARDQRSRATLYALSHAVNNYHKTNPTGTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQTMEQLTPELKSSILKCVQSTKPSLMIQKAAIQALRKMEPKDKDQEVLLQTFLDDASPGDKRLAAYLMLMRSPSQADINKIVQILPWEQNEQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPTVMDFRKFSRNYQLYKSVSLPSLDPASAKIEGNLIFDPNNYLPKESMLKTTLTAFGFASADLIEIGLEGKGFEPTLEALFGKQGFFPDSVNKALYWVNGQVPDGVSKVLVDHFGYTKDDKHEQDMVNGIMLSVEKLIKDLKSKEVPEARAYLRILGEELGFASLHDLQLLGKLLLMGARTLQGIPQMIGEVIRKGSKNDFFLHYIFMENAFELPTGAGLQLQISSSGVIAPGAKAGVKLEVANMQAELVAKPSVSVEFVTNMGIIIPDFARSGVQMNTNFFHESGLEAHVALKAGKLKFIIPSPKRPVKLLSGGNTLHLVSTTKTEVIPPLIENRQSWSVCKQVFPGLNYCTSGAYSNASSTDSASYYPLTGDTRLELELRPTGEIEQYSVSATYELQREDRALVDTLKFVTQAEGAKQTEATMTFKYNRQSMTLSSEVQIPDFDVDLGTILRVNDESTEGKTSYRLTLDIQNKKITEVALMGHLSCDTKEERKIKGVISIPRLQAEARSEILAHWSPAKLLLQMDSSATAYGSTVSKRVAWHYDEEKIEFEWNTGTNVDTKKMTSNFPVDLSDYPKSLHMYANRLLDHRVPQTDMTFRHVGSKLIVAMSSWLQKASGSLPYTQTLQDHLNSLKEFNLQNMGLPDFHIPENLFLKSDGRVKYTLNKN"
    elif cycle == "2":
        contigs="15-25"
        target_sequence="GKIDFLNNYALFLSPSAQQASWQVSARFNQYKYNQNFSAGNNENIMEAHVGINGEANLDFLNIPLTIPEMRLPYTIITTPPLKDFSLWEKTGLKEFLKTTKQSFDLSVKAQYKKNKHRHSITNPLAVLCEFISQSIKSFDRHFEKNRNNALDFVTKSYNETKIKFDKYKAEKSHDELPRTFQIPGYTVPVVNVEVSPFTIEMSAFGYVFPKAVSMPSFSILGSDVRVPSYTLILPSLELPVLHVPRNLKLSLPDFKELCTISHIFIPAMGNITYDFSFKSSVITLNTNAELFNQSDIVAHLLSSSSSVIDALQYKLEGTTRLTRKRGLKLATALSLSNKFVEGSHNSTVSLTTKNMEVSVATTTKAQIPILRMNFKQELNGNTKSKPTVSSSMEFKYDFNSSMLYSTAKGAVDHKLSLESLTSYFSIESSTKGDVKGSVLSREYSGTIASEANTYLNSKSTRSSVKLQGTSKIDDIWNLEVKENFAGEATLQRIYSLWEHSTKNHLQLEGLFFTNGEHTSKATLELSPWQMSALVQVHASQPSSFHDFPDLGQEVALNANTKNQKIRWKNEVRIHSGSFQSQVELSNDQEKAHLDIAGSLEGHLRFLKNIILPVYDKSLWDFLKLDVTTSIGRRQHLRVSTAFVYTKNPNGYSFSIPVKVLADKFIIPGLKLNDLNSVLVMPTFHVPFTDLQVPSCKLDFREIQIYKKLRTSSFALNLPTLPEVKFPEVDVLTKYSQPEDSLIPFFEITVPESQLTVSQFTLPKSVSDGIAALDLNAVANKIADFELPTIIVPEQTIEIPSIKFSVPAGIVIPSFQALTARFEVDSPVYNATWSASLKNKADYVETVLDSTCSSTVQFLEYELNVLGTHKIEDGTLASKTKGTFAHRDFSAEYEEDGKYEGLQEWEGKAHLNIKSPAFTDLHLRYQKDKKGISTSAASPAVGTVGMDMDEDDDFSKWNFYYSPQSSPDKKLTIFKTELRVRESDEETQIKVNWEEEAASGLLTSLKDNVPKATGVLYDYVNKYHWEHTGLTLREVSSKLRRNLQNNAEWVYQGAIRQIDDIDVRFQKAASGTTGTYQEWKDKAQNLYQELLTQEGQASFQGLKDNVFDGLVRVTQEFHMKVKHLIDSLIDFLNFPRFQFPGKPGIYTREELCTMFIREVGTVLSQVYSKVHNGSEILFSYFQDLVITLPFELRKHKLIDV"
        precomputed_pdb_path = "/home/ubuntu/cycle2_alphafold2_output.pdb" # 

    elif "1" in cycle:
        target_sequence="MDPPRPALLALLALPALLLLLLAGARAEEEMLENVSLVCPKDATRFKHLRKYTYNYEAESSSGVPGTADSRSATRINCKVELEVPQLCSFILKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYELKLAIPEGKQVFLYPEKDEPTYILNIKRGIISALLVPPETEEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLGQCDRFKPIRTGISPLALIKGMTRPLSTLISSSQSCQYTLDAKRKHVAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKINSRFFGEGTKKMGLAFESTKSTSPPKQAEAVLKTLQELKKLTISEQNIQRANLFNKLVTELRGLSDEAVTSLLPQLIEVSSPITLQALVQCGQPQCSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREIFNMARDQRSRATLYALSHAVNNYHKTNPTGTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQTMEQLTPELKSSILKCVQSTKPSLMIQKAAIQALRKMEPKDKDQEVLLQTFLDDASPGDKRLAAYLMLMRSPSQADINKIVQILPWEQNEQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPTVMDFRKFSRNYQLYKSVSLPSLDPASAKIEGNLIFDPNNYLPKESMLKTTLTAFGFASADLIEIGLEGKGFEPTLEALFGKQGFFPDSVNKALYWVNGQVPDGVSKVLVDHFGYTKDDKHEQDMVNGIMLSVEKLIKDLKSKEVPEARAYLRILGEELGFASLHDLQLLGKLLLMGARTLQGIPQMIGEVIRKGSKNDFFLHYIFMENAFELPTGAGLQLQISSSGVIAPGAKAGVKLEVANMQAELVAKPSVSVEFVTNMGIIIPDFARSGVQMNTNFFHESGLEAHVALKAGKLKFIIPSPKRPVKLLSGGNTLHLVSTTKTEVIPPLIENRQSWSVCKQVFPGLNYCTSGAYSNASSTDSASYYPLTGDTRLELELRPTGEIEQYSVSATYELQREDRALVDTLKFVTQAEGAKQTEATMTFKYNRQSMTLSSEVQIPDFDVDLGTILRVNDESTEGKTSYRLTLDIQNKKITEVALMGHLSCDTKEERKIKGVISIPRLQAEARSEILAHWSPAKLLLQMDSSATAYGSTVSKRVAWHYDEEKIEFEWNTGTNVDTKKMTSNFPVDLSDYPKSLHMYANRLLDHRVPQTDMTFRHVGSKLIVAMSSWLQKASGSLPYTQTLQDHLNSLKEFNLQNMGLPDFHIPENLFLKSDGRVKYTLNKN",
        precomputed_pdb_path = "/home/ubuntu/cycle1_alphafold2_output.pdb" # 
        
        # if cycle == "1A": 
        #     contigs="A400-430/15-25"
        # if cycle == "1B":
        #     contigs="A450-480/15-25"
        # if cycle == "1C":
        #     contigs="A500-540/15-25"
        # if cycle == "1D":
        #     contigs="A560-600/15-25"
        if cycle == "1A": 
            contigs="430-450/15-25"
        if cycle == "1B":
            contigs="A480-510/15-25"
        if cycle == "1C":
            contigs="A520-545/15-25"
        if cycle == "1D":
            contigs="A575-600/15-25"
            
    elif "2" in cycle:
        target_sequence="GKIDFLNNYALFLSPSAQQASWQVSARFNQYKYNQNFSAGNNENIMEAHVGINGEANLDFLNIPLTIPEMRLPYTIITTPPLKDFSLWEKTGLKEFLKTTKQSFDLSVKAQYKKNKHRHSITNPLAVLCEFISQSIKSFDRHFEKNRNNALDFVTKSYNETKIKFDKYKAEKSHDELPRTFQIPGYTVPVVNVEVSPFTIEMSAFGYVFPKAVSMPSFSILGSDVRVPSYTLILPSLELPVLHVPRNLKLSLPDFKELCTISHIFIPAMGNITYDFSFKSSVITLNTNAELFNQSDIVAHLLSSSSSVIDALQYKLEGTTRLTRKRGLKLATALSLSNKFVEGSHNSTVSLTTKNMEVSVATTTKAQIPILRMNFKQELNGNTKSKPTVSSSMEFKYDFNSSMLYSTAKGAVDHKLSLESLTSYFSIESSTKGDVKGSVLSREYSGTIASEANTYLNSKSTRSSVKLQGTSKIDDIWNLEVKENFAGEATLQRIYSLWEHSTKNHLQLEGLFFTNGEHTSKATLELSPWQMSALVQVHASQPSSFHDFPDLGQEVALNANTKNQKIRWKNEVRIHSGSFQSQVELSNDQEKAHLDIAGSLEGHLRFLKNIILPVYDKSLWDFLKLDVTTSIGRRQHLRVSTAFVYTKNPNGYSFSIPVKVLADKFIIPGLKLNDLNSVLVMPTFHVPFTDLQVPSCKLDFREIQIYKKLRTSSFALNLPTLPEVKFPEVDVLTKYSQPEDSLIPFFEITVPESQLTVSQFTLPKSVSDGIAALDLNAVANKIADFELPTIIVPEQTIEIPSIKFSVPAGIVIPSFQALTARFEVDSPVYNATWSASLKNKADYVETVLDSTCSSTVQFLEYELNVLGTHKIEDGTLASKTKGTFAHRDFSAEYEEDGKYEGLQEWEGKAHLNIKSPAFTDLHLRYQKDKKGISTSAASPAVGTVGMDMDEDDDFSKWNFYYSPQSSPDKKLTIFKTELRVRESDEETQIKVNWEEEAASGLLTSLKDNVPKATGVLYDYVNKYHWEHTGLTLREVSSKLRRNLQNNAEWVYQGAIRQIDDIDVRFQKAASGTTGTYQEWKDKAQNLYQELLTQEGQASFQGLKDNVFDGLVRVTQEFHMKVKHLIDSLIDFLNFPRFQFPGKPGIYTREELCTMFIREVGTVLSQVYSKVHNGSEILFSYFQDLVITLPFELRKHKLIDV",
        precomputed_pdb_path = "/home/ubuntu/cycle2_alphafold2_output.pdb" # 

        if cycle == "2A":
            contigs="A110-135/15-25"
        if cycle == "2B":
            contigs="A150-175/15-25"
        if cycle == "2C":
            contigs="A200-230/15-25"
        if cycle == "2D":
            contigs="A250-275/15-25"
    else:
        raise ValueError("Invalid cycle number.")


# Set up variables part 2
name = f"cycle{cycle}_{num_seq}seqs_{diffusion}diff_{temp}temp"
outdir = f"{root}/{diffusion}diff_{temp}temp_{num_seq}numseq"
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

# Save scores and probs to files
# scores = proteinmpnn_response["scores"]
# with open(f"{outdir}/3_{name}_proteinmpnn_scores.txt", "w") as scores_file:
#     for i, score in enumerate(scores):
#         scores_file.write(f"Sequence {i+1}: Score = {score}\n")

probs = proteinmpnn_response["probs"]
with open(f"{outdir}/3_{name}_proteinmpnn_probs.txt", "w") as probs_file:
    for i, prob_matrix in enumerate(probs):
        probs_file.write(f"Sequence {i+1}:\n")
        for position_probs in prob_matrix:
            probs_file.write(",".join(map(str, position_probs)) + "\n")
        probs_file.write("\n")

##############################################################
# 4. AlphaFold2-Multimer
# binder_target_pair contains two sequences: [binder_sequence, target_sequence]
# can predict interactions between multiple chains (e.g., protein-protein complexes or protein-peptide complexes).
# AlphaFold2-Multimer should attempt to predict the structure of the entire complex, 
# # including all 5 components (target + 4 peptides).
##############################################################
# print()
# print(f"Loading AlphaFold-Multimer...")
# print()

# # Load binder_target_pairs from the JSON file
# # os.chdir(outdir)
# binder_target_pairs_path = f"3_{name}_proteinmpnn_pairs.json" # 3_cycle1A_5seqs_25diff_0.1temp_proteinmpnn_pairs.json
# with open(binder_target_pairs_path, "r") as json_file:
#     binder_target_pairs = json.load(json_file)
# print(binder_target_pairs[:2])  # Print the first 2 pairs for preview

# # Intiate variables
# n_processed = 0 # Tracks how many binder-target pairs have been processed.
# multimer_response_codes = [0 for i in binder_target_pairs] # A list to store response codes (rc) from each query to the AlphaFold2-Multimer API. A value of 0 typically indicates success
# multimer_results = [None for i in binder_target_pairs] # A list to store the results from each pair

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
# os.chdir(outdir)# change directory to outdir
# with open(f"{outdir}/4_{name}_multimer.txt", "w") as results_file:
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


# # Save PDB structures for each prediction
# for i, result in enumerate(multimer_results):
#     if result is not None:
#         for prediction_idx, prediction in enumerate(result):
#             pdb_file_path = f"{root}/4_{name}_{i+1}_prediction{prediction_idx+1}.pdb"
#             with open(pdb_file_path, "w") as pdb_file:
#                 pdb_file.write(prediction)  # Assuming prediction is in PDB format

# ####################################################################
# # PART 5: VALIDATION
# ####################################################################
# print()
# print(f"Loading validation analyses...")
# print()

# multimer_path = f"4_multimer.json" # 3_cycle1A_5seqs_25diff_0.1temp_proteinmpnn_pairs.json
# with open(multimer_path, "r") as json_file:
#     multimer_results = json.load(json_file)
# print(multimer_results[:2])  # Print the first 2 pairs for preview


# # Function to calculate average pLDDT over all residues 
# def calculate_average_pLDDT(pdb_string):
#     total_pLDDT = 0.0
#     atom_count = 0
#     pdb_lines = pdb_string.splitlines()
#     for line in pdb_lines:
#         # PDB atom records start with "ATOM"
#         if line.startswith("ATOM"):
#             atom_name = line[12:16].strip() # Extract atom name
#             if atom_name == "CA":  # Only consider atoms with name "CA"
#                 try:
#                     # Extract the B-factor value from columns 61-66 (following PDB format specifications)
#                     pLDDT = float(line[60:66].strip())
#                     total_pLDDT += pLDDT
#                     atom_count += 1
#                 except ValueError:
#                     pass  # Skip lines where B-factor can't be parsed as a float

#     if atom_count == 0:
#         return 0.0  # Return 0 if no N atoms were found

#     average_pLDDT = total_pLDDT / atom_count
#     return average_pLDDT

# # Run 
# plddts = []
# for idx in range(0, len(multimer_results)):
#     if multimer_results[idx] is not None:
#         plddts.append(calculate_average_pLDDT(multimer_results[idx][0]))

# ## Combine the results with their pLDDTs
# binder_target_results = list(zip(binder_target_pairs, multimer_results, plddts))

# ## Sort the results by plddt
# sorted_binder_target_results = sorted(binder_target_results, key=lambda x : x[2])

# ## print the top 5 results
# for i in range(0, len(sorted_binder_target_results)):
#     print("-"*80)
#     print(f"rank: {i}")
#     print(f"binder: {sorted_binder_target_results[i][0][0]}")
#     print(f"target: {sorted_binder_target_results[i][0][1]}")
#     print(f"pLDDT: {sorted_binder_target_results[i][2]}")
#     print("-"*80)

# # Save all results to a .txt file
# with open(f"{outdir}/5_{name}_assess.txt", "w") as results_file:
#     # Section 1: Average pLDDT values
#     results_file.write("### Average pLDDT Values ###\n")
#     for idx, plddt in enumerate(plddts):
#         results_file.write(f"Result {idx + 1}: Average pLDDT = {plddt:.2f}\n")
#     results_file.write("\n\n")

#     # Section 2: Binder-Target Results with Raw Data
#     results_file.write("### Binder-Target Results (Raw Data) ###\n")
#     for idx, (pair, result, plddt) in enumerate(binder_target_results):
#         results_file.write(f"Result {idx + 1}:\n")
#         results_file.write(f"Binder: {pair[0]}\n")
#         results_file.write(f"Target: {pair[1]}\n")
#         results_file.write(f"Average pLDDT: {plddt:.2f}\n")
#         results_file.write(f"Raw Multimer Result: {result}\n")
#         results_file.write("-" * 80 + "\n")
#     results_file.write("\n\n")

#     # Section 3: Sorted Binder-Target Results by pLDDT
#     results_file.write("### Sorted Binder-Target Results by pLDDT ###\n")
#     for rank, (pair, result, plddt) in enumerate(sorted_binder_target_results):
#         results_file.write(f"Rank {rank + 1}:\n")
#         results_file.write(f"Binder: {pair[0]}\n")
#         results_file.write(f"Target: {pair[1]}\n")
#         results_file.write(f"Average pLDDT: {plddt:.2f}\n")
#         results_file.write("-" * 80 + "\n")

# print(f"------------------------------------------------")
# print()