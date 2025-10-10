from Bio.PDB import PDBParser, PDBIO, Superimposer, Select
import os
import sys

######################################################################
# Use BioPython's Superimposer module to align the peptide binder 
# onto the target protein and create a combined PDB file

# Outputs the combined structure as combined.pdb, with:
#       - Chain A for the target protein.
#       - Chain B for the peptide binder.

# Usage: python modify_pdb_chain.py <target_file> <binder_file> <diffusion> <temp> <DIR_OUT>
######################################################################

# Parse command-line arguments
if len(sys.argv) != 7:
    print("Usage: python align_complex.py <target_pdb> <binder_pdb> <diffusion> <temp> <DIR_OUT> <raw>")
    sys.exit(1)

# Parse command-line arguments
target_pdb = sys.argv[1]
binder_pdb = sys.argv[2]
diffusion = sys.argv[3]
temp = sys.argv[4]
DIR_OUT = sys.argv[5]
raw = sys.argv[6]  # Not used for now

new_chain = "B"
cycle = os.path.splitext(os.path.basename(target_pdb))[0]

file_name = os.path.basename(binder_pdb)  # Extract the file name from the full path
file_part = "_".join(file_name.split("_")[1:])  # Remove everything before and including the first "_"
outname = f"5_{file_part}"
if raw == "raw":
    outname = outname.replace(".pdb", "_raw.pdb")  # Add "_raw" before the ".pdb"
    outfile = os.path.join(DIR_OUT, outname)
else:
    outfile = os.path.join(DIR_OUT, outname)
print(f"\noutfile:", outfile)

# define output file name
if not os.path.exists(DIR_OUT):
    os.makedirs(DIR_OUT)
output_pdb = os.path.join(DIR_OUT, outname)

# Parse the PDB files
parser = PDBParser(QUIET=True)
target_structure = parser.get_structure("target", target_pdb)
binder_structure = parser.get_structure("binder", binder_pdb)

# Step 1: Filter the target structure to keep only chain A
# =========================================================
filtered_structure = target_structure.copy()  # Copy the structure
filtered_model = filtered_structure[0]  # Select the first model

# Remove all chains except chain "A"
for chain in list(filtered_model.get_chains()):
    if chain.id != "A":  # Retain only chain A
        filtered_model.detach_child(chain.id)

# Now, `filtered_structure` contains only chain A from the target

# Step 2: Extract CA atoms for alignment
# =======================================
# Select chain A from both structures
target_chain = filtered_model["A"]  # Chain A from the filtered target
binder_chain = binder_structure[0]["A"]  # Chain A from the binder

# Select CA atoms for alignment
target_atoms = [atom for atom in target_chain.get_atoms() if atom.get_id() == "CA"]
binder_atoms = [atom for atom in binder_chain.get_atoms() if atom.get_id() == "CA"]

# Ensure the atom lists are of the same size by using only the overlapping residues
if len(target_atoms) != len(binder_atoms):
    min_len = min(len(target_atoms), len(binder_atoms))
    target_atoms = target_atoms[:min_len]
    binder_atoms = binder_atoms[:min_len]

# Align the binder to the filtered target structure
super_imposer = Superimposer()
super_imposer.set_atoms(target_atoms, binder_atoms)
super_imposer.apply(binder_chain.get_atoms())  # Apply the transformation to the binder

# Step 3: Save the combined structure
# =====================================
# Rename the binder chain to "B" to avoid conflicts
binder_chain.id = "B"

# Add the binder chain to the filtered target structure
filtered_model.add(binder_chain)

# Save the combined structure
io = PDBIO()
io.set_structure(filtered_structure)
io.save(output_pdb)
print(f"Combined PDB file saved as: {output_pdb}")

# ############################################################################
# # Align binder to target (using CA atoms of the first chain)
# target_atoms = []
# binder_atoms = []

# for target_chain, binder_chain in zip(target_model, binder_model):
#     # Only use residue alpha carbons (CA) to align
#     target_atoms += [atom for atom in target_chain.get_atoms() if atom.get_id() == "CA"]
#     binder_atoms += [atom for atom in binder_chain.get_atoms() if atom.get_id() == "CA"]

# # Perform the alignment
# super_imposer = Superimposer()
# super_imposer.set_atoms(target_atoms, binder_atoms)

# # Apply the transformation to the binder structure
# super_imposer.apply(binder_model.get_atoms())

# # Save the combined structure
# class ChainSelect(Select):
#     def __init__(self, chain_id):
#         self.chain_id = chain_id

#     def accept_chain(self, chain):
#         return chain.id == self.chain_id

# io = PDBIO()
# io.set_structure(target_structure)
# io.save(output_pdb, select=ChainSelect("A"))  # Save target chain as chain A

# binder_chain = list(binder_model.get_chains())[0]
# binder_chain.id = "B"  # Change binder chain to B
# target_structure[0].add(binder_chain)  # Add binder chain to target structure

# io.set_structure(target_structure)
# io.save(output_pdb)

# print(f"\nCombined PDB file: {output_pdb}\n")