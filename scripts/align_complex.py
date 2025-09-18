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
if len(sys.argv) != 6:
    print("Usage: python align_complex.py <target_pdb> <binder_pdb> <diffusion> <temp> <DIR_OUT>")
    sys.exit(1)

target_pdb = sys.argv[1]
binder_pdb = sys.argv[2]
diffusion = sys.argv[3]
temp = sys.argv[4]
DIR_OUT = sys.argv[5]


new_chain = "B"
cycle = os.path.splitext(os.path.basename(target_pdb))[0]

file_name = os.path.basename(binder_pdb)  # Extract the file name from the full path
file_part = "_".join(file_name.split("_")[1:])  # Remove everything before and including the first "_"
outname = f"5_{file_part}"

# define output file name
if not os.path.exists(DIR_OUT):
    os.makedirs(DIR_OUT)
output_pdb = os.path.join(DIR_OUT, outname)

print("==================== ", cycle, "======================")

# Input files
# target_pdb = "target.pdb"  # Target protein
# binder_pdb = "binder.pdb"  # Peptide binder
# output_pdb = "combined.pdb"  # Output combined file

# Parse the PDB files
parser = PDBParser(QUIET=True)
target_structure = parser.get_structure("target", target_pdb)
binder_structure = parser.get_structure("binder", binder_pdb)

# Select the first model and chain from each structure
target_model = target_structure[0]
binder_model = binder_structure[0]

# Select chains (assume chain A for both structures)
target_chain = target_model["A"]  # Replace "A" with the correct chain ID for the target
binder_chain = binder_model["A"]  # Replace "A" with the correct chain ID for the binder

# Select CA atoms for alignment
target_atoms = [atom for atom in target_chain.get_atoms() if atom.get_id() == "CA"]
binder_atoms = [atom for atom in binder_chain.get_atoms() if atom.get_id() == "CA"]

# Debug: Print number of CA atoms selected
print(f"Number of CA atoms in target: {len(target_atoms)}")
print(f"Number of CA atoms in binder: {len(binder_atoms)}")

# Ensure the atom lists are of the same size by using only the overlapping residues
if len(target_atoms) != len(binder_atoms):
    min_len = min(len(target_atoms), len(binder_atoms))
    target_atoms = target_atoms[:min_len]
    binder_atoms = binder_atoms[:min_len]

# Perform the alignment
super_imposer = Superimposer()
super_imposer.set_atoms(target_atoms, binder_atoms)  # Align binder to target
super_imposer.apply(binder_model.get_atoms())  # Apply transformation to the binder

# Save the combined structure
class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

io = PDBIO()
io.set_structure(target_structure)
io.save("target_only.pdb", select=ChainSelect("A"))  # Save target chain as chain A

# Rename binder chain to B and add it to the target structure
binder_chain.id = "B"
target_structure[0].add(binder_chain)

# Save the combined structure
io.set_structure(target_structure)
io.save(output_pdb)
print(f"\nCombined PDB file saved as: {output_pdb}")

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