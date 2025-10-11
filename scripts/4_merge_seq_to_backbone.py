######################################################################
# Align the peptide binder onto the target protein and create a combined PDB file

# Inputs:
# 1. Backbone PDB (from RFdiffusion, containing the peptide chain with placeholder residues, usually “ALA” or poly-GLY).
# 2. ProteinMPNN FASTA file (with the designed peptide sequence).
# 3. Outputs a new PDB where the peptide residues are replaced with the correct sequence.

# Usage: python 4_merge_seq_to_backbone.py ${REPO_DIR} B peptide.fasta peptide_merged.pdb

# Where:

# rf_backbone.pdb → RFdiffusion output structure
# B → chain ID of the peptide
# peptide.fasta → ProteinMPNN-designed sequence
# peptide_merged.pdb → final PDB with correct peptide sequence
######################################################################
import sys
import time
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO
import openmm as mm
from openmm import app, unit
from openmm.app import PDBFile
from pdbfixer import PDBFixer

REPO_DIR = sys.argv[1]
chain_id = sys.argv[2]   # e.g. "B"
iterations = sys.argv[3]            # Replace with the number of iterations
num_seq = sys.argv[4]                 # Replace with the number of sequences per iteration
name = sys.argv[5]  # e.g. target_A10_20
params = sys.argv[6]  # e.g. 1diff_0.5
use_solvent = "--solvent" in sys.argv
time_ns = 1.0  # default 1 ns

num_seq = int(num_seq)
iterations = int(iterations)

# REPO_DIR = "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/protein-binder-design"
# REPO_DIR = "/home/shadeform/protein-binder-design"
# chain_id="A"
# iterations=1
# num_seq=1
# name="target_A60_90"
# params="50diff_0.3temp"

if len(sys.argv) < 6:
    print("Usage: python 4_merge_seq_to_backbone.py ${REPO_DIR} B peptide.fasta peptide_merged.pdb iterations num_seq params --solvent")

# Map 1-letter → 3-letter codesd
aa_map = {
    'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE',
    'G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU',
    'M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG',
    'S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'
}

import os
os.chdir(REPO_DIR)

num_structures=iterations * num_seq

print(f"\n########################################")
print(f"4. Merging sequences onto backbone")
print(f"########################################\n")
print("Current directory:", os.getcwd())
print(f"Total structures to process: {num_structures}")
total_start_time = time.time()
count=1

if not os.path.exists(f"{REPO_DIR}/{name}/2_{name}_{params}_i1.pdb"):
    raise FileNotFoundError(f"File not found: {REPO_DIR}/{name}/2_{name}_{params}_i1.pdb")
if not os.path.exists(f"{REPO_DIR}/{name}/3_{name}_{params}_i1.fasta"):
    raise FileNotFoundError(f"File not found: {REPO_DIR}/{name}/3_{name}_{params}_i1.fasta")

for iteration in range(1, iterations + 1):
    for num in range(1, num_seq + 1):

        print(f"\n========== Structure {count} of {num_structures} ==========\n")
        count+=1

        # 1. Align ProteinPMNN's peptide sequence to RFdiffusion's backbone structure
        ###############################################################################

        # Input RFDiffusion and proteinMPNN files
        backbone_pdb=f"{REPO_DIR}/{name}/2_{name}_{params}_i{iteration}.pdb"
        pmnn_file=f"{REPO_DIR}/{name}/3_{name}_{params}_i{iteration}.fasta"

        # Output PDB file
        output_pdb=f"{REPO_DIR}/{name}/4_{name}_{params}_i{iteration}_{num}_temp.pdb"

        try:
            with open(pmnn_file, "r") as file:
                lines = file.readlines()
                line_number = (num * 2) + 2 - 1  # Adjust to 0-based index
                new_seq = lines[line_number].strip()  # Get the sequence and strip whitespace
                print(f"ProteinMPNN sequence: {new_seq}")
        except FileNotFoundError:
            print(f"File not found: {pmnn_file}")
        except IndexError:
            print(f"Invalid line range for Sequence {num} in file: {pmnn_file}")
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("model", backbone_pdb)
        # Replace residue names in the peptide chain
        for model in structure:
            for chain in model: # chain corresponds to a sequence of residues with the same chain ID
                if chain.id == chain_id:
                    res_list = [res for res in chain if res.id[0] == " "]
                    print(f"Merging seq to backbone. Chain {chain_id} has {len(res_list)} residues.")
                    if len(res_list) != len(new_seq): # New sequence (15 residues)
                        raise ValueError(f"\nProteinPMNN seq length {len(new_seq)} doesn't match number of residues {len(res_list)} in chain {chain_id}\n")
                    for res, aa in zip(res_list, new_seq):
                        res.resname = aa_map[aa]
        # Write new PDB
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)

        # 2. Fix PDB file by adding missing side chains and hydrogens 
        #######################################################################

        if not os.path.exists(output_pdb):
            raise FileNotFoundError(f"Output PDB file not found: {output_pdb}")
        fixer = PDBFixer(filename=f"{output_pdb}")

        # output fixed file
        output_fixed_pdb=f"{REPO_DIR}/{name}/4_{name}_{params}_i{iteration}_{num}.pdb"

        # Add missing side-chain atoms
        fixer.findMissingResidues()  # Looks for missing residues
        fixer.findMissingAtoms()     # Detects missing atoms
        fixer.addMissingAtoms()      # Adds missing atoms, including side chains
        fixer.addMissingHydrogens(7.0)  # Adds hydrogens at pH 7.0

        with open(output_fixed_pdb, 'w') as output:
            PDBFile.writeFile(fixer.topology, fixer.positions, output)

        with open(output_fixed_pdb, "r") as infile:
            lines = infile.readlines()  # Read all lines from the file
        with open(output_fixed_pdb, "w") as outfile:
            for line in lines:
                if not line.startswith("REMARK"):  # Skip lines starting with "REMARK"
                    outfile.write(line)  # Write the line, but don't implicitly print anything
        print(f"Merged PDB:{output_fixed_pdb}")

        if os.path.exists(output_pdb):
            os.remove(output_pdb)

        # 3. Performs energy minimization and molecular dynamics (MD) simulations
        #######################################################################

        start_time = time.time()
        for arg in sys.argv:
            if arg.startswith("--time_ns"):
                try:
                    time_ns = float(arg.split()[1])  # --time_ns 0.5
                except Exception:
                    pass

        # Load structure
        pdb = app.PDBFile(output_fixed_pdb) # from step 2
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

        optimized_pdb=f"{REPO_DIR}/{name}/5_{name}_{params}_i{iteration}_{num}_complex.pdb"

        if not use_solvent:
            print("\n>> Running vacuum minimization...")
            system = forcefield.createSystem(pdb.topology,nonbondedMethod=app.NoCutoff,constraints=None)
            integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
            simulation = app.Simulation(pdb.topology, system, integrator)
            simulation.context.setPositions(pdb.positions)
            simulation.minimizeEnergy(maxIterations=1000)
        else:
            print("\n>> Running explicit solvent relaxation...")
            # Create solvated system in a water box
            modeller = app.Modeller(pdb.topology, pdb.positions)
            modeller.addSolvent(forcefield, model='tip3p', padding=1.0*unit.nanometer, ionicStrength=0.15*unit.molar)
            system = forcefield.createSystem(modeller.topology,nonbondedMethod=app.PME,nonbondedCutoff=1.0*unit.nanometer,constraints=app.HBonds)
            integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
            simulation = app.Simulation(modeller.topology, system, integrator)
            simulation.context.setPositions(modeller.positions)
            # Minimization
            print("  > Minimizing...")
            simulation.minimizeEnergy(maxIterations=1000)
            # Equilibration (100 ps NVT + 100 ps NPT)
            simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
            print("  > Equilibrating (200 ps)...")
            simulation.step(int(200000/2))  # 200 ps, dt=2 fs
            # Short MD production run
            nsteps = int((time_ns*1000) / 0.002)  # ns → steps (dt=2 fs)
            print(f"  > Running production MD: {time_ns} ns ({nsteps} steps)...")
            simulation.reporters.append(app.StateDataReporter(sys.stdout, 50000, 
                                        step=True,potentialEnergy=True, temperature=True, progress=True,
                                        speed=True, totalSteps=nsteps, separator='\t')
                                        )
            simulation.step(nsteps)
            # Take final coordinates
            positions = simulation.context.getState(getPositions=True).getPositions()
            # uses the topology of the solvated system, including water molecules and ions)
            pdb = app.PDBFile.writeFile(modeller.topology, positions, open(optimized_pdb, 'w'))
            print(f">> Saved solvated + relaxed complex to {optimized_pdb}\n")
            # sys.exit(0)

        # optimized_pdb=f"{REPO_DIR}/{name}/5_{name}_{params}_i{iteration}_{num}_complex_vacuum.pdb"

        # For vacuum run, save minimized coords
        # positions = simulation.context.getState(getPositions=True).getPositions()
        # with open(optimized_pdb, "w") as f:
        #     app.PDBFile.writeFile(pdb.topology, positions, f)
        # print(f">> Saved minimized complex to {optimized_pdb}")

        # with open(optimized_pdb, "r") as infile:
        #     lines = infile.readlines()  # Read all lines from the file

        # with open(optimized_pdb, "w") as outfile:
        #     for line in lines:
        #         # Skip lines starting with "REMARK" or "CRYST1"
        #         if not (line.startswith("REMARK") or line.startswith("CRYST1")):
        #             outfile.write(line)  # Write the line if it doesn't match the condition

        # print(f"Cleaned PDB file saved to {optimized_pdb}")

        # record end time
        end_time = time.time()
        total_time = end_time - start_time
        minutes = total_time / 60
        print(f"Total time: {minutes:.2f} mins")

# record end time
total_end_time = time.time()
total_time = total_end_time - total_start_time
minutes = total_time / 60
print(f"\nTotal time to process {num_structures} structures: {minutes:.2f} mins\n")