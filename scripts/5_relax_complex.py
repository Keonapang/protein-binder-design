######################################################################
# Usage: python 5_relax_complex.py peptide_merged.pdb complex_minimized.pdb --solvent

# If too slow, we can use a faster (but not as good option):
######################################################################

import sys
from openmm import app
import openmm as mm
from openmm import unit

if len(sys.argv) < 3:
    print("Usage: python relax_complex.py peptide_merged.pdb complex_minimized.pdb [--solvent] [--time_ns N]")
    sys.exit(1)

input_pdb = sys.argv[1]
output_pdb = sys.argv[2]

use_solvent = "--solvent" in sys.argv
time_ns = 1.0  # default 1 ns
for arg in sys.argv:
    if arg.startswith("--time_ns"):
        try:
            time_ns = float(arg.split()[1])  # --time_ns 0.5
        except Exception:
            pass

# Load structure
pdb = app.PDBFile(input_pdb)
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

if not use_solvent:
    print(">> Running vacuum minimization...")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None
    )
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy(maxIterations=1000)

else:
    print(">> Running explicit solvent relaxation...")
    # Create solvated system in a water box
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(forcefield, model='tip3p', padding=1.0*unit.nanometer, ionicStrength=0.15*unit.molar)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometer,
        constraints=app.HBonds
    )

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
    simulation.reporters.append(app.StateDataReporter(sys.stdout, 50000, step=True,
                                                     potentialEnergy=True, temperature=True, progress=True,
                                                     speed=True, totalSteps=nsteps, separator='\t'))
    simulation.step(nsteps)

    # Take final coordinates
    positions = simulation.context.getState(getPositions=True).getPositions()
    pdb = app.PDBFile.writeFile(modeller.topology, positions, open(output_pdb, 'w'))
    print(f">> Saved solvated + relaxed complex to {output_pdb}")
    sys.exit(0)

# For vacuum run, save minimized coords
positions = simulation.context.getState(getPositions=True).getPositions()
with open(output_pdb, "w") as f:
    app.PDBFile.writeFile(pdb.topology, positions, f)
print(f">> Saved minimized complex to {output_pdb}")
