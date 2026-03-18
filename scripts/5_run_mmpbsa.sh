#!/bin/bash
#!/usr/bin/env bash
######################################################
# g_mmpbsa calculates relative binding free energy using the MM-PBSA method 

# Inputs: aligned PDB complex, MM/PBSA settings file (.in)
# Output: 
# MMPBSA_results.dat → total ΔG_bind (kcal/mol)
# MMPBSA_energy_per_frame.csv → per-frame energies
# If idecomp=1, you’ll also get per-residue contributions for interface analysis.

# Usage: chmod +x 6_run_mmpbsa.sh
# ./6_run_mmpbsa.sh

# Set up Environment (all open source)
########################################

REPO_DIR=$1
name=$2
params=$3
i=$4
num=$5

# --------- example variables ---------
# REPO_DIR="/home/shadeform/protein-binder-design"
# name="target_A60_90"
# params="50diff_0.3temp"
# i="1"
# num="1"
# wc -l ${REPO_DIR}/${name}/5_${name}_${params}_i${i}_${num}_complex.pdb

# export REPO_DIR=/home/shadeform/protein-binder-design
# export name="target_A60_90"
# export params="50diff_0.3temp"
# export i=1
# export num=1
# export FF=amber99sb-ildn
# export WATER=tip3p
# export BOX_D=1.0
# export ION_CONC=0.15
# export MD_TIME_NS=25
# export FRAME_EVERY_PS=100
# export RECEPTOR_CHAIN=A
# export LIGAND_CHAIN=B
# ------------------------------------

set -euo pipefail

# ==== USER SETTINGS ====
FF=amber99sb-ildn       # or charmm36-mar2023
WATER=tip3p             # (for charmm36 use tip3p)
BOX_D=1.0               # nm padding
ION_CONC=0.15           # M
MD_TIME_NS=25            # production length (increase to 5–20 ns for robustness)
FRAME_EVERY_PS=100      # snapshot spacing for MM/PBSA
RECEPTOR_CHAIN=A
LIGAND_CHAIN=B

start_time=$(date +%s)

input_pdb=${REPO_DIR}/${name}/5_${name}_${params}_i${i}_${num}_complex.pdb

echo ""
echo "####################################################"
echo "Running g_mmpbsa to calculate binding free energy"
echo "####################################################"
echo "Input aligned, optimized file: ${input_pdb}"
echo ""
cd $REPO_DIR

if [ ! -d "${REPO_DIR}/mdp" ]; then
    echo "Error: ${REPO_DIR}/mdp folder not found!"
fi

# =================================================================

log() { echo -e "\n[$(date '+%H:%M:%S')] $*"; }

log "1) pdb2gmx: build topology"
# gmx pdb2gmx -f ${REPO_DIR}/${name}/5_target_A60_90_50diff_0.3temp_i1_1_complex.pdb -o complex_processed.gro -p topol.top -i posre.itp -ff $FF -water $WATER -ter <<EOF
gmx pdb2gmx -f ${input_pdb} -o complex_processed.gro -p topol.top -i posre.itp -ff $FF -water $WATER -ter <<EOF
1
1
EOF
# Above answers may prompt termini states; adjust as desired (e.g., capped vs charged peptide).

log "2) Define box"
gmx editconf -f complex_processed.gro -o boxed.gro -c -d $BOX_D -bt cubic

log "3) Solvate"
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

log "4) Add ions (prepare EM tpr)"
gmx grompp -f mdp/ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1
# neutralize + set ionic strength
echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc $ION_CONC

log "5) Energy minimization"
gmx grompp -f mdp/minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em

log "6) NVT equilibration"
gmx grompp -f mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

log "7) NPT equilibration"
gmx grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

log "8) Production MD (${MD_TIME_NS} ns)"
# Set nsteps in mdp/md.mdp to match MD_TIME_NS; or override via sed:
NSTEPS=$(( MD_TIME_NS * 1000000 / 2 )) # 2 fs timestep → steps per ns = 500k
sed "s/^nsteps.*/nsteps = ${NSTEPS}/" mdp/md.mdp > mdp/md_runtime.mdp
gmx grompp -f mdp/md_runtime.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md

log "9) Center, make trajectory, strip PBC for analysis"
# center on protein+peptide and output every FRAME_EVERY_PS
echo "Protein System" | gmx trjconv -s md.tpr -f md.xtc -o md_centered.xtc -pbc mol -center
gmx trjconv -s md.tpr -f md_centered.xtc -o md_strided.xtc -dt $FRAME_EVERY_PS <<EOF
0
EOF

log "10) Build index groups: Receptor (chain $RECEPTOR_CHAIN), Ligand (chain $LIGAND_CHAIN)"
echo -e "chain $RECEPTOR_CHAIN\nname 3 Receptor\nchain $LIGAND_CHAIN\nname 4 Ligand\nq" | gmx make_ndx -f md.gro -o index.ndx

log "11) MM/PBSA calculation (fast path via gmx_MMPBSA)"

# If you installed gmx_MMPBSA via pip:
gmx_MMPBSA -O \
  -i mmpbsa.in \
  -cs md.tpr \
  -ci index.ndx \
  -cg 3 4 \
  -ct md_strided.xtc \
  -cp topol.top \
  -o MMPBSA_results.dat \
  -eo MMPBSA_energy_per_frame.csv

# If instead you have the classic binary 'g_mmpbsa', uncomment and use:
# g_mmpbsa -f md_strided.xtc -s md.tpr -n index.ndx -i mmpbsa.in -decomp
# python MmPbSaStat.py -m energy_MM.xvg -p energy_PBSA.xvg -o FINAL_RESULTS_MMPBSA.dat

log "Outputs:"
echo " - MMPBSA_results.dat (summary ΔGbind)"
echo " - MMPBSA_energy_per_frame.csv (framewise energies)"

# record end time
end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
minutes=$(( elapsed / 60 ))
echo ""
echo "Total time: $minutes mins "
echo ""
