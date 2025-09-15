#!/usr/bin/env bash
# Command-line tool for PRODIGY

$1=directory_with_molecules # /example/complex

# directory_with_molecules="/path/to/PDB_dir"
prodigy ${directory_with_molecules}

# multi_model_file="/binder_target_complex.pdb" # (an ensemble)
# np=3 # number of processors to use 
# prodigy ${multi_model_file} -np ${np}