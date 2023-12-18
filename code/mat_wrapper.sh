#!/bin/bash

# Bash script for running power analysis in the SLURM

# Change to directory with analysis code
ROOT_DIR="/base/directory/with/ANALYSIS_FOLDER"
ANALYSIS_FOLDER="pyFR_stim_analysis"

cd $ROOT_DIR
cd $ANALYSIS_FOLDER
cd code

# Load MATLAB module
module load matlab/2020a

# Run MATLAB code without display
matlab -nodesktop -nodisplay -singleCompThread -r "n18_power_analysis('$ROOT_DIR', '$ANALYSIS_FOLDER'); exit;"
