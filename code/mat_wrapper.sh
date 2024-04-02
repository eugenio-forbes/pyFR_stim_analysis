#!/bin/bash

# Node start time, parent directory, and custom preferences for rereferencing,
# denoising, clustering

ROOT_DIR="/directory/to/pyFR_stim_analysis"

cd $ROOT_DIR
cd code
module load matlab

# Run MATLAB code without display
matlab -nodesktop -nodisplay -singleCompThread -r "n18_power_analysis('$ROOT_DIR'); exit;"

