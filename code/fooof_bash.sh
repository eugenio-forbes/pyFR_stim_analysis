#!/bin/bash

#Bash script for running fooof algorithm wrapper on different sets of data:
#-encoding or retrieval period
#-hippocampus or stim_site data
#-fooof (all trials) or bosc (trial mean) power spectra

#Change directory to folder containing analyisis code
cd /directory/with/pyFR_stim_analysis/main

#Activate conda environment with fooof installation
source /directory/with/anaconda/etc/profile.d/conda.sh
conda activate fooof_env

#Run fooof wrapper with all combinations of input arguments (task period, site, power spectra type)
python3 n15_fooof_processing.py encoding hippocampus fooof
python3 n15_fooof_processing.py encoding hippocampus bosc
python3 n15_fooof_processing.py retrieval hippocampus fooof
python3 n15_fooof_processing.py retrieval hippocampus bosc
python3 n15_fooof_processing.py encoding stim_site fooof
python3 n15_fooof_processing.py encoding stim_site bosc
python3 n15_fooof_processing.py retrieval stim_site fooof
python3 n15_fooof_processing.py retrieval stim_site bosc
