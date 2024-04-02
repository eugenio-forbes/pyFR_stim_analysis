#!/bin/bash

cd /directory/to/pyFR_stim_analysis/code

source /directory/to/python/3.6.4-anaconda/etc/profile.d/conda.sh

conda activate fooofpy

python3 n15_fooof_processing.py encoding hippocampus fooof
python3 n15_fooof_processing.py encoding hippocampus bosc
python3 n15_fooof_processing.py retrieval hippocampus fooof
python3 n15_fooof_processing.py retrieval hippocampus bosc
python3 n15_fooof_processing.py encoding stim_site fooof
python3 n15_fooof_processing.py encoding stim_site bosc
python3 n15_fooof_processing.py retrieval stim_site fooof
python3 n15_fooof_processing.py retrieval stim_site bosc
