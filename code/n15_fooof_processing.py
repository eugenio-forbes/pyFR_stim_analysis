import sys
import os
import json
import pandas as pd
import numpy as np
import h5py
from scipy.io import loadmat,savemat
from fooof import FOOOF, FOOOFGroup
from random import randint
import time

# To fill in manually or receive input arguments
if len(sys.argv) <= 1:
    event_type = 'encoding'
    region_of_interset = 'hippocampus'
    fooof_mode = 'fooof'
else:
    event_type = sys.argv[1]
    region_of_interest = sys.argv[2]
    fooof_mode = sys.argv[3]

# Modify below to root directory this script should operate on
root_directory = '/project/TIBIR/Lega_lab/s427026/pyFR_stim_' + folder_name

# Determine other directories based on input arguments. Create new directories if needed for saving locks and 
list_directory = root_directory + '/lists'
list_text_file = list_directory + '/' + event_type + '_fooof_directories_' + region_of_interest + '.txt'
parent_directory = root_directory + '/locks'
lock_directory = parent_directory + '/' + event_type + '_' + fooof_mode + '_locks_' + region_of_interest
if not os.path.isdir(parent_directory):
    os.mkdir(parent_directory)
if not os.path.isdir(lock_directory):
    os.mkdir(lock_directory)

# Read .txt file with list of directories with data to be processed
list_file = open(list_text_file,'r')
electrode_directories = list(list_file.read().split('\n'))

# Initialize counter to create lock files with this number (index of processed directory)
count = 0

# Loop through directories with data to be processed
for electrode_directory in electrode_directories:
    
    # Create filename for lock, done, and error, files
    count = count +1
    lock_file = lock_directory + '/py_lock_' + str(count) + '_LOCK'
    done_file = lock_directory + '/py_lock_' + str(count) + '_DONE'
    error_file = lock_directory + '/py_lock_' + str(count) + '_ERROR'
    
    # Sleep a random number of seconds (between 1 and 5) to separate 
    # multiple instances of code execution run simulatenously
    time.sleep(randint(1,5))
    
    # Only process data if lock file for it has not been created
    if not os.path.exists(lock_file):
        f = open(lock_file,'w')
        f.close()
        
        try:
            data_path = electrode_directory + '/' +event_type + '_' + fooof_mode 
            mat_file = data_path + '/data.mat'

            print(count)
            print(mat_file)

            start_time = time.time()

            with h5py.File(mat_file, 'r') as file:
                frequencies = np.array(file['frequencies']).squeeze().astype('float')

            if fooof_mode == 'fooof':
                with h5py.File(mat_file, 'r') as file:
                    power_spectra = np.array(file['power_spectra']).squeeze().transpose().astype('float')
                fg = FOOOFGroup(aperiodic_mode='knee',peak_width_limits = [0.3, 24],peak_threshold = 1,verbose = False)
                fg.fit(frequencies,power_spectra,[1,45],n_jobs=-1)
            else:
                with h5py.File(mat_file, 'r') as file:
                    mean_power_spectrum = np.array(file['mean_power_spectrum']).squeeze().astype('float')
                fg = FOOOF(aperiodic_mode='knee',peak_width_limits = [0.3, 24],peak_threshold = 1,verbose = False)
                fg.fit(frequencies,mean_power_spectrum,[1,45])
                     
            os.chdir(data_path)
            
            fg.save('fooof_results',save_results = True)
            
            end_time = time.time()
            print(end_time-start_time)

            f = open(done_file,'w')
            f.close()
        except Exception as e:
            with open(error_file,'w') as file:
                file.write(f"An error occurred: {e}\n")
