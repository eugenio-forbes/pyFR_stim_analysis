# pyFR_stim_analysis
Code for performing analysis of UPenn's Restoring Active Memory free recall deep brain stimulation experiments.

## Requirements

### Software
- MATLAB (compatible with MATLAB 2020a and later versions)
- Python
- Anaconda

### MATLAB Toolboxes
- Machine Learning and Statistics Toolbox: Required for machine learning algorithms and statistical analysis.
- Parallel Computing Toolbox: Necessary for parallel processing to improve computation efficiency.
- Signal Processing Toolbox: Used for advanced signal processing tasks.

Please ensure that these toolboxes are installed and properly configured in your MATLAB environment to use all the functionalities of this project.

## Installation

### Step 1: Clone the Repository
Begin by cloning this repository to your local machine. Cloning this repository will also include code from other repositories that are necessary for the proper functioning of this code.
A copy of these repositories has been included to avoid issues with updates. 
For more information on these repositories, visit:
- [fooof-tools/fooof](https://github.com/fooof-tools/fooof)
- [jkosciessa/eBOSC](https://github.com/jkosciessa/eBOSC)
- [neurofractal/fBOSC](https://github.com/neurofractal/fBOSC)

### Step 2: Set Up the Conda Environment
To ensure a consistent environment with all necessary dependencies, please create a Conda environment using the provided `fooofpy.yml` file located in `/code/base_code`:
```
conda env create -f /code/base_code/fooofpy.yml
conda activate fooofpy
```
After completing these steps, your environment should be set up with all the necessary dependencies to run fooof within the analysis code. If you encounter any issues, please feel free to open an issue in this repository.
