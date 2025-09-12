
![toc_pacs_large](https://github.com/user-attachments/assets/e6e0acbe-2009-410c-afbd-6bcbe83171a0)

Peptide drugs make up a growing proportion of the pharmaceutical market. However, computational prediction of peptide-protein binding affinity remains difficult, which limits the use of computational methods to assist the optimization of bioactive peptides. We have developed a novel method for modelling peptide unbinding, designated contact parallel cascade selection molecular dynamics (cPaCS-MD), which, when combined with Markov state models, accurately predict peptide binding affinities. We applied cPaCS-MD to a diverse set of twelve protein-peptide complexes and found that it predicts experimental peptide binding free energy, with a strong correlation (R<sup>2</sup> = 0.84) and high accuracy (mean absolute error of 2.7 kJ/mol and root mean squared error and 3.4 kJ/mol). We compared cPaCS-MD to the widely used steered molecular dynamics-umbrella sampling (SMD-US) method and found that PaCS-MD is more computationally efficient and of much greater accuracy. This work is the first benchmarking study of MD methods for predicting peptide binding affinities that uses a diverse set of peptide ligands.
 
# Executing PaCS-MD Simulations

PaCS-MD simulations are executed using the `cpacs_script.py` provided. Below are the steps required to set up your environment and run the simulations.

## Prerequisites

1. **Gromacs Installation:**
   PaCS-MD requires Gromacs for running and analyzing simulations. Please ensure that Gromacs is installed and the appropriate directory is linked in the script.

2. **PaCS Environment Installation:**
   To set up the PaCS environment, use the `pacs.yml` file provided in the repository. It is recommended to use Conda for environment management.

## Installation

1. **Install Conda (if not installed):**
   If Conda is not already installed, download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/individual).

2. **Create the Conda environment:**
   In the terminal, navigate to the directory containing the `pacs.yml` file and run the following command to create the environment:

   ```bash
   conda env create -f pacs.yml

3. **Activate the Conda environemnt:**
   Once the environment is created, activate it with:
    ```bash
   conda activate pacs

## Running Simulations ##
After setting up the environment, navigate to the directory containing cpacs_script.py, your and your equilibrated gromacs input file, as well as the appropriate .mdp file, and run the simulation with the following command:
 ```bash
python cpacs_script.py
