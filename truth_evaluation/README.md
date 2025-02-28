
# Overview

A pipeline to rerun the truth evaluation in the Amira paper on 32 *E. colI* samples.

# Installation

The pipeline assumes the Amira singularity container is available in the directory, which can be installed [here](https://github.com/Danderson123/amira). You also need to have conda installed. You will need to build the `kma` binary to run ResFinder. You can do this by running:
```{bash}
cd software/kma && nake && cd ../..
```
The remaining dependencies for the *E. coli* evaluation can be installed with:
```{bash}
conda env create -f envs/truth_env.yaml && conda activate truth_env
```
The *E. coli* panRG needs to be located in this directory and it can be downloaded from [here](https://drive.google.com/file/d/13c_bUXnBEs9iEPPobou7-xEgkz_t08YP/view?usp=sharing). You will also need to download the AMRFinderPlus database `v2024-01-31.1` from [here](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/2024-01-31.1/) and specify the path to it in by modifying the `amrfinder_db` parameter in the `Snakefile`.

# Running the evaluation
```{bash}
snakemake --cores 12 --use-conda --use-singularity --nolock --rerun-incomplete --keep-going 
```