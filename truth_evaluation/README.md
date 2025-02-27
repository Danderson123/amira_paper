
# Overview

A pipeline to rerun the truth evaluation in the Amira paper on 32 *E. colI* samples.

# Installation

The pipeline assumes the Amira singularity container is available in the directory, which can be installed [here](https://github.com/Danderson123/amira). You also need to have conda installed.
The dependencies for the *E. coli* evaluation can be installed with:
```{bash}
conda env create -f envs/truth_env.yaml && conda activate truth_env
```
The *E. coli* panRG needs to be located in this directory and it can be downloaded from [here](https://drive.google.com/file/d/13c_bUXnBEs9iEPPobou7-xEgkz_t08YP/view?usp=sharing). You will also need to install the AMRFinderPlus database.

# Running the evaluation
```{bash}
snakemake --cores 12 --use-conda --nolock --rerun-incomplete --keep-going 
```