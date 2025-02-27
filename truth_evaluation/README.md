
# Overview


# Installation

The pipeline assumes the Amira singularity container is available in the directory, which can be installed [here](https://github.com/Danderson123/amira). You also need to have conda installed.
The dependencies for the *E. coli* evaluation can be installed with:
```{bash}
conda env create -f Escherichia_coli/envs && conda activate E_coli_env
```
The *E. coli* panRG can be downloaded from [here](https://drive.google.com/file/d/13c_bUXnBEs9iEPPobou7-xEgkz_t08YP/view?usp=sharing).
You will also need to install the AMRFinderPlus database.

# Running the evaluation
