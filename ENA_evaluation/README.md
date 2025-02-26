
# Overview

Pipelines to run Amira, Amira `--no-filtering` and AMRFinderPlus with Flye on a large number of *E. coli* and *K. pneumoniae* reads in the ENA. 

# Installation

The pipeline assumes the Amira singularity container is available in the directory, which can be installed [here](https://github.com/Danderson123/amira). You also need to have conda installed.
The dependencies for the *E. coli* evaluation can be installed with:
```{bash}
conda env create -f Escherichia_coli/envs && conda activate E_coli_env
```
The dependencies for the *K. pneumoniae* evaluation can be installed with:

```{bash}
conda env create -f Klebsiella_pneumoniae/envs/K_pneumoniae_env.yaml && conda activate K_pneumoniae_env
```
The *E. coli* panRG can be downloaded from [here](https://drive.google.com/file/d/13c_bUXnBEs9iEPPobou7-xEgkz_t08YP/view?usp=sharing).
The *K. pneumoniae* panRG can be downloaded from [here](https://drive.google.com/file/d/1DYG3QW3nrQfSckIX9Vjbhbqz5bRd9W3j/view?usp=drive_link).
You will also need to install the AMRFinderPlus database.

# Running the *E. coli* evaluation

The *E. coli* evaluation can be run with:
```{bash}
cd Escherichia_coli && snakemake --cores 12 --use-conda --nolock --rerun-incomplete --keep-going
```
Some samples will fail to assemble so the result plots have to be generated separately with:
```{bash}
python3 scripts/make_result_plots.py
```

# Running the *K. pneumoniae* evaluation

The *K. pneumoniae* evaluation can be run with:
```{bash}
cd Klebsiella_pneumoniae && snakemake --cores 12 --use-conda --nolock --rerun-incomplete --keep-going
```
Some samples will fail to assemble so the result plots have to be generated separately with:
```{bash}
python3 scripts/make_result_plots.py
```