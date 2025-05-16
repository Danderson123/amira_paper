
# Overview

Pipelines to run Amira and AMRFinderPlus with Flye on a large number of *E. coli*, *K. pneumoniae* and *E. faceium" samples with long reads from the ENA.

# Installation

The pipeline assumes the Amira singularity container is available in the directory and is called `amira.v0.9.2.img`. This can be installed [here](https://github.com/Danderson123/amira). You also need to have conda installed.
The dependencies for the *E. coli* evaluation can be installed with:
```{bash}
conda env create -f Escherichia_coli/envs/E_coli_env.yaml && conda activate E_coli_env
```
The dependencies for the *K. pneumoniae* evaluation can be installed with:

```{bash}
conda env create -f Klebsiella_pneumoniae/envs/K_pneumoniae_env.yaml && conda activate K_pneumoniae_env
```
The dependencies for the *E. faecium* evaluation can be installed with:

```{bash}
conda env create -f Enterococcus_faecium/envs/E_faecium_env.yaml && conda activate E_faecium_env
```
The pipeline assumes the *E. coli* panRG is located in `Escherichia_coli`. It can be downloaded from [here](https://figshare.com/ndownloader/files/54318899).
The pipeline assumes the *K. pneumoniae* panRG is located in `Klebsiella_pneumoniae`. It can be downloaded from [here](https://figshare.com/ndownloader/files/53398349).
The pipeline assumes the *E. faecium* panRG is located in `Enterococcus_faecium`. It can be downloaded from [here](https://figshare.com/ndownloader/files/53395052).
You will also need to install the AMRFinderPlus database.

# Running the *E. coli* evaluation

The *E. coli* evaluation can then be run with:
```{bash}
cd Escherichia_coli && python3 scripts/download_reads.py && snakemake --cores 12 --use-conda --nolock --rerun-incomplete --keep-going
```
Some samples will fail to assemble so the result plots have to be generated separately. You will need to make a separate conda environment for the plotting dependencies. This can be done with:
```{bash}
conda env create -f envs/plot_results.yaml && conda activate E_coli_plot_env
```
The plots can then be generated with:
```{bash}
python3 scripts/make_result_plots.py
```
# Running the *K. pneumoniae* evaluation

The *K. pneumoniae* evaluation can then be run with:
```{bash}
cd Klebsiella_pneumoniae && python3 scripts/download_reads.py && snakemake --cores 12 --use-conda --nolock --rerun-incomplete --keep-going
```
Some samples will fail to assemble so the result plots have to be generated separately. You will need to make a separate conda environment for the plotting dependencies. This can be done with:
```{bash}
conda env create -f envs/plot_results.yaml && conda activate K_pneumoniae_plot_env
```
The plots can then be generated with:
```{bash}
python3 scripts/make_result_plots.py
```

# Running the *E. faecium* evaluation

The *E. faecium* evaluation can then be run with:
```{bash}
cd Enterococcus_faecium && python3 scripts/download_reads.py && snakemake --cores 12 --use-conda --nolock --rerun-incomplete --keep-going
```
Some samples will fail to assemble so the result plots have to be generated separately. You will need to make a separate conda environment for the plotting dependencies. This can be done with:
```{bash}
conda env create -f envs/plot_results.yaml && conda activate E_faecium_plot_env
```
The plots can then be generated with:
```{bash}
python3 scripts/make_result_plots.py
```