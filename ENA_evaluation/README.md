
# Overview


# Installation

Amira can be installed from [here](https://github.com/Danderson123/amira). You also need to have conda installed.
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

The *E. coli* reads can be downloaded with:
```{bash}
cd Escherichia_coli && python3 scripts/download_reads.py
```
The *E. coli* evaluation can then be run with:
```{bash}
snakemake --cores 12 --use-conda --nolock --rerun-incomplete --keep-going
```

# Running the *K. pneumoniae* evaluation

The *K. pneumoniae* reads can be downloaded with:
```{bash}
cd Klebsiella_pneumoniae && python3 scripts/download_reads.py
```
The *K. pneumoniae* evaluation can then be run with:
```{bash}
nakemake --cores 12 --use-conda --nolock --rerun-incomplete --keep-going
```