
# Overview

Simulations of *E. coli*, using the *Escherichia coli str. K-12 substr. MG1655* reference genome (Accession CU00096.3) and the *Escherichia coli ATCC 11775* plasmid (Accession CP033091.2). The AMR_contexts directory contain fastas of modules that are inserted into the genomes to create synthetic references, from which reads are simulated for evaluation of Amira.

# Installation

The pipeline assumes the Amira singularity container is available in the directory and is called `amira.v0.9.3.img`. This can be installed [here](https://github.com/Danderson123/amira). You also need to have conda and singularity installed. All other dependencies can be installed via conda with:

```{bash}
conda env create -f envs/sims_env.yaml && conda activate simulation_env
```

The pipeline assumes the panRG is avilable in this directory and it can be downloaded from [here](https://figshare.com/ndownloader/files/54318899).

You will also need to install the AMRFinderPlus database.

# Running the simulations

There are four Snakefiles in this directory, each corresponding to simulation for different mean read lengths (5kb, 10kb, 20kb and 40kb). For example, to rerun the 5kb simulations run:
```{bash}
snakemake --snakefile Snakefile_5kb --cores 12 --use-conda --nolock --rerun-incomplete --keep-going 
```

# Generating results

The final plots can be generated with this command after running the simulations for all read lengths:
```{bash}
python3 scripts/make_sim_results.py
```