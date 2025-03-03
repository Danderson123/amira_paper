# Amira paper analyses

This GitHub repository contains all of the code to r-erun each of the three analyses from the Amira paper. Each directory contains its own README that describes how to install the dependencies for and run that analysis. To get started you need to clone this repository using:
```{bash}
git clone --recursive https://github.com/Danderson123/amira_paper
```
The analyses take a while to run so it is highly recommended that you use a cluster computing environment if you want to run the pipelines. Tools to install Snakemake profiles for your favourite cluster computing environment can be found [here](https://github.com/Snakemake-Profiles). Once installed, you can then subsitute the `--cores <CORES>` command line option for each example snakemake command with e.g. `--profile slurm`.