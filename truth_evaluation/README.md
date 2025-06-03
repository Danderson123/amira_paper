
# Overview

A pipeline to rerun the truth evaluation in the Amira paper on 32 *E. coli* samples.

# Data

The pipeline assumes you have downloaded the truth evaluation data and it is available in this directory in a subdirectory called `amira_paper_data`. The reference assemblies can be downloaded from [here](https://figshare.com/ndownloader/articles/28958849/versions/1). You will also need to download the Illumina and Nanopore reads and place them in subdirectories called `amira_paper_data/illumina_reads` and `amira_paper_data/nanopore_reads` respectively. The *E. coli* panRG needs to be located in this directory as well and it can be downloaded from [here](https://figshare.com/ndownloader/files/54318899). You will also need to download the AMRFinderPlus database `v2024-01-31.1` from [here](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/2024-01-31.1/) and specify the path to it in by modifying the `amrfinder_db` parameter in the `Snakefile`.

# Installation

The pipeline assumes the Amira singularity container is available in the directory and is called `amira.v0.9.3.img`. This can be installed [here](https://github.com/Danderson123/amira). You also need to have conda installed. You will need to build the `kma` binary to run ResFinder. You can do this by running:
```{bash}
cd software/kma && make && cd ../..
```
The remaining dependencies for the *E. coli* evaluation can be installed with:
```{bash}
conda env create -f envs/truth_env.yaml && conda activate truth_env
```

# Sub-sampling the reads

The pipeline expects the sub-sampled Nanopore reads to be available in `amira_paper_data/subsampled_nanopore_reads`. This can be done automatically after downloaded the un-sampled Nanopore reads (and placing them in `amira_paper_data/nanopore_reads`) by running:
`python3 software/downsample_reads.py`

# Running the evaluation
```{bash}
snakemake --cores 12 --use-conda --use-singularity --nolock --rerun-incomplete --keep-going 
```
