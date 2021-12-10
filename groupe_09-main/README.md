# Hackathon

## Installation
Before launching the pipeline, please enter the following commands:
```bash

# 1) To obtain the correct version of snakemake
conda install -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
#if "conda activate snakemake" does not work try "conda activate"

# 2) Check if you have the snakemake latest stable version : (should be 5.27.4 or 5.29.0)
snakemake --version

# 3) Download singularity
conda install singularity=3.6.3 -y
#if not working then try :
conda config --set allow_conda_downgrades true
#then try "conda install singularity=3.6.3" again

#4) Check that the launcher file is executable. If not :
chmod 777 launcher.sh

#5) Launch the pipeline
./launcher.sh
#This will install all the remaining necessary tools and launch the pipeline
#BEWARE : see Usage for details !
```
NB : (1) the command sudo apt install snakemake does not give the up-to-date version. Update from that command leads to bugs. (2) sudo apt install singularity corresponds to a game named singularity coded in Python.


## Usage
1. launcher.sh : file with commands necessary to run the pipeline correctly
2. pipeline.smk : file containing all steps (rules) needed to study if there is a differential expression between two types of melanoma samples (4 WT SF3B1, 4 mutated SF3B1).
3. deseq2_annot_visualisation.R : R script run into rule deseq of pipeline.

> All 3 files need to be downloaded in the working directory where the launcher is run.
> The pipeline has been designed to suit any types of machines. We tested the pipeline with option --cores 2. One can adapt the number of cores/CPUs in the launcher pipeline script. 

/!\ BEWARE : Because of a force interaction in rule configure, user has to wait for the job to appear (sra-configuration blue panel will pop-up). Press 'x' key (or ctrl+x if 'x' not sufficient) and the pipeline will continue running automatically.
/!\ BEWARE 2 : the pipeline requires some files to be in the same directory :
- launcher.sh (creation of some necessary empty files, and downloads of dependencies)
- deseq2_annot_visualisation.R (R script called for the analysis in the pipeline)

Estimated execution time of whole pipeline (with default option of --cores 2) : 3,8-4.2 hours

## Contributing
* Anaïs CHOPIN
* Anna KARPUKHINA
* Karim RAQBI
* Ombeline LAMER

## Licence
Permission is hereby granted, free of charge, to any person obtaining a copy of this project to deal without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute copies of the project.

Please quote the contributors and do not sell if it is a strict copy of the project.
