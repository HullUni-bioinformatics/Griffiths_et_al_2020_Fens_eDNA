# Griffiths *et al.* (2020)
Data processing workflow and supplementary info for:
Griffiths *et al.* (2020) - **Environmental DNA metabarcoding provides enhanced detection of the European eel Anguilla anguilla and fish community structure in pumped river catchments**

Permanently archived at:

## Contents

## Instructions to set up dependencies for data processing and analyses
To facilitate full reproducibility of our analyses, we provide Jupyter notebooks illustrating our workflow and all necessary associated data in this repository.

Raw llumina outputs were processed (from raw reads to taxonomic assignment) using the [metaBEAT](https://github.com/HullUni-bioinformatics/metaBEAT) pipeline. The pipeline relies on a range of open bioinformatics tools, which we have wrapped up in a self-contained docker image that includes all necessary dependencies [here](https://hub.docker.com/r/chrishah/metabeat/).

## Setting up the environment
In order to retrieve scripts and associated data (reference sequences, sample metadata etc.), start by cloning this repository to your current directory:

```
git clone --recursive https://github.com/HullUni-bioinformatics/Griffiths_et_al_2020_Fens.git
```
In order to make use of our self contained analysis environment, you will have to install [Docker](https://www.docker.com/) on your computer. 
Docker is compatible with all major operating systems, but see the [Docker documenation](https://docs.docker.com/) for details. On Ubuntu, installing Docker should be as easy as:
```
sudo apt-get install docker.io
```

Once Docker is installed, you can enter the environment by typing, e.g.:
```
sudo docker run -i -t --net=host --name metaBEAT -v $(pwd):/home/working chrishah/metabeat /bin/bash
```
This will download the metaBEAT image (if not yet present on your computer) and enter the 'container' i.e. the self contained environment (**NB:** ```sudo``` may be necessary in some cases). With the above command, the container's directory ```/home/working``` will be mounted to your current working directory (as instructed by ```$(pwd)```). In other words, anything you do in the container's ```/home/working``` directory will be synced with your current working directory on your local machine.


## Data processing workflow as Jupyter notebooks
Raw illumina data has been deposited on the NCBI SRA:
- BioProject: 
- BioSample accessions: 
- SRA accessions: 

The sample specific accessions can be found [here](https://github.com/HullUni-bioinformatics). Before following the workflow for data processing, you'll need to download the raw reads from the SRA. To download the raw read data, you can follow the steps in this [Jupyter notebook](https://github.com/HullUni-bioinformatics).

With the data in place, you should be able to fully reproduce our analyses by following the steps outlined in the [Jupyter notebook](https://github.com/HullUni-bioinformatics).

The workflow illustrated in the notebooks assumes that the raw Illumina data is present in a directory ```raw_reads``` at the base of the repository structure and that the files are named according to the following convention: 'sampleID', followed by '_R1' or '_R2' to identify the forward/reverse read file respectively. SampleID must correspond to the first column in the file ```Sample_accessions.tsv``` [here](https://github.com/HullUni-bioinformatics).
