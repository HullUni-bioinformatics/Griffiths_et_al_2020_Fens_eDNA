# Griffiths *et al.* (2020)
Data processing workflow and supplementary info for:
Griffiths *et al.* (2020) - **Environmental DNA metabarcoding provides enhanced detection of the European eel *Anguilla anguilla* and fish community structure in pumped river catchments**

Permanently archived at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3951418.svg)](https://doi.org/10.5281/zenodo.3951418)

## Contents

- NCBI Sequence Read Archive (SRA) accession numbers for raw Illumina data [(here)](https://github.com/NPGriffiths/Griffiths_et_al_2020_Fens/blob/master/Data/Sample_accessions.tsv)

- Curated reference databases used in analyses (GenBank/fasta format) [(here)](https://github.com/NPGriffiths/Griffiths_et_al_2020_Fens/tree/master/Reference_DBs)

- Notebook to run metaBEAT pipeline [(here)](https://github.com/NPGriffiths/Griffiths_et_al_2020_Fens/blob/master/Jupyter_notebooks/Fens_2017_metaBEAT.ipynb)

- Raw output data [(here)](https://github.com/NPGriffiths/Griffiths_et_al_2020_Fens/blob/master/Data/Raw_Data)

- R scripts with associated .csv files used to analyse metaBEAT output and produce figures [(here)](https://github.com/NPGriffiths/Griffiths_et_al_2020_Fens/tree/master/R_scripts)

- Supplementary info [(here)](https://github.com/NPGriffiths/Griffiths_et_al_2020_Fens/tree/master/Data/Supp_Info)


## Instructions to set up dependencies for data processing and analyses
To facilitate full reproducibility of our analyses, we provide Jupyter notebooks illustrating our workflow and all necessary associated data in this repository.

llumina outputs were processed (from raw reads to taxonomic assignment) using the [metaBEAT](https://github.com/HullUni-bioinformatics/metaBEAT) pipeline. The pipeline relies on a range of open bioinformatics tools, which we have wrapped up in a self-contained docker image that includes all necessary dependencies [here](https://hub.docker.com/r/chrishah/metabeat/).

## Setting up the environment
In order to retrieve scripts and associated data (reference sequences, sample metadata etc.), start by cloning this repository to your current directory:

```
git clone --recursive https://github.com/HullUni-bioinformatics/Griffiths_et_al_2020_Fens.git
```
In order to make use of our self contained analysis environment, you will have to install [Docker](https://www.docker.com/) on your computer. 
Docker is compatible with all major operating systems, but see the [Docker documentation](https://docs.docker.com/) for details. On Ubuntu, installing Docker should be as easy as:
```
sudo apt-get install docker.io
```

Once Docker is installed, you can enter the environment by typing for example:
```
sudo docker run -i -t --net=host --name metaBEAT -v $(pwd):/home/working chrishah/metabeat /bin/bash
```
This will download the metaBEAT image (if not yet present on your computer) and enter the 'container' i.e. the self contained environment (**NB:** ```sudo``` may be necessary in some cases). With the above command, the container's directory ```/home/working``` will be mounted to your current working directory (as instructed by ```$(pwd)```). In other words, anything you do in the container's ```/home/working``` directory will be synced with your current working directory on your local machine.


## Data processing workflow as Jupyter notebooks
Raw illumina data has been deposited on the NCBI SRA:
- Study: SRP272137
- BioProject: PRJNA646357
- BioSample accessions: SAMN15541168 - SAMN15541275
- SRA accessions: SRR12232432 - SRR12232539 

The sample specific accessions can be found [here](https://github.com/NPGriffiths/Griffiths_et_al_2020_Fens/blob/master/Data/Sample_accessions.tsv). Before following the workflow for data processing, you'll need to download the raw reads from the SRA. To download the raw read data, you can follow the steps in this [Jupyter notebook](https://github.com/NPGriffiths/Griffiths_et_al_2020_Fens/blob/master/Jupyter_notebooks/How_to_download_from_SRA.ipynb).

With the data in place, you should be able to fully reproduce our analyses by following the steps outlined in the [Jupyter notebook](https://github.com/NPGriffiths/Griffiths_et_al_2020_Fens/blob/master/Jupyter_notebooks/Fens_2017_metaBEAT.ipynb).

The workflow illustrated in the notebooks assumes that the raw Illumina data is present in a directory ```raw_reads``` at the base of the repository structure and that the files are named according to the following convention: 'sampleID_marker', followed by '_R1' or '_R2' to identify the forward/reverse read file respectively. SampleID must correspond to the first column in the file ```Sample_accessions.tsv``` [here](https://github.com/NPGriffiths/Griffiths_et_al_2020_Fens/blob/master/Data/Sample_accessions.tsv).
