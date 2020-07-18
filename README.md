# Griffiths_et_al_2020_Fens

Data processing workflow and supp info for:

"INSERT TITLE" 

Permanently archived at:

# Contents 




# Bioinformatics processing and analyses
To facilitate full reproducibility of our analyses, we provide Jupyter notebooks illustrating our workflow and all necessary associated data in this repository.

Illumina output was processed (from raw reads to taxonomic assignment) using the metaBEAT pipeline. 
The pipeline relies on a range of open bioinformatics tools, which we have wrapped up in a self-contained docker image that includes all necessary dependencies.

# Data processing workflow

Raw illumina data has been deposited on the NCBI SRA:
BioProject:
BioSample accessions:
SRA accessions:
The sample specific accessions can be found here. 
Before following the workflow for data processing, you'll need to download the raw reads from the SRA.

With the data in place, you should be able to fully reproduce our analyses by following the steps outlined in the Jupyter notebook. 

The workflow illustrated in the notebooks assumes that the raw Illumina data is present in a directory raw_reads at the base of the repository structure and that the files are named according to the following convention: 'sampleID', followed by 'R1' or 'R2' to identify the forward/reverse read file respectively. SampleID must correspond to the first column in the file Sample_accessions.tsv.

