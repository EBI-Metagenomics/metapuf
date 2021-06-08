# MetaPUF
# Task : generating protein databases
##  Setup
# To use metapuf the following dependencies need to be satisfied
- python 3.7 or above
- Sourmash 4.1
- Prodigal
- mmseqs2

# you can create a virtual environment and install the dependencies using metapuf.yml

## The script has upper level ftp  location hard coded, which can be edited

## The folder structure is expected to be similar to the one in ftop location

## metadata file contains all the information about the genomes and can be found at http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/

## To install metapuf, run:
- git clone
- cd metapuf
-  if you are using conda make sure to activate the environment

## Usage
The signature file for UHGG catalogue can be found /ebi/ftp/pub/databases/metagenomics/metapuf/UHGG_catalogue
#### run scripts/sourmash_matches.py to query the metagenome against UHGG calatogue, generating a .csv file

- sourmash_matches.py --ref_dir {path for UHGG_catalogue genomes} --dest_dir {path for saving ref genome signatures} --query_dir {path for query genomes} --k_size {default 31} --scale {default 1000}

#### run genomes_matches.py to copy matching genomes and pan-genomes from ftp location

- genomes_matches.py --wdir {path for all .csv files} --contig_dir {path for contigs directory} --metadata {path for metadata}
