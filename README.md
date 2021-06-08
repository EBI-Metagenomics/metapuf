# MetaPUF
## Task : generating protein databases
##  Setup
### To use metapuf the following dependencies need to be satisfied
- python 3.7 or above
- Sourmash 4.1
- Prodigal
- mmseqs2
- Biopython

If you don't have conda on your computer you can install the Miniconda Python3 distribution from this link (https://conda.io/en/latest/miniconda.html). Make sure you answer yes to the question whether conda shall be put into your PATH. You will need the following channels
a). bioconda,
b). conda-forge
you can create a virtual environment and install the dependencies using metapuf.yml
- conda env create -n metapuf -f metapuf.yml

The script has upper level ftp location hard coded, which can be edited. The folder structure is expected to be similar to the one in ftp location. for e.g. ../human-gut/v1.0/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00002/genome/

Metadata file contains all the information about the genomes and can be found at http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/

#### To install metapuf, run:
- git clone https://github.com/EBI-Metagenomics/metapuf.git
- cd metapuf
-  if you are using conda make sure to activate the environment

### Usage
The signature file for UHGG catalogue can be found /ebi/ftp/pub/databases/metagenomics/metapuf/UHGG_catalogue. The workflow assumes that the reference signatures are already present.
#### run run_task.py to generate a protein database with unique sequences
- run_task.py  --ref_dir {path of ref genome signatures} --query_dir {path for query genomes} --k_size {default 31} --scale {default 1000} --contig_dir {path for contigs directory} --metadata {path for metadata}

#### run mmseqs-wf.sh to generate clusters. Create an output foler to save reults

- mmseqs_wf.sh -t {no of threads} -f {input protein file} -i {amino acid identity} -c {sequence identity} -o {output_folder}"

- The representative sequences can be found at output_folder/mmseqs_cluster_rep.fa
