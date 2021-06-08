# -*- coding: utf-8 -*-
# pylint: disable=line-too-long, missing-module-docstring, unused-variable, undefined-variable, too-many-locals, invalid-name
import os
import sys
import time
from traceback import format_exc
from argparse import ArgumentParser
import shutil
import logging
from collections import defaultdict
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

def gen_match_list(f_input: str, w_dir: str) -> list:
    """
    Generates contigs names to search in database
    :param f_input:  input sample file which is a .csv file
    :param w_dir: path for the working directory

    :return : a list
    """
    file_name = os.path.join(w_dir, f_input)
    logging.info("genrating matched from %s", file_name)
    data = pd.read_csv(file_name, sep=',')
    unique=set()
    match_dir ={}
    for i in range(len(data)):
        if round(data["f_match"][i], 2) >=0.30:
            unique.add(data["name"][i])
            name_list = list(unique)
    if len(unique) >= 1:
        match_list=  [item.split("/")[-1] for item in name_list]
        logging.info("matches found %s", match_list)
        return match_list
    else:
        raise Exception("No matches found")

def get_genomes_from_ftp(species_names: list,  file_name:str, d_dir: str):
    """
    copy genomes from ftp location
    :param species_names:  list of mtches found
    :param file_name: path for the metadata file
    :param d_dir: path of the destination directory
    """
    uniq_sp_name=set()
    for item in species_names:
        sp_name = (item.split('.')[0]).strip()
        print(sp_name)
        ftp_dir_path = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgg_catalogue"
        # ftp_dir_path = "/Users/kaurs/Desktop/project_files/ERP104047"
        data = pd.read_csv(file_name,  dtype=str, sep=',')
        for i in range(len(data)):
            if sp_name == data['Species_rep'][i]:
                if not sp_name in uniq_sp_name:
                    uniq_sp_name.add(sp_name)
                    MGnify_acc = data['MGnify_accession'][i]
                    genome_dir_path = ftp_dir_path+"/"+MGnify_acc[:13]+"/"+MGnify_acc
                    print(genome_dir_path)
                    if os.path.isdir(genome_dir_path):
                        pan_genomes_dir = os.path.join(genome_dir_path,"pan-genome")
                        genome_dir = os.path.join(genome_dir_path,"genome")
                        if os.path.isdir(pan_genomes_dir):
                            pan_genome_file = os.path.join(pan_genomes_dir, "pan-genome.faa")
                            print(pan_genome_file)
                            output_file = os.path.join(d_dir, MGnify_acc+"_pg.faa")
                            dest = shutil.copy(pan_genome_file, output_file)
                            genome_file = os.path.join(genome_dir, MGnify_acc+".faa")
                            out_file = os.path.join(d_dir, MGnify_acc+"_sp.faa")
                            dest_file = shutil.copy(genome_file, out_file)
                        else:
                            if os.path.isdir(genome_dir):
                                protein_file = os.path.join(genome_dir, MGnify_acc+".faa")
                                output_file = os.path.join(d_dir, MGnify_acc+".faa")
                                dest = shutil.copy(protein_file, output_file)
def get_genomes_from_ebi(species_names: list,  file_name:str, d_dir: str):
    """
    copy genomes from EBI filesystem
    :param species_names:  list of mtches found
    :param file_name: path for the metadata file
    :param d_dir: path of the destination directory
    """
    distinct_sp_name = set()
    for sp_name in species_names:
        sp_name = sp_name.strip()
        with open(file_name) as fin:
            for line in fin:
                line = line.rstrip('\n')
                species_name = (line.split("/")[8:9][0]).strip()
                genome_name = ((line.split("/")[-1]).split(".")[0]).strip()
                if sp_name == species_name:
                    if not sp_name in distinct_sp_name:
                        distinct_sp_name.add(sp_name)
                        src_dir = "/".join(line.split("/")[:9])
                        roary_dir = os.path.join(src_dir,"roary")
                        if os.path.isdir(roary_dir):
                            print("copying pan-genome")
                            pan_protein_file = os.path.join(roary_dir, "pan_genome_reference.faa")
                            output_file = os.path.join(d_dir, species_name+"_pg.faa")
                            dest = shutil.copy(pan_protein_file, output_file)
                if ((sp_name == species_name) and (species_name == genome_name)):
                    species_protein_path = line.split(".")[0]+"_prokka"
                    protein_file = os.path.join(species_protein_path, species_name+".faa")
                    sp_outfile = os.path.join(d_dir, species_name+"_sp.faa")
                    print("copying matching species_rep genomes")
                    dest1 = shutil.copy(protein_file, sp_outfile)

def contig_prod(c_dir: str, assembly: str, sm_proteins: str):
    """
    predict ORFs from contig assemblies
    :param c_dir:  path of directory containing contigs
    :param assembly: basename of an assembly
    :param sm_proteins: path for concatenated sourmash protein file
    """
    os.chdir(c_dir)
    for file in os.listdir(c_dir):
        contig_file_name = (file.split(".")[0]).strip()
        if (contig_file_name == assembly) and (file.endswith(".fa") or file.endswith(".fasta")):
            file_name = os.path.join(c_dir, file)
            cmd_prod = "  ".join(["prodigal -p meta -a ", contig_file_name+".faa", " -c -d ", contig_file_name+".fna", " -m -f gff -o ", contig_file_name+".gff"," -i ", file_name])
            subprocess.call(cmd_prod, shell=True)
            with open(sm_proteins, 'r') as fp:
                data = fp.read()
            contig_protein = os.path.join(c_dir, contig_file_name+".faa")
            final_file = os.path.join(c_dir, "completed_"+contig_file_name+".faa")
            with open(contig_protein,'r') as fc:
                data2 = fc.read()
            data += "\n"
            data += data2
            with open(final_file,'a') as fout:
                fout.write(data)

def concat_proteins(protein_dir: str):
    """
    Concatenate all the matching protein files from UHGG catalogue
    :param protein_dir:  path of directory containing protein files
    """
    prot_folder_name = protein_dir.split("/")[-1]
    outfilename = os.path.join(protein_dir, "all_"+prot_folder_name+".faa")
    with open(outfilename, 'wb') as outfile:
        for file in os.listdir(protein_dir):
            if file.endswith(".faa"):
                filename = os.path.join(protein_dir, file)
                if not filename == outfilename:
                    with open(filename, 'rb') as readfile:
                        shutil.copyfileobj(readfile, outfile)
                else:
                    continue
def uniq_proteins(c_dir: str, assembly: str):
    """
    generate a file with unique protein sequences
    :param c_dir:  path of directory containing contigs
    :param assembly: basename of an assembly
    """
    unique_records = defaultdict(list)
    protein_file = os.path.join(c_dir, "completed_"+assembly+".faa")
    uniq_seqs = os.path.join(c_dir, "unique_"+assembly+".faa")
    for record in SeqIO.parse(protein_file, "fasta"):
        unique_records[str(record.seq)].append(record.id)
    final_seq = (SeqRecord(Seq(seqi), id="|".join(gi), name='', description='') for seqi, gi in unique_records.items())
    SeqIO.write(final_seq, uniq_seqs, 'fasta')


def main(argv=None):
    program_name = os.path.basename(sys.argv[0])
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    try:
        parser = ArgumentParser(description='extract proteins for all matching genomes from UHGG catalogue')

        parser.add_argument('--wdir',  required=True, help='enter the path of the directory for all .csv files')

        parser.add_argument('--contig_dir', required=True, help='path for the contigs_dir')

        parser.add_argument('--metadata', required=True, help='path for file that has all genome paths, a .txt file')

        args = parser.parse_args()

        starttime = time.time()
        for file in os.listdir(args.wdir):
            if file.endswith(".csv"):
                unique=set()
                basename = (file.split(".")[0]).strip()
                pan_genome_dir = os.path.join(args.wdir, basename)
                contig_protein = os.path.join(args.contig_dir, basename+".faa")
                all_protein_file=os.path.join(pan_genome_dir, "all_"+basename+".faa")
                concatenated_file = os.path.join(args.contig_dir, "completed_"+basename+".faa")
                print(pan_genome_dir )
                if not os.path.isdir(pan_genome_dir):
                    p=subprocess.Popen(' '.join(['mkdir', pan_genome_dir]), shell = True)
                try:
                    matched_sp_names = gen_match_list(file, args.wdir)
                    for match in matched_sp_names:
                        gen_name = match.split(".")[0]
                        unique.add(gen_name)
                    print(unique)
                    genome_path = input ("Enter the source of genomes metadata as FTP_data or EBI_internal: ")
                    if genome_path == "FTP_data":
                        get_genomes_from_ftp(unique, args.metadata, pan_genome_dir)
                    else:
                        get_genomes_from_ebi(unique, args.metadata, pan_genome_dir)
                    if not os.path.isfile(all_protein_file):
                        concat_proteins(pan_genome_dir)
                    if not os.path.isfile(contig_protein):
                        contig_prod(args.contig_dir, basename, all_protein_file)
                    if os.path.isfile(concatenated_file):
                        uniq_proteins(args.contig_dir, basename)
                except Exception:
                    continue
        print('runtime is {} seconds'.format(time.time() - starttime))
    except Exception as error:
        print(program_name + ": " + repr(error) + '\n' + format_exc() + '\n')
        raise error


if __name__ == "__main__":
    sys.exit(main())
