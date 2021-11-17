# -*- coding: utf-8 -*-
# pylint: disable=line-too-long, missing-module-docstring, unused-variable, undefined-variable, too-many-locals, invalid-name
import os
import sys
import glob
import time
from traceback import format_exc
from argparse import ArgumentParser
import subprocess
from scripts import sourmash_matches as sm
from scripts import genomes_matches as gm

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def main(argv=None):
    """
    Generates protein database containing non-redundant sequences from assembled contigs and matched genomes from UHGG catalogue
    """
    program_name = os.path.basename(sys.argv[0])
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    try:
        parser = ArgumentParser(
            description='Generate protein sequence datanase')
        parser.add_argument('--source',  type=str, required=True,
                            help='enter the source of genomes metadata as FTP_data or EBI_internal')
        parser.add_argument('--ref_dir',  type=dir_path, required=True,
                            help='full path of the directory with all the reference signatures')
        parser.add_argument('--query_dir',  type=dir_path, required=True,
                            help='full path of the directory with all the sample genomes')
        parser.add_argument('--k_size', type=str, required=False,
                            help='enter any k-size 21,31 or 51', default='31')
        parser.add_argument('--scale', type=str, required=False, help='enter any scaling factor ', default='1000')
        parser.add_argument('--metadata', type=str, required=True, help='path for file that has all genome paths, a .txt file')

        starttime = time.time()
        args = parser.parse_args()

        query_sig = os.path.join(args.query_dir, "signatures")
        if not os.path.isdir(query_sig):
            p=subprocess.Popen(' '.join(['mkdir ', query_sig]), shell = True)
        for file in os.listdir(args.query_dir):
            if (file.endswith(".fasta")) or (file.endswith(".fasta.gz")) or (file.endswith(".fa")):
                query_name = os.path.join(args.query_dir, file)
                basename = os.path.splitext(file)[0]
                out_query= os.path.join(query_sig, basename+".sig")
                if not os.path.isfile(out_query):
                    sm.sourmash_sig(query_name, out_query, args.k_size, args.scale)
                else:
                    print("{} is already present".format(out_query))
                    continue
        if glob.glob(args.ref_dir + "/**/*.sbt.json", recursive = True):
            print("signature index file is present")
            index_file = (glob.glob(args.ref_dir + "/**/*.sbt.json", recursive = True)[0])
            for file in os.listdir(query_sig):
                q_file = os.path.join(query_sig, file)
                if file.endswith(".sig"):
                    sm.sourmash_gather(query_sig, args.k_size, q_file, index_file)
        else:
            genomes_list_len = len(os.listdir(args.ref_dir))
            index_file = os.path.join(args.ref_dir, "genome_index_"+str(genomes_list_len)+".sbt.json")
            os.chdir(args.ref_dir)
            sm.signature_index(args.k_size, str(genomes_list_len), args.ref_dir)
            for file in os.listdir(query_sig):
                q_file = os.path.join(query_sig, file)
                if file.endswith(".sig"):
                    sm.sourmash_gather(query_sig, args.k_size, q_file, index_file)
        #working with .csv files
        genome_path = args.source
        for file in os.listdir(query_sig):
            if file.endswith(".csv"):
                unique=set()
                basename = os.path.splitext(file)[0]
                pan_genome_dir = os.path.join(query_sig, basename)
                contig_protein = os.path.join(args.query_dir, basename+".faa")
                all_protein_file=os.path.join(pan_genome_dir, "all_"+basename+".faa")
                concatenated_file = os.path.join(args.query_dir, "completed_"+basename+".faa")
                print("This is the folder with uhgg matches ", pan_genome_dir )
                if not os.path.isdir(pan_genome_dir):
                    p=subprocess.Popen(' '.join(['mkdir', pan_genome_dir]), shell = True)
                try:
                    matched_sp_names = gm.gen_match_list(file, query_sig)
                    for match in matched_sp_names:
                        gen_name = match.split(".")[0]
                        unique.add(gen_name)
                    print(" UHGG matches",unique)
                    if genome_path == "FTP_data":
                        print("collecting metadata from ftp")
                        gm.get_genomes_from_ftp(unique, args.metadata, pan_genome_dir)
                    elif genome_path == "EBI_internal":
                        gm.get_genomes_from_ebi(unique, args.metadata, pan_genome_dir)
                    else:
                        print("Genome metadata not found. Please enter the correct source")
                    if not os.path.isfile(all_protein_file):
                        gm.concat_proteins(pan_genome_dir)
                    if not os.path.isfile(contig_protein):
                        gm.contig_prod(args.query_dir, basename, all_protein_file)
                    if os.path.isfile(concatenated_file):
                        gm.uniq_proteins(args.query_dir, basename)
                except Exception:
                    continue
        print('runtime is {} seconds'.format(time.time() - starttime))
    except Exception as error:
        print(program_name + ": " + repr(error) + '\n' + format_exc() + '\n')
        raise error


if __name__ == "__main__":
    sys.exit(main())
