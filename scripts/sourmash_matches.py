# -*- coding: utf-8 -*-
# pylint: disable=line-too-long, missing-module-docstring, unused-variable, undefined-variable, too-many-locals, invalid-name
import os
import sys
import subprocess
import glob
from traceback import format_exc
from argparse import ArgumentParser



def sourmash_sig(file_name: str, out_file: str, k_size: str, scale: str):
    """
    Compute sourmash signatures for the .fasta files
    :param file_name: path of the file
    :param out_dir: path of the directory where signatures are stored
    :param k_size: size of k-mer, the default value being 31
    :param scale: scaling factor, the default value is 1000

    """
    print("computing sourmash signatures ....")
    cmd_sig = "  ".join(["sourmash compute -k ", k_size, " --scaled ", scale, " --track-abundance -o ", out_file, file_name])
    subprocess.call(cmd_sig, shell=True)

def signature_index(k_size: str, file_list: list, file_list_len: str):
    """
    Compute index for reference genomes
    :param file_list: list of genome names
    :param k_size: size of k-mer, the default value being 31
    """
    print("computing sourmash index....")
    cmd_index = "  ".join(["sourmash index -k ", k_size, " genome_index_"+file_list_len , file_list])
    subprocess.call(cmd_index, shell=True)

def sourmash_gather(w_dir: str, k_size: str, query: str, index_file: str):
    """
    find the best match of query against reference geneomes
    :param w_dir: path of output .csv and .sm files
    :param k_size: size of k-mer, the default value being 31
    :param query: file path of the query signature file
    :param index_file: file path for index of all reference genomes
    """
    print("finding best matches using sourmash......")
    csv_out = os.path.join(w_dir, query+".csv")
    sm_out = os.path.join(w_dir, query+".sm")
    cmd_gather = "  ".join(["sourmash gather -k ", k_size, query, index_file, " -o ", csv_out , " > ", sm_out])
    subprocess.call(cmd_gather, shell=True)

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def main(argv=None):
    """
    Computes sourmash signatures and performs search for best match in reference databasesa
    """
    program_name = os.path.basename(sys.argv[0])
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    try:
        parser = ArgumentParser(
            description='Find matches from reference database based on sourmash signatures')

        parser.add_argument('--ref_dir',  type=dir_path, required=True,
                            help='path of the directory with all the reference genomes')
        parser.add_argument('--query_dir',  type=dir_path, required=True,
                            help='path of the directory with all the sample genomes')
        parser.add_argument('--k_size', type=str, required=False,
                            help='enter any k-size 21,31 or 51', default='31')
        parser.add_argument('--scale', type=str, required=False, help='enter any scaling factor ', default='1000')

        args = parser.parse_args()
        #directory that conatins all reference signatures
        signature_list=[]
        sig_dir = os.path.join(args.ref_dir, "signatures")
        if not os.path.isdir(sig_dir):
            p=subprocess.Popen(' '.join(['mkdir ', sig_dir]), shell = True)
        #directory that contains all query signatures
        query_sig = os.path.join(args.query_dir, "signatures")
        if not os.path.isdir(query_sig):
            p=subprocess.Popen(' '.join(['mkdir ', query_sig]), shell = True)
        #reference signatures
        for file in os.listdir(args.ref_dir):
            if file.endswith(".fa"):
                ref_name = os.path.join(args.ref_dir, file)
                basename = file.split(".")[0]
                out_ref = os.path.join(sig_dir, basename+".sig")
                if not os.path.isfile(out_ref):
                    sourmash_sig(ref_name, out_ref,  args.k_size, args.scale)
                else:
                    print("{} is already present".format(out_ref))
                    continue
        # query signatures
        for file in os.listdir(args.query_dir):
            if file.endswith(".fasta"):
                query_name = os.path.join(args.query_dir, file)
                basename = file.split(".")[0]
                out_query= os.path.join(query_sig, basename+".sig")
                if not os.path.isfile(out_query):
                    sourmash_sig(query_name, out_query, args.k_size, args.scale)
                else:
                    print("{} is already present".format(out_query))
                    continue
        if glob.glob(args.ref_dir + "/**/*.sbt.json", recursive = True):
            print("signature index file is present")
            index_file = (glob.glob(args.ref_dir + "/**/*.sbt.json", recursive = True)[0])
            for file in os.listdir(query_sig):
                q_file = os.path.join(query_sig, file)
                if q_file.endswith(".sig"):
                    sourmash_gather(query_sig, args.k_size, q_file, index_file)
        else:
            for f_in in os.listdir(sig_dir):
                f_name = os.path.join(sig_dir,f_in)
                if f_name.endswith(".sig"):
                    signature_list.append(f_name)
            if len(signature_list) < 2:
                print("Less than 2 input signatures provided...")
                sys.exit(1)
            else:
                genome_list = " ".join(signature_list[:])
            genome_list_len = str(len(signature_list))
            os.chdir(sig_dir)
            index_file = os.path.join(sig_dir, "genome_index_"+genome_list_len+".sbt.json")
            signature_index(args.k_size, genome_list, genome_list_len)
            for file in os.listdir(query_sig):
                q_file = os.path.join(query_sig, file)
                if q_file.endswith(".sig"):
                    sourmash_gather(query_sig, args.k_size, q_file, index_file)

    except Exception as e:
        print(program_name + ": " + repr(e) + '\n' + format_exc() + '\n')
        raise(e)


if __name__ == "__main__":
    sys.exit(main())
