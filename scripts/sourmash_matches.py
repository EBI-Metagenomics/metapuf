# -*- coding: utf-8 -*-
# pylint: disable=line-too-long, missing-module-docstring, unused-variable, undefined-variable, too-many-locals, invalid-name
import os
import sys
import subprocess
import glob
import time


def sourmash_sig(file_name: str, out_file: str, k_size: str, scale: str):
    """
    Compute sourmash signatures for the .fasta files
    :param file_name: path of the file
    :param out_dir: path of the directory where signatures are stored
    :param k_size: size of k-mer, the default value being 31
    :param scale: scaling factor, the default value is 1000

    """
    cmd_sig = "  ".join(["sourmash compute -k ", k_size, " --scaled ", scale, " --track-abundance --name-from-first -o ", out_file, file_name])
    subprocess.call(cmd_sig, shell=True)

def signature_index(k_size: str, file_list_len: str, sig_dir: str):
    """
    Compute index for reference genomes
    :param k_size: size of k-mer, the default value being 31
    :param file_list_len: number of files in signature directory
    :param sig_dir: path of the signature directories
    """
    cmd_index = "  ".join(["sourmash index  -k ", k_size, " genome_index_"+str(file_list_len) , " --traverse-directory  ", sig_dir])
    subprocess.call(cmd_index, shell=True)

def sourmash_gather(w_dir: str, k_size: str, query: str, index_file: str):
    """
    find the best match of query against reference geneomes
    :param w_dir: path of output .csv and .sm files
    :param k_size: size of k-mer, the default value being 31
    :param query: file path of the query signature file
    :param index_file: file path for index of all reference genomes
    """
    query_name = os.path.splitext(query)
    print(query_name)
    csv_out = os.path.join(w_dir, query_name[0]+".csv")
    sm_out = os.path.join(w_dir, query_name[0]+".sm")
    cmd_gather = "  ".join(["sourmash gather -k ", k_size, query, index_file, " -o ", csv_out , " > ", sm_out])
    subprocess.call(cmd_gather, shell=True)
