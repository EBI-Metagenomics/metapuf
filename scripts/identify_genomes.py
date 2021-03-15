# -*- coding: utf-8 -*-
# pylint: disable=line-too-long, missing-module-docstring, unused-variable, undefined-variable, too-many-locals, invalid-name
import os
import sys
import threading
import queue
import time
import argparse
from traceback import format_exc
import pandas as pd
from Bio import SeqIO


def gen_match_list(f_input: str, wdir: str, ref_dir: str):
    """
    Generates file queue for all files in the directory
    :param f_input:  input sample file which is a .csv file
    :param wdir: path for the working directory
    :param ref_dir: path for the parent directory containing all reference genomes

    :return : a list and a file queue
    """
    file_name = os.path.join(wdir, f_input )
    print(file_name)
    data = pd.read_csv(file_name, sep=',')
    unique=set()
    match_dir ={}
    for i in range(len(data)):
        if round(data["f_match"][i], 2) >=0.30:
            unique.add(data["name"][i])
            name_list = list(unique)
    match_list= [item.split() for item in name_list]
    # generate file queue for the reference genome files
    for file in os.listdir(ref_dir):
        filename=os.path.join(ref_dir,file)
        if file.endswith(".fa"):
            print(filename)
            files_queue.put(filename)
    print("files in the queue are {}".format(files_queue))
    return match_list, files_queue

def search(files_queue: list, res_queue: list, match_list: list):
    """
    Performs serach in reference genome files and appends the matches to a list
    :param files_queue: a queue of all files in the reference genome directory
    :param res_queue: a queue storing list of results
    :param match_list: list of matches
    """
    while not files_queue.empty():
        file = files_queue.get(block=False)
        print(file)
        match=[]
        with open(file,'r') as fin:
            record_list = list(SeqIO.parse(fin, "fasta"))
            sequence_id = [str(record.description) for record in record_list]
            sequence_list = [str(record.seq) for record in record_list]
        sequences =dict(zip(sequence_id,sequence_list))
        for i in range(len(match_list)):
            for k,v in sequences.items():
                word= k.split()[0:2]
                if (match_list[i][0] == word[0] or match_list[i][0] == word[1] or match_list[i][0] == word[-1]):
                    match.append(">"+k)
                    match.append(sequences[k])
        res_queue.put(match)

if __name__ == '__main__':
    starttime = time.time()
    program_name = os.path.basename(sys.argv[0])
    try:
        parser = argparse.ArgumentParser(description='Extract genomes from reference database based on matches')
        parser.add_argument('--f_input', required=True, help='enter the f_input filename, which is a csv')
        parser.add_argument('--output', required=True, help='enter the output filename, which is .fa')
        parser.add_argument('--wdir',  required=True, help=' enter location of files for which genome info in required')
        parser.add_argument('--ref_dir', required=True,help='location of the reference database')

        args = vars(parser.parse_args())
        files_queue = queue.Queue()
        results_queue = queue.Queue()
        threads=[]
        genome_matches,  ref_files_queue = gen_match_list(args['f_input'], args['wdir'], args['ref_dir'])
        for _ in range(4):
            thread = threading.Thread(target=search, args=[ref_files_queue,  results_queue, genome_matches])
            thread.start()
            threads.append(thread)
        for thread in threads:
            thread.join()
        final_results = []
        temp=set()
        while not results_queue.empty():
            for item in results_queue.get():
                final_results.append( str(item))
        print("final results are {}".format(final_results))
        result_file=os.path.join(args['wdir'], args['output'])
        with open(result_file,'w') as fout:
            for item in final_results:
                if not item in temp:
                    fout.write(("%s\n" % item))
                    temp.add(item)
                else:
                    continue
        print('That took {} seconds'.format(time.time() - starttime))
    except Exception as error:
        print(program_name + ": " + repr(error) + '\n' + format_exc() + '\n')
        raise error
