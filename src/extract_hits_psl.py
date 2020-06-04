# extract_hits_psl.py
# This script is to extract the highest hits from psl files after BLAT
# Made by: Elfy Ly
# Date: 22 May 2020

import os, fnmatch
import re
from Bio import SearchIO


def create_fout(path_to_fout):
    f = open(path_to_fout, "w+")
    f.write("target_name\tquery_name\tpercent_ID\tPSL_score\n")
    f.close()


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def create_dictionary_start_end(path_to_mapped_contigs):
    dictionary = {}
    with open(path_to_mapped_contigs) as f:
        for line in f:
            (contig_name, contig_start, contig_end) = line.split()
            dictionary[contig_name] = contig_start, contig_end
    print(dictionary['>Contig1'])


# checks if psl file is empty, if no: reads psl file
def read_psl(path_to_psl):
    with open(path_to_psl, 'rt') as myfile:
        nlines = 0
        for every_line in myfile:
            nlines += 1
        if nlines <= 5:
            print(path_to_psl + " is empty")
        else:
            t_name_highest, q_name_highest, id_pct_highest, max_score = extract_hits(path_to_psl)
            write_fout(path_to_fout, t_name_highest, q_name_highest, id_pct_highest, max_score)


# extract hits and returns the t_name, q_name, id_pct and score for highest hit from psl file per contig per species
def extract_hits(path_to_psl):
    qresult = SearchIO.read(path_to_psl, "blat-psl")
    n_hits = len(qresult)

    t_name_highest = ""
    q_name_highest = ""
    id_pct_highest = 0
    max_score = 0
    for hit in range(n_hits):  # loops over hits
        n_hsp_hits = len(qresult[hit])
        for n in range(n_hsp_hits):  # loops over hsp
            blat_hsp = qresult[hit][n]
            if blat_hsp.score > max_score:
                t_name_highest = qresult.id
                q_name_highest = blat_hsp.hit_id
                id_pct_highest = blat_hsp.ident_pct
                max_score = blat_hsp.score

    return t_name_highest, q_name_highest, id_pct_highest, max_score


# add highest BLAT scored exons hits per contig per species in highest_hits.txt file
def write_fout(path_to_fout, t_name_highest, q_name_highest, id_pct_highest, max_score):
    f = open(path_to_fout, "a+")
    f.write(t_name_highest + "\t" + q_name_highest + "\t" + str(id_pct_highest) + "\t" + str(max_score) + "\n")
    f.close()


# Code starts here:
# creates a new empty highest_hits.txt file
path_to_fout = "./results/blat/SRR8528336/highest_hits.txt"
create_fout(path_to_fout)

# moet nog loopen over species
path_to_psl_dir = "./results/blat/SRR8528336/"
list_in_psl_dir = os.listdir(path_to_psl_dir)
list_in_psl_dir_sorted = natural_sort(list_in_psl_dir)
pattern = "*.psl"

# creates dictionary for mapped contigs with their start and end position
path_to_mapped_contigs = './results/assembled_exons/SRR8528336/mapped_contigs.txt'
create_dictionary_start_end(path_to_mapped_contigs)

# loops over psl files and read them one by one
# parse the highest target name, query name, ID percentage and score
# paste them in a file called highest_hits.txt for every species
for entry in list_in_psl_dir_sorted:
    if fnmatch.fnmatch(entry, pattern):
        psl_file = entry
        path_to_psl = path_to_psl_dir + psl_file
        read_psl(path_to_psl)
