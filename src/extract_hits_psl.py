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


def create_dictionary_start_end(path_to_list, dictionary):
    with open(path_to_list) as f:
        for line in f:
            (name, start, end) = line.split()
            dictionary[name] = start, end
    f.close()
    return dictionary


def check_overlap(sorted_dictionary_contigs, sorted_dictionary_exons):
    for contig in sorted_dictionary_contigs:
        contig_start = dictionary_contigs[contig][0]
        contig_end = dictionary_contigs[contig][1]

        for exon in sorted_dictionary_exons:
            exon_start = dictionary_exons[exon][0]
            exon_end = dictionary_exons[exon][1]

            # checks if overlap, if yes then writes contig exon pair in contig_exon_match_list.txt
            if exon_start <= contig_start <= exon_end or exon_start <= contig_end <= exon_end or \
                    contig_start <= exon_start <= contig_end or contig_start <= exon_end <= contig_end:
                f = open(path_to_contig_exon_match, "a+")
                f.write(contig + ' ' + exon + '\n')
                f.close()


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


# # # Code starts here:
# creates a new empty highest_hits.txt file
path_to_fout = "./results/blat/SRR8528336/highest_hits.txt"
create_fout(path_to_fout)

# creates separate dictionary for mapped contigs with their start and end position
path_to_mapped_contigs = './results/assembled_exons/SRR8528336/mapped_contigs.txt'
dictionary_contigs = {}
dictionary_contigs = create_dictionary_start_end(path_to_mapped_contigs, dictionary_contigs)
sorted_dictionary_contigs = natural_sort(dictionary_contigs)

# creates separate dictionary for target exons with their start and end position
path_to_exons_enum = './data/exons/AT_exon_enum.txt'
dictionary_exons = {}
dictionary_exons = create_dictionary_start_end(path_to_exons_enum, dictionary_exons)
sorted_dictionary_exons = natural_sort(dictionary_exons)

# create an empty contig_exon_match_list.txt
path_to_contig_exon_match = "./results/assembled_exons/SRR8528336/contig_exon_match_list.txt"
f = open(path_to_contig_exon_match, "w+")
f.close()

# checks which mapped contig starting and ending position matches with exons starting and ending position
# and puts the matching pairs in contig_exon_match_list.txt
check_overlap(sorted_dictionary_contigs, sorted_dictionary_exons)

# loops over psl files in the psl dir
# moet nog loopen over species
path_to_psl_dir = "./results/blat/SRR8528336/"
list_in_psl_dir = os.listdir(path_to_psl_dir)
list_in_psl_dir_sorted = natural_sort(list_in_psl_dir)
pattern = "*.psl"

# read the psl files one by one
# parse the highest target name, query name, ID percentage and score
# paste them in a file called highest_hits.txt for every species
for entry in list_in_psl_dir_sorted:
    if fnmatch.fnmatch(entry, pattern):
        psl_file = entry
        path_to_psl = path_to_psl_dir + psl_file
        read_psl(path_to_psl)


