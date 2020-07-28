# extract_hits_psl_seq_genomes.py
# This script is to extract the highest hits from psl files after BLAT on the sequence genomes
# Can be executed by: $ python extract_hits_psl_seq_genomes.py [REF_NAME],
# for example: $ python extract_hits_psl_seq_genomes.py ref-tla
# Made by: Elfy Ly
# Date: 14 July 2020

import sys
import os
import fnmatch
import re
from Bio import SearchIO

PSL_HEADER_LINES = 5
SEQ_GENOME = sys.argv[1]
# At least half of the target and query have to be matched
TCUTOFF = 0.5
QCUTOFF = 0.5


# Create directory of given path if it doesn't exist
def create_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        print("Directory ", path, " Created ")
    else:
        print("Directory ", path, " already exists")


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def create_fhits(path_to_fhits):
    f = open(path_to_fhits, "w+")
    f.write("target_name\tquery_name\tpercent_ID\tPSL_score\tqratio\ttratio\ttstart\ttend\n")
    f.close()


# checks if psl file is empty, if no: reads psl file
def read_psl(path_to_psl, nempty, nexon_match, nexf):
    with open(path_to_psl, 'rt') as myfile:
        nlines = 0
        for every_line in myfile:
            nlines += 1
        if nlines <= PSL_HEADER_LINES:
            print(path_to_psl + " is empty")
            nempty += 1
        else:
            nexon_match += 1
            t_name_highest, q_name_highest, id_pct_highest, max_score, qsize, tsize, tstart, tend, tbase = \
                extract_hits(path_to_psl)
            qratio = float(max_score) / qsize
            tratio = float(max_score) / tsize
            print("check if tbase: " + str(tbase) + "is smaller than qsize: " + str(qsize))

            if (tratio > TCUTOFF) and (qratio > QCUTOFF) and (tbase < qsize):
                nexf += 1
                write_fhighest_hits(path_to_fhighest_hits, t_name_highest, q_name_highest, id_pct_highest, max_score,
                                    qratio, tratio, tstart, tend)
            else:
                print("does not satisfy cutoff")

        return nempty, nexon_match, nexf


# extract hits and returns the t_name, q_name, id_pct and score for highest hit from psl file per contig per species
def extract_hits(path_to_psl):
    qresult = SearchIO.read(path_to_psl, "blat-psl")
    n_hits = len(qresult)

    t_name_highest = ""
    q_name_highest = ""
    id_pct_highest = 0
    max_score = 0
    for hit in range(n_hits):
        n_hsp_hits = len(qresult[hit])
        for n in range(n_hsp_hits):
            blat_hsp = qresult[hit][n]
            if blat_hsp.score > max_score:
                t_name_highest = blat_hsp.hit_id
                q_name_highest = qresult.id
                id_pct_highest = blat_hsp.ident_pct
                max_score = blat_hsp.score
                qsize = blat_hsp.query_span
                tsize = blat_hsp.hit_span
                tstart = blat_hsp.hit_start
                tend = blat_hsp.hit_end
                tbase = blat_hsp.hit_gap_num

    return t_name_highest, q_name_highest, id_pct_highest, max_score, qsize, tsize, tstart, tend, tbase


# add highest BLAT scored exons hits per contig per SAMPLE_NAME in highest_hits.txt file
def write_fhighest_hits(path_to_fhighest_hits, t_name_highest, q_name_highest, id_pct_highest, max_score, qratio,
                        tratio, tstart, tend):
    f = open(path_to_fhighest_hits, "a+")
    f.write(t_name_highest + "\t" + q_name_highest + "\t" + str(id_pct_highest) + "\t" + str(max_score) + "\t" +
            str("%.2f" % qratio) + "\t" + str("%.2f" % tratio) + "\t" + str(tstart) + "\t" + str(tend) + "\n")
    f.close()


# Code starts here
path_to_blat_dir = "results/B01_identified_contigs_blat/" + SEQ_GENOME + "/"
path_to_blat_stats = path_to_blat_dir + "stats/"
create_dir(path_to_blat_stats)

path_to_fhighest_hits = path_to_blat_stats + "highest_hits.txt"
create_fhits(path_to_fhighest_hits)


list_in_psl_dir = os.listdir(path_to_blat_dir)
list_in_psl_dir_sorted = natural_sort(list_in_psl_dir)
pattern = "*.psl"
nexf = 0  # aantal exonen filtered
nexon_match = 0
nempty = 0
for file in list_in_psl_dir_sorted:
    if fnmatch.fnmatch(file, pattern):
        psl_file = file
        path_to_psl = path_to_blat_dir + psl_file
        nexon_match, nempty, nexf = read_psl(path_to_psl, nexon_match, nempty, nexf)
print("there are " + str(nempty) + " empty files and " + str(nexon_match) + " exons with match")
print("after filtering only " + str(nexf) + " are left")


