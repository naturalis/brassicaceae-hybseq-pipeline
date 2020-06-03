# extract_hits_psl.py
# This script is to extract the highest hits from psl files after BLAT
# Made by: Elfy Ly
# Date: 22 May 2020

import os, fnmatch
from Bio import SearchIO


def create_fout(path_to_fout):
    f = open(path_to_fout, "w+")
    f.write("target_name\tquery_name\tpercent_ID\tPSL_score\n")
    f.close()


# reads psl file and returns the t_name, q_name, id_pct and score for highest hit per contig
def read_psl(path_to_psl):
    # checks if empty
    nlines = 0
    with open(path_to_psl, 'rt') as myfile:
        for myline in myfile:
            nlines += 1
        if nlines > 1:
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

            # highest BLAT hits are written per contig per species
            write_fout(path_to_fout, t_name_highest, q_name_highest, id_pct_highest, max_score)
        else:
            print(path_to_psl + " is empty")


# add highest scored exons per contig in highest_hits.txt file
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
pattern = "*.psl"

# loops over psl files and read them one by one
# parse the highest target name, query name, ID percentage and score
# paste them in a file called highest_hits.txt for every species
for entry in list_in_psl_dir:
    if fnmatch.fnmatch(entry, pattern):
        psl_file = entry
        path_to_psl = path_to_psl_dir + psl_file
        read_psl(path_to_psl)
