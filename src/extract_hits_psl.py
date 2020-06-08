# extract_hits_psl.py
# This script is to extract the highest hits from psl files after BLAT
# Made by: Elfy Ly
# Date: 22 May 2020

import os, fnmatch
import re
from Bio import SearchIO

# moet nog loopen over species
SPECIES = 'SRR8528336'
ID_PCT_CUTOFF = 75
SCORE_CUTOFF = 20


def create_fhits(path_to_fhits):
    f = open(path_to_fhits, "w+")
    f.write("target_name\tquery_name\tpercent_ID\tPSL_score\n")
    f.close()


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def create_dictionary_start_end(path_to_list):
    dictionary = {}
    with open(path_to_list) as f:
        for line in f:
            (name, start, end) = line.split()
            dictionary[name] = start, end
    f.close()
    return dictionary


def check_overlap(dictionary_contigs_sorted, dictionary_exons_sorted):
    for contig in dictionary_contigs_sorted:
        contig_start = int(dictionary_contigs[contig][0])
        contig_end = int(dictionary_contigs[contig][1])

        for exon in dictionary_exons_sorted:
            exon_start = int(dictionary_exons[exon][0])
            exon_end = int(dictionary_exons[exon][1])

            # checks if overlap, if yes then writes contig exon pair in contig_exon_match_list.txt
            if (exon_start <= contig_start <= exon_end) or (exon_start <= contig_end <= exon_end) or \
                    (contig_start <= exon_start <= contig_end) or (contig_start <= exon_end <= contig_end):
                f = open(path_to_contig_exon_match, "a+")
                f.write(contig + '\t' + exon + '\n')
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
            write_fhighest_hits(path_to_fhighest_hits, t_name_highest, q_name_highest, id_pct_highest, max_score)


def check_cutoff(path_to_fhighest_hits, path_to_fhighest_hits_filtered, path_to_fbelow_cutoff):
    with open(path_to_fhighest_hits, 'r') as myfile:
        for line in myfile:
            if not line.startswith('target_name'):
                t_name, q_name, id_pct, score = line.split("\t")
                if float(id_pct) >= ID_PCT_CUTOFF and float(score) > SCORE_CUTOFF:
                    filtered_hits_file = open(path_to_fhighest_hits_filtered, "a+")
                    filtered_hits_file.write(line)
                    filtered_hits_file.close()
                else:
                    fbelow_cutoff = open(path_to_fbelow_cutoff, "a+")
                    fbelow_cutoff.write(t_name + "\t" + q_name)
                    fbelow_cutoff.close()


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
def write_fhighest_hits(path_to_fhighest_hits, t_name_highest, q_name_highest, id_pct_highest, max_score):
    f = open(path_to_fhighest_hits, "a+")
    f.write(t_name_highest + "\t" + q_name_highest + "\t" + str(id_pct_highest) + "\t" + str(max_score) + "\n")
    f.close()


# Create directory of given path if it doesn't exist
def create_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        print("Directory ", path, " Created ")
    else:
        print("Directory ", path, " already exists")


# # creates new file
# def create_ftxt(path):
#     txt_file = open(path + '.txt', "w+")
#     txt_file.close()
#     print(path + ".txt is created")


# checks if contig-exon pairs in highest_hits_filtered.txt are present in contig_exon_match_list.txt
# if yes, creates new exon_name.fasta file with contig consensus sequence for MAFFT
def check_overlap_hit_pairs(dictionary_hits_sorted, dictionary_match_sorted, dictionary_hits, dictionary_match,
                            path_to_mafft_species_dir, path_to_consensus_dir, list_consensus_dir, path_to_fseq_exons,
                            path_to_fno_match_pairs):
    for contig_name_fhit in dictionary_hits_sorted:
        for contig_name_fmatch in dictionary_match_sorted:
            if contig_name_fhit == contig_name_fmatch:
                exon_name_fhit = dictionary_hits[contig_name_fhit].strip()
                exon_name_fmatch = dictionary_match[contig_name_fmatch].strip()
                if exon_name_fhit == exon_name_fmatch:
                    path_to_mafft_fexon = path_to_mafft_species_dir + exon_name_fhit
                    # create_ftxt(path_to_mafft_fexon)
                    # read and write matching consensus sequence in new text file
                    write_ffasta(list_consensus_dir, contig_name_fhit, path_to_consensus_dir, path_to_mafft_fexon)
                    # write all contig-exon pairs in fseq_exons
                    fseq_exons = open(path_to_fseq_exons, "a+")
                    fseq_exons.write(contig_name_fhit + "\t" + exon_name_fhit + "\n")
                    fseq_exons.close()
                else:
                    # print pair die niet present is in path_to_fno_match_pairs
                    fno_match_pairs = open(path_to_fno_match_pairs, "a+")
                    fno_match_pairs.write(contig_name_fhit + "\t" + exon_name_fhit)
                    fno_match_pairs.close()


def write_ffasta(list_consensus_dir, contig_name_fhit, path_to_consensus_dir, path_to_mafft_fexon):
    for contig_name_file in list_consensus_dir:
        contig_name, txt = contig_name_file.split('.')

        if contig_name_fhit == contig_name:
            contig_consensus_file = open(path_to_consensus_dir + contig_name + '.txt', 'rt')
            exon_file = open(path_to_mafft_fexon, "a+")
            for line in contig_consensus_file:
                exon_file.write(line)
            contig_consensus_file.close()
            exon_file.close()


# # # Code starts here:
# creates a new empty highest_hits.txt file
path_to_fhighest_hits = "./results/blat/" + SPECIES + "/highest_hits.txt"
create_fhits(path_to_fhighest_hits)

# creates separate dictionary for mapped contigs with their start and end position
path_to_mapped_contigs = "./results/assembled_exons/" + SPECIES + "/mapped_contigs.txt"
dictionary_contigs = create_dictionary_start_end(path_to_mapped_contigs)
dictionary_contigs_sorted = natural_sort(dictionary_contigs)

# creates separate dictionary for target exons with their start and end position
path_to_exons_enum = "./data/exons/AT_exon_enum.txt"
dictionary_exons = create_dictionary_start_end(path_to_exons_enum)
dictionary_exons_sorted = natural_sort(dictionary_exons)

# create an empty contig_exon_match_list.txt
path_to_contig_exon_match = "./results/assembled_exons/" + SPECIES + "/contig_exon_match_list.txt"
f = open(path_to_contig_exon_match, "w+")
f.close()

# checks which mapped contig starting and ending position matches with exons starting and ending position
# and puts the matching pairs in contig_exon_match_list.txt
check_overlap(dictionary_contigs_sorted, dictionary_exons_sorted)

# loops over psl files in the psl dir
path_to_psl_dir = "./results/blat/" + SPECIES + "/"
list_in_psl_dir = os.listdir(path_to_psl_dir)
list_in_psl_dir_sorted = natural_sort(list_in_psl_dir)
pattern = "*.psl"

# create new file for filtered highest hits based on cutoffs
path_to_fhighest_hits_filtered = path_to_psl_dir + "highest_hits_filtered.txt"
create_fhits(path_to_fhighest_hits_filtered)

# read the psl files one by one
# parse the highest target name, query name, ID percentage and score
# paste them in a file called highest_hits.txt for every species
for entry in list_in_psl_dir_sorted:
    if fnmatch.fnmatch(entry, pattern):
        psl_file = entry
        path_to_psl = path_to_psl_dir + psl_file
        read_psl(path_to_psl)

# create file with matches with PID/score < cutoffs
path_to_fbelow_cutoff = path_to_psl_dir + "below_cutoff_pairs.txt"
# create_ftxt(path_to_fbelow_cutoff)
check_cutoff(path_to_fhighest_hits, path_to_fhighest_hits_filtered, path_to_fbelow_cutoff)

# create dictionary for highest_hits_filtered contig_exon match pairs
with open(path_to_fhighest_hits_filtered, 'rt') as myfile:
    dictionary_hits = {}
    for line in myfile:
        if not line.startswith('target_name'):
            (contig_name, exon_name, percent_ID, PSL_score) = line.split('\t')
            dictionary_hits[contig_name] = exon_name
dictionary_hits_sorted = natural_sort(dictionary_hits)

# create dictionary for contig_exon_match pairs
with open(path_to_contig_exon_match, 'rt') as myfile:
    dictionary_match = {}
    for line in myfile:
        (contig_name, exon_name) = line.split('\t')
        dictionary_match[contig_name] = exon_name
dictionary_match_sorted = natural_sort(dictionary_match)

# create new MAFFT dirs for input and stats
path_to_mafft = './results/mafft/'
path_to_mafft_species_dir = path_to_mafft + SPECIES + "/"
path_to_mafft_stats = path_to_mafft_species_dir + "stats/"
create_dir(path_to_mafft)
create_dir(path_to_mafft_species_dir)
create_dir(path_to_mafft_stats)

# create statistical files
path_to_fseq_exons = path_to_mafft_stats + "sequenced_exons.txt"
path_to_fno_match_pairs = path_to_mafft_stats + "no_match_pairs.txt"
# create_ftxt(path_to_fseq_exons)
# create_ftxt(path_to_fno_match_pairs)

# prepare fasta files for MAFFT
path_to_consensus_dir = "./results/consensus/" + SPECIES + "/"
list_consensus_dir = natural_sort(os.listdir(path_to_consensus_dir))
check_overlap_hit_pairs(dictionary_hits_sorted, dictionary_match_sorted, dictionary_hits, dictionary_match,
                        path_to_mafft_species_dir, path_to_consensus_dir, list_consensus_dir,
                        path_to_fseq_exons, path_to_fno_match_pairs)
