# extract_hits_psl.py
# This script is to extract the highest hits from psl files after BLAT
# Made by: Elfy Ly
# Date: 22 May 2020

import os, fnmatch
import re
from Bio import SearchIO

PSL_HEADER_LINES = 5
ID_PCT_CUTOFF = 75
SCORE_CUTOFF = 20


def create_fhits(path_to_fhits):
    f = open(path_to_fhits, "w+")
    f.write("target_name\tquery_name\tpercent_ID\tPSL_score\n")
    f.close()


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def create_ftxt(path):
    f = open(path, "w+")
    f.close()
    print(path + " is created")


def create_dictionary_start_end(path_to_list):
    dictionary = {}
    with open(path_to_list) as f:
        for line in f:
            (name, start, end) = line.split()
            dictionary[name] = start, end
    f.close()
    return dictionary


def check_overlap(dictionary_contigs_lsorted, dictionary_exons_lsorted, path_to_contig_exon_match):
    for contig in dictionary_contigs_lsorted:
        contig_start = int(dictionary_contigs[contig][0])
        contig_end = int(dictionary_contigs[contig][1])

        for exon in dictionary_exons_lsorted:
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
        if nlines <= PSL_HEADER_LINES:
            print(path_to_psl + " is empty")
        else:
            t_name_highest, q_name_highest, id_pct_highest, max_score = extract_hits(path_to_psl)
            write_fhighest_hits(path_to_fhighest_hits, t_name_highest, q_name_highest, id_pct_highest, max_score)


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
                t_name_highest = qresult.id
                q_name_highest = blat_hsp.hit_id
                id_pct_highest = blat_hsp.ident_pct
                max_score = blat_hsp.score

    return t_name_highest, q_name_highest, id_pct_highest, max_score


# add highest BLAT scored exons hits per contig per species_name in highest_hits.txt file
def write_fhighest_hits(path_to_fhighest_hits, t_name_highest, q_name_highest, id_pct_highest, max_score):
    f = open(path_to_fhighest_hits, "a+")
    f.write(t_name_highest + "\t" + q_name_highest + "\t" + str(id_pct_highest) + "\t" + str(max_score) + "\n")
    f.close()


def check_cutoff(path_to_fhighest_hits, path_to_fhighest_hits_filtered, path_to_fbelow_cutoff):
    with open(path_to_fhighest_hits, 'r') as myfile:
        for line in myfile:
            if not line.startswith('target_name'):
                t_name_highest, q_name_highest, id_pct_highest, max_score = line.split("\t")
                if float(id_pct_highest) >= ID_PCT_CUTOFF and float(max_score) > SCORE_CUTOFF:
                    write_line(path_to_fhighest_hits_filtered, line)
                else:
                    write_line(path_to_fbelow_cutoff, line)


def write_line(path, line):
    f = open(path, "a+")
    f.write(line)
    f.close()


# create dictionary for highest_hits_filtered contig_exon match pairs
def create_dictionary_highest_hits_filtered(path_to_fhighest_hits_filtered):
    dictionary_hits = {}
    with open(path_to_fhighest_hits_filtered, 'rt') as myfile:
        for line in myfile:
            if not line.startswith('target_name'):
                (contig_name, exon_name, percent_ID, PSL_score) = line.split('\t')
                dictionary_hits[contig_name] = exon_name
    return dictionary_hits


# create dictionary for contig_exon_match pairs
def create_dictionary_match_pairs(path_to_contig_exon_match):
    dictionary_match = {}
    with open(path_to_contig_exon_match, 'rt') as myfile:
        for line in myfile:
            (contig_name, exon_name) = line.split('\t')
            dictionary_match[contig_name] = exon_name
    return dictionary_match


# Create directory of given path if it doesn't exist
def create_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        print("Directory ", path, " Created ")
    else:
        print("Directory ", path, " already exists")


def write_fstat(path, contig_name, exon_name):
    f = open(path, "a+")
    f.write(contig_name + "\t" + exon_name + "\n")
    f.close()


# if match: creates new exon_name.fasta file for MAFFT
# and write pair in sequenced_exons.txt. If not a match: write pair in no_match_pairs.txt
def create_mafft_ffasta(dictionary_hits_lsorted, dictionary_match_lsorted, dictionary_hits, dictionary_match,
                        path_to_mafft_species_dir, path_to_fseq_exons, path_to_fno_match):
    for hit_contig in dictionary_hits_lsorted:
        for match_contig in dictionary_match_lsorted:
            if hit_contig == match_contig:
                hit_exon = dictionary_hits[hit_contig].strip()
                match_exon = dictionary_match[match_contig].strip()
                if hit_exon == match_exon:
                    path_to_fexon = path_to_mafft_species_dir + hit_exon + '.txt'
                    create_ftxt(path_to_fexon)
                    write_fstat(path_to_fseq_exons, hit_contig, hit_exon)
                else:
                    write_fstat(path_to_fno_match, hit_contig, hit_exon)


def append_mafft_ffasta(dictionary_hits_lsorted, dictionary_match_lsorted, dictionary_hits, dictionary_match,
                    path_to_mafft_species_dir, sorted_list_consensus_dir, path_to_consensus_dir):
    for hit_contig in dictionary_hits_lsorted:
        for match_contig in dictionary_match_lsorted:
            if hit_contig == match_contig:
                hit_exon = dictionary_hits[hit_contig].strip()
                match_exon = dictionary_match[match_contig].strip()
                if hit_exon == match_exon:
                    path_to_fexon = path_to_mafft_species_dir + hit_exon + '.txt'
                    for contig_name_file in sorted_list_consensus_dir:
                        contig_name, txt = contig_name_file.split('.')
                        if hit_contig == contig_name:
                            path_to_fcontig_consensus = path_to_consensus_dir + contig_name + '.txt'
                            contig_consensus_file = open(path_to_fcontig_consensus, 'rt')
                            exon_file = open(path_to_fexon, "a+")
                            for line in contig_consensus_file:
                                exon_file.write(line)
                            contig_consensus_file.close()
                            exon_file.close()


# # # Code starts here:
path_to_assembled_exons_dir = "results/assembled_exons/"
dirs = os.listdir(path_to_assembled_exons_dir)
for species_name in dirs:

    # # STEP 1: Checks all fragment overlaps between contigs and exons

    # creates separate dictionary for mapped contigs and target exons with their start and end position
    path_to_assembled_species_dir = "results/assembled_exons/" + species_name + "/"
    path_to_mapped_contigs = path_to_assembled_species_dir + "mapped_contigs.txt"
    path_to_exons_enum = "./data/exons/AT_exon_enum.txt"
    dictionary_contigs = create_dictionary_start_end(path_to_mapped_contigs)
    dictionary_exons = create_dictionary_start_end(path_to_exons_enum)
    dictionary_contigs_lsorted = natural_sort(dictionary_contigs)
    dictionary_exons_lsorted = natural_sort(dictionary_exons)

    # add matching pairs in contig_exon_match_list.txt if mapped contig fragments matches/overlaps with the exon fragments
    path_to_blat_species_dir = "./results/blat/" + species_name + "/"
    path_to_blat_stats = path_to_blat_species_dir + "stats/"
    create_dir(path_to_blat_stats)

    path_to_contig_exon_match = path_to_blat_stats + "contig_exon_match_list.txt"
    create_ftxt(path_to_contig_exon_match)
    check_overlap(dictionary_contigs_lsorted, dictionary_exons_lsorted, path_to_contig_exon_match)


    # # STEP 2: Reads all PSL files one by one and parse info in highest_hits.txt for every contig
    path_to_fhighest_hits = path_to_blat_stats + "highest_hits.txt"
    create_fhits(path_to_fhighest_hits)

    list_in_psl_dir = os.listdir(path_to_blat_species_dir)
    list_in_psl_dir_sorted = natural_sort(list_in_psl_dir)
    pattern = "*.psl"
    for file in list_in_psl_dir_sorted:
        if fnmatch.fnmatch(file, pattern):
            psl_file = file
            path_to_psl = path_to_blat_species_dir + psl_file
            read_psl(path_to_psl)


    # # STEP 3: Filter on cutoffs for score and ID%
    # if > cutoff: paste information in highest_hits_filtered.txt; if < cutoff: copy in below_cutoff_pairs.txt
    path_to_fhighest_hits_filtered = path_to_blat_stats + "highest_hits_filtered.txt"
    path_to_fbelow_cutoff = path_to_blat_stats + "below_cutoff_pairs.txt"
    create_fhits(path_to_fhighest_hits_filtered)
    create_ftxt(path_to_fbelow_cutoff)
    check_cutoff(path_to_fhighest_hits, path_to_fhighest_hits_filtered, path_to_fbelow_cutoff)


    # # STEP 4: Check if highest_hits_filtered.txt matches are also overlapped fragments in contig_exon_match_list.txt
    # # If yes: creates new .fasta file for MAFFT

    # create dictionaries for highest_hits_filtered.txt and contig_exon_match_list.txt
    dictionary_hits = create_dictionary_highest_hits_filtered(path_to_fhighest_hits_filtered)
    dictionary_match = create_dictionary_match_pairs(path_to_contig_exon_match)
    dictionary_hits_lsorted = natural_sort(dictionary_hits)
    dictionary_match_lsorted = natural_sort(dictionary_match)

    # create new MAFFT dir for input
    path_to_mafft = './results/mafft/'
    path_to_mafft_species_dir = path_to_mafft + species_name + "/"
    create_dir(path_to_mafft)
    create_dir(path_to_mafft_species_dir)

    # creates new stats directory and stats files sequenced_exons.txt and no_match_pairs.txt
    path_to_mafft_stats = path_to_mafft_species_dir + "stats/"
    create_dir(path_to_mafft_stats)
    path_to_fseq_exons = path_to_mafft_stats + "sequenced_exons.txt"
    path_to_fno_match = path_to_mafft_stats + "no_match_pairs.txt"
    create_ftxt(path_to_fseq_exons)
    create_ftxt(path_to_fno_match)

    create_mafft_ffasta(dictionary_hits_lsorted, dictionary_match_lsorted, dictionary_hits, dictionary_match,
                        path_to_mafft_species_dir, path_to_fseq_exons, path_to_fno_match)

    # append all contig consensus sequence to the exon_name.fasta file
    path_to_consensus_dir = "./results/consensus/" + species_name + "/"
    list_consensus_dir = os.listdir(path_to_consensus_dir)
    sorted_list_consensus_dir = natural_sort(list_consensus_dir)

    append_mafft_ffasta(dictionary_hits_lsorted, dictionary_match_lsorted, dictionary_hits, dictionary_match,
                        path_to_mafft_species_dir, sorted_list_consensus_dir, path_to_consensus_dir)
