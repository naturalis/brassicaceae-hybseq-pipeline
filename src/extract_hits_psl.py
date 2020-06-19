# extract_hits_psl.py
# This script is to extract the highest hits from psl files after BLAT
# Can be executed by: $ python extract_hits_psl.py [SAMPLE_NAME], for example: $ python extract_hits_psl.py SRR8528336
# Made by: Elfy Ly
# Date: 22 May 2020

import sys
import os, fnmatch
import re
from Bio import SearchIO
from Bio import SeqIO

SAMPLE_NAME = sys.argv[1]
PSL_HEADER_LINES = 5
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


# add highest BLAT scored exons hits per contig per SAMPLE_NAME in highest_hits.txt file
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


# check match exon pairs and if match: creates new exon_name.fasta file for MAFFT input and append original exon
# sequence and contig if file does not exist yet. If already exist: only appends contigs.
# Writes all exons with multiple contigs in multiple_contigs.txt.
# and write pair in sequenced_exons.txt. If not a match: write pair in no_match_pairs.txt
def create_exon_ffasta(dictionary_hits_lsorted, dictionary_match_lsorted, dictionary_hits, dictionary_match,
                       path_to_mapped_exons_species_dir, path_to_fseq_exons, path_to_fexons_seq,
                       sorted_list_consensus_dir, path_to_consensus_dir, path_to_fmultiple_contigs, path_to_fno_match):
    temporary_exon = ""
    for hit_contig in dictionary_hits_lsorted:
        for match_contig in dictionary_match_lsorted:
            if hit_contig == match_contig:
                hit_exon = dictionary_hits[hit_contig].strip()
                match_exon = dictionary_match[match_contig].strip()
                if hit_exon == match_exon:
                    path_to_fexon = path_to_mapped_exons_species_dir + hit_exon + '.fasta'
                    if match_exon != temporary_exon:
                        create_ftxt(path_to_fexon)
                        write_fstat(path_to_fseq_exons, hit_contig, hit_exon)
                        for record in SeqIO.parse(path_to_fexons_seq, "fasta"):
                            if hit_exon == record.id:
                                fexon = open(path_to_fexon, "a+")
                                fexon.write(">" + record.id + "\n")
                                fexon.write(str(record.seq) + "\n")
                        append_seq(sorted_list_consensus_dir, hit_contig, path_to_consensus_dir, path_to_fexon)
                        temporary_exon = match_exon
                    elif match_exon == temporary_exon:
                        append_seq(sorted_list_consensus_dir, hit_contig, path_to_consensus_dir, path_to_fexon)
                        write_fstat(path_to_fmultiple_contigs, hit_contig, hit_exon)
                else:
                    write_fstat(path_to_fno_match, hit_contig, hit_exon)


def append_seq(sorted_dir_list, match_name, path_to_seq_dir, path_to_fexon):
    for name_file in sorted_dir_list:
        name, txt = name_file.split('.')
        if match_name == name:
            path_to_fseq = path_to_seq_dir + name + '.txt'
            fseq = open(path_to_fseq, 'rt')
            exon_file = open(path_to_fexon, "a+")
            for line in fseq:
                exon_file.write(line)
            fseq.close()
            exon_file.close()


def create_YAML(path):
    f = open(path, "w+")
    f.write("exons:\n")
    f.close()
    print(path + " is created")


# # # Code starts here:
path_to_mapped_contigs_dir = "results/4_mapped_contigs/"

'''STEP 1: Checks all fragment overlaps between contigs and exons'''
# creates separate dictionary for mapped contigs and target exons with their start and end position
path_to_assembled_species_dir = path_to_mapped_contigs_dir + SAMPLE_NAME + "/"
path_to_fmapped_contigs = path_to_assembled_species_dir + "mapped_contigs.txt"
path_to_exons_enum = "./data/exons/AT_exon_enum.txt"
dictionary_contigs = create_dictionary_start_end(path_to_fmapped_contigs)
dictionary_exons = create_dictionary_start_end(path_to_exons_enum)
dictionary_contigs_lsorted = natural_sort(dictionary_contigs)
dictionary_exons_lsorted = natural_sort(dictionary_exons)

# add matching pairs in contig_exon_match_list.txt if mapped contig fragments matches/overlaps
# with the exon fragments
path_to_blat_species_dir = "./results/6_identified_contigs_blat/" + SAMPLE_NAME + "/"
path_to_blat_stats = path_to_blat_species_dir + "stats/"
create_dir(path_to_blat_stats)

path_to_contig_exon_match = path_to_blat_stats + "contig_exon_match_list.txt"
create_ftxt(path_to_contig_exon_match)
check_overlap(dictionary_contigs_lsorted, dictionary_exons_lsorted, path_to_contig_exon_match)

'''STEP 2: Reads all PSL files one by one and parse info in highest_hits.txt for every contig'''
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

'''STEP 3: Filter on cutoffs for score and ID%'''
# if > cutoff: paste information in highest_hits_filtered.txt; if < cutoff: copy in below_cutoff_pairs.txt
path_to_fhighest_hits_filtered = path_to_blat_stats + "highest_hits_filtered.txt"
path_to_fbelow_cutoff = path_to_blat_stats + "below_cutoff_pairs.txt"
create_fhits(path_to_fhighest_hits_filtered)
create_ftxt(path_to_fbelow_cutoff)
check_cutoff(path_to_fhighest_hits, path_to_fhighest_hits_filtered, path_to_fbelow_cutoff)

'''STEP 4: Check if highest_hits_filtered.txt matches are also overlapped fragments in contig_exon_match_list.txt
If yes: creates new .fasta file for MAFFT input'''
# create dictionaries for highest_hits_filtered.txt and contig_exon_match_list.txt
dictionary_hits = create_dictionary_highest_hits_filtered(path_to_fhighest_hits_filtered)
dictionary_match = create_dictionary_match_pairs(path_to_contig_exon_match)
dictionary_hits_lsorted = natural_sort(dictionary_hits)
dictionary_match_lsorted = natural_sort(dictionary_match)

# create new mapped_exons dir for input
path_to_mapped_exons = './results/7_mapped_exons/'
path_to_mapped_exons_species_dir = path_to_mapped_exons + SAMPLE_NAME + "/"
create_dir(path_to_mapped_exons)
create_dir(path_to_mapped_exons_species_dir)

# creates new stats directory and stats files sequenced_exons.txt and no_match_pairs.txt
path_to_mapped_exons_stats = path_to_mapped_exons_species_dir + "stats/"
create_dir(path_to_mapped_exons_stats)
path_to_fseq_exons = path_to_mapped_exons_stats + "sequenced_exons.txt"
path_to_fmultiple_contigs = path_to_mapped_exons_stats + "multiple_contigs.txt"
path_to_fno_match = path_to_mapped_exons_stats + "no_match_pairs.txt"
create_ftxt(path_to_fseq_exons)
create_ftxt(path_to_fno_match)
create_ftxt(path_to_fmultiple_contigs)

# create path to append original exon and consensus sequences to the exon fasta files as MAFFT input
path_to_fexons_seq = "./data/exons/exons_AT.fasta"
path_to_consensus_dir = "./results/5_consensus_contigs/" + SAMPLE_NAME + "/"
list_consensus_dir = os.listdir(path_to_consensus_dir)
sorted_list_consensus_dir = natural_sort(list_consensus_dir)

create_exon_ffasta(dictionary_hits_lsorted, dictionary_match_lsorted, dictionary_hits, dictionary_match,
                   path_to_mapped_exons_species_dir, path_to_fseq_exons, path_to_fexons_seq,
                   sorted_list_consensus_dir, path_to_consensus_dir, path_to_fmultiple_contigs, path_to_fno_match)

'''STEP 5: Creates configuration files for MAFFT'''
# creates configuration YAML file in envs MAFFT dir
path_to_MAFFT_configs = "./envs/MAFFT/"
create_dir(path_to_MAFFT_configs)
path_to_fYAML = path_to_MAFFT_configs + SAMPLE_NAME + ".yaml"
create_YAML(path_to_fYAML)

with open(path_to_fseq_exons, "rt") as fseq_exons:
    for line in fseq_exons:
        contig_name, exon_name = line.split("\t")
        exon_name = exon_name.strip()
        fYAML = open(path_to_fYAML, "a+")
        fYAML.write("    " + exon_name + ": " + path_to_mapped_exons_species_dir + exon_name + ".fasta\n")
        fYAML.close()
fseq_exons.close()
