# extract_hits_psl.py
# This script is to extract the highest hits from psl files after BLAT
# Can be executed by: $ python extract_hits_psl.py [SAMPLE], for example: $ python extract_hits_psl.py SRR8528336
# Made by: Elfy Ly
# Date: 22 May 2020

import sys
import os
import fnmatch
import re
from Bio import SearchIO
from Bio import SeqIO

PSL_FILE = sys.argv[1]
PSL_HEADER_LINES = 5
ID_PCT_CUTOFF = 75
SCORE_CUTOFF = 20


# Create directory of given path if it doesn't exist
def create_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        print("Directory ", path, " Created ")
    else:
        print("Directory ", path, " already exists")


def create_ftxt(path):
    if not os.path.exists(path):
        f = open(path, "w+")
        f.close()
        print(path + " is created")
    else:
        print("File ", path, " already exists")


def create_dictionary_start_end(path_to_list):
    dictionary = {}
    with open(path_to_list) as f:
        for line in f:
            (name, start, end) = line.split()
            dictionary[name] = start, end
    f.close()
    return dictionary


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def get_start_end_contig(path_fmapped_contigs, contig_nr):
    f = open(path_to_fmapped_contigs, "rt")
    search_contig = "Contig" + contig_nr + "\t"
    for line in f:
        if search_contig in line:
            contig_code, contig_start, contig_end = line.split("\t")
            return contig_start, contig_end


def check_overlap(contig_start, contig_end, dictionary_exons_lsorted, path_to_contig_exon_match, contig_nr):
    for exon in dictionary_exons_lsorted:
        exon_start = int(dictionary_exons[exon][0])
        exon_end = int(dictionary_exons[exon][1])
        contig_start = int(contig_start)
        contig_end = int(contig_end)

        contig_name = "Contig" + contig_nr
        # checks if overlap, if yes then writes contig exon pair in contig_exon_match_list.txt
        if (exon_start <= contig_start <= exon_end) or (exon_start <= contig_end <= exon_end) or \
                (contig_start <= exon_start <= contig_end) or (contig_start <= exon_end <= contig_end):
            f = open(path_to_contig_exon_match, "a+")
            f.write(contig_name + '\t' + exon + '\n')
            f.close()


def string_in_file(file_name, string_to_search):
    """ Check if any line in the file contains given string """
    with open(file_name, 'r') as read_obj:
        for line in read_obj:
            if string_to_search in line:
                return True
    return False


def create_fhits(path_to_fhits):
    if not os.path.exists(path_to_fhits):
        f = open(path_to_fhits, "w+")
        f.write("target_name\tquery_name\tpercent_ID\tPSL_score\n")
        f.close()
    else:
        print(path_to_fhits + " already exists")


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


# extract hits and returns the t_name, q_name, id_pct and score for highest hit from psl file
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


# add highest BLAT scored exons hits per contig per SAMPLE in highest_hits.txt file
def write_fhighest_hits(path_to_fhighest_hits, t_name_highest, q_name_highest, id_pct_highest, max_score):
    f = open(path_to_fhighest_hits, "a+")
    f.write(t_name_highest + "\t" + q_name_highest + "\t" + str(id_pct_highest) + "\t" + str(max_score) + "\n")
    f.close()


def check_cutoff(path_to_fhighest_hits, search_contig, path_to_fhighest_hits_filtered, path_to_fbelow_cutoff):
    with open(path_to_fhighest_hits, 'r') as myfile:
        for line in myfile:
            if line.startswith(search_contig):
                t_name_highest, q_name_highest, id_pct_highest, max_score = line.split("\t")
                if not string_in_file(path_to_fhighest_hits_filtered, search_contig):
                    if float(id_pct_highest) >= ID_PCT_CUTOFF and float(max_score) > SCORE_CUTOFF:
                        write_line(path_to_fhighest_hits_filtered, line)
                    else:
                        write_line(path_to_fbelow_cutoff, line)


def write_line(path, line):
    f = open(path, "a+")
    f.write(line)
    f.close()


# check match exon pairs and if match: creates new exon_name.fasta file for MAFFT input and append original exon
# sequence and contig if file does not exist yet. If already exist: only appends contigs.
# Writes all exons with multiple contigs in multiple_contigs.txt.
# and write pair in sequenced_exons.txt. If not a match: write pair in no_match_pairs.txt
def create_exon_ffasta(search_contig, path_to_contig_exon_match, path_to_fhighest_hits_filtered,
                       path_to_mapped_exons_species_dir, path_to_fseq_exons, path_to_fexons_seq,
                       sorted_list_consensus_dir, path_to_consensus_dir, path_to_fmultiple_contigs,
                       path_to_fno_match):
    contig_match = ""
    exon_match = ""
    fmatch = open(path_to_contig_exon_match, "rt")
    for line in fmatch:
        if line.startswith(search_contig):
            contig_match, exon_match = line.split("\t")

    contig_hit = ""
    exon_hit = ""
    fhits = open(path_to_fhighest_hits_filtered, "rt")
    for line in fhits:
        if line.startswith(search_contig):
            contig_hit, exon_hit, id_perc, psl_score = line.split("\t")

    exon_hit = exon_hit.strip()
    exon_match = exon_match.strip()
    if exon_hit == exon_match:
        path_to_fexon = path_to_mapped_exons_species_dir + exon_hit + '.fasta'
        create_ftxt(path_to_fexon)
        append_exon_seq(path_to_fexons_seq, exon_hit, path_to_fexon)
        append_contig_seq(sorted_list_consensus_dir, contig_hit, path_to_consensus_dir, path_to_fexon)
        write_fhits(path_to_fseq_exons, contig_hit, exon_hit)
        if multiple_contigs(path_to_fexon):
            write_fmultiple_hits(path_to_fmultiple_contigs, exon_hit)
    else:
        write_fhits(path_to_fno_match, contig_hit, exon_hit)


def write_fhits(path, contig_hit, exon_name):
    if not string_in_file(path_to_fseq_exons, contig_hit):
        f = open(path, "a+")
        f.write(contig_hit + "\t" + exon_name + "\n")
        f.close()


def append_exon_seq(path_to_fexons_seq, hit_exon, path_to_fexon):
    if not string_in_file(path_to_fexon, hit_exon):
        for record in SeqIO.parse(path_to_fexons_seq, "fasta"):
            if hit_exon == record.id:
                fexon = open(path_to_fexon, "a+")
                fexon.write(">" + record.id + "\n")
                fexon.write(str(record.seq) + "\n")


def append_contig_seq(sorted_list_consensus_dir, hit_contig, path_to_consensus_dir, path_to_fexon):
    if not string_in_file(path_to_fexon, hit_contig):
        for name_file in sorted_list_consensus_dir:
            name, fasta = name_file.split('.')
            if hit_contig == name:
                path_to_fseq = path_to_consensus_dir + name + '.fasta'
                exon_file = open(path_to_fexon, "a+")
                for record in SeqIO.parse(path_to_fseq, "fasta"):
                    exon_file.write(">" + record.id + "\n")
                    exon_file.write(str(record.seq) + "\n")
                exon_file.close()


def multiple_contigs(path_to_fexon):
    f = open(path_to_fexon, "rt")
    ncontigs = 0
    for line in f:
        if line.startswith(">Contig"):
            ncontigs += 1
    if ncontigs > 1:
        return True


def write_fmultiple_hits(path, exon_hit):
    if not string_in_file(path, exon_hit):
        f = open(path, "a+")
        f.write(exon_hit + "\n")
        f.close()


# Code starts here:
# STEP 1: Checks all fragment overlaps between contigs and exons
# creates separate dictionary for mapped contigs and target exons with their start and end position
results_dir, A_dir, sample, contig = PSL_FILE.split("/")
contig, at = contig.split("_")
contig, contig_nr = contig.split("contig")

path_to_exons_start_end = "./data/exons/AT_exon_enum.txt"
dictionary_exons = create_dictionary_start_end(path_to_exons_start_end)
dictionary_exons_lsorted = natural_sort(dictionary_exons)

path_to_fmapped_contigs = "results/A04_mapped_contigs/" + sample + "/mapped_contigs.txt"
contig_start, contig_end = get_start_end_contig(path_to_fmapped_contigs, contig_nr)

path_to_blat_stats = "./results/A06_identified_contigs_blat/" + sample + "/stats/"
create_dir(path_to_blat_stats)
path_to_contig_exon_match = path_to_blat_stats + "contig_exon_match_list.txt"
create_ftxt(path_to_contig_exon_match)

search_contig = "Contig" + contig_nr + "\t"
if not string_in_file(path_to_contig_exon_match, search_contig):
    check_overlap(contig_start, contig_end, dictionary_exons_lsorted, path_to_contig_exon_match, contig_nr)

# STEP 2: Reads all PSL files one by one and parse info in highest_hits.txt for every contig
path_to_fhighest_hits = path_to_blat_stats + "highest_hits.txt"
create_fhits(path_to_fhighest_hits)
if not string_in_file(path_to_fhighest_hits, search_contig):
    read_psl(PSL_FILE)

# STEP 3: Filter on cutoffs for score and ID%
# if > cutoff: paste information in highest_hits_filtered.txt; if < cutoff: copy in below_cutoff_pairs.txt
path_to_fhighest_hits_filtered = path_to_blat_stats + "highest_hits_filtered.txt"
path_to_fbelow_cutoff = path_to_blat_stats + "below_cutoff_pairs.txt"
create_fhits(path_to_fhighest_hits_filtered)
create_ftxt(path_to_fbelow_cutoff)
check_cutoff(path_to_fhighest_hits, search_contig, path_to_fhighest_hits_filtered, path_to_fbelow_cutoff)

# STEP 4: Check if highest_hits_filtered.txt matches are also overlapped fragments in contig_exon_match_list.txt
# If yes: creates new .fasta file for MAFFT input
path_to_mapped_exons = './results/A07_mapped_exons/'
path_to_mapped_exons_species_dir = path_to_mapped_exons + sample + "/"
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

path_to_fexons_seq = "./data/exons/exons_AT.fasta"
path_to_consensus_dir = "./results/A05_consensus_contigs/" + sample + "/"
list_consensus_dir = os.listdir(path_to_consensus_dir)
sorted_list_consensus_dir = natural_sort(list_consensus_dir)

create_exon_ffasta(search_contig, path_to_contig_exon_match, path_to_fhighest_hits_filtered,
                   path_to_mapped_exons_species_dir, path_to_fseq_exons, path_to_fexons_seq,
                   sorted_list_consensus_dir, path_to_consensus_dir, path_to_fmultiple_contigs,
                   path_to_fno_match)


# STEP 5: Creates configuration files for MAFFT
list_exons = os.listdir(path_to_mapped_exons_species_dir)
sorted_list_exons = natural_sort(list_exons)

path_to_exon_env = "./envs/config_exons.yaml"
if not sample_exist(path_to_exon_env):
    append_YAML(path_to_exon_env, sorted_list_exons)



'''
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


def sample_exist(path):
    fYAML = open(path, "rt")
    for line in fYAML:
        if SAMPLE in line:
            return True


def append_YAML(path, sorted_list_contigs):
    fYAML = open(path, "a+")
    fYAML.write(SAMPLE + ":\n")
    pattern = "*.fasta"
    for file in sorted_list_contigs:
        if fnmatch.fnmatch(file, pattern):
            exon_ffasta = file
            exon_name, ffasta = exon_ffasta.split(".fasta")
            fYAML.write("    - " + str(exon_name) + "\n")
    fYAML.close()


path_to_mapped_contigs_dir = "results/A04_mapped_contigs/"

# STEP 1: Checks all fragment overlaps between contigs and exons
# creates separate dictionary for mapped contigs and target exons with their start and end position
path_to_assembled_species_dir = path_to_mapped_contigs_dir + SAMPLE + "/"
path_to_fmapped_contigs = path_to_assembled_species_dir + "mapped_contigs.txt"
path_to_exons_enum = "./data/exons/AT_exon_enum.txt"
dictionary_contigs = create_dictionary_start_end(path_to_fmapped_contigs)
dictionary_exons = create_dictionary_start_end(path_to_exons_enum)
dictionary_contigs_lsorted = natural_sort(dictionary_contigs)
dictionary_exons_lsorted = natural_sort(dictionary_exons)

# add matching pairs in contig_exon_match_list.txt if mapped contig fragments matches/overlaps
# with the exon fragments
path_to_blat_species_dir = "./results/A06_identified_contigs_blat/" + SAMPLE + "/"
path_to_blat_stats = path_to_blat_species_dir + "stats/"
create_dir(path_to_blat_stats)

path_to_contig_exon_match = path_to_blat_stats + "contig_exon_match_list.txt"
create_ftxt(path_to_contig_exon_match)
check_overlap(dictionary_contigs_lsorted, dictionary_exons_lsorted, path_to_contig_exon_match)

# STEP 2: Reads all PSL files one by one and parse info in highest_hits.txt for every contig
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

# STEP 3: Filter on cutoffs for score and ID%
# if > cutoff: paste information in highest_hits_filtered.txt; if < cutoff: copy in below_cutoff_pairs.txt
path_to_fhighest_hits_filtered = path_to_blat_stats + "highest_hits_filtered.txt"
path_to_fbelow_cutoff = path_to_blat_stats + "below_cutoff_pairs.txt"
create_fhits(path_to_fhighest_hits_filtered)
create_ftxt(path_to_fbelow_cutoff)
check_cutoff(path_to_fhighest_hits, path_to_fhighest_hits_filtered, path_to_fbelow_cutoff)


# STEP 4: Check if highest_hits_filtered.txt matches are also overlapped fragments in contig_exon_match_list.txt
# If yes: creates new .fasta file for MAFFT input
# create dictionaries for highest_hits_filtered.txt and contig_exon_match_list.txt
dictionary_hits = create_dictionary_highest_hits_filtered(path_to_fhighest_hits_filtered)
dictionary_match = create_dictionary_match_pairs(path_to_contig_exon_match)
dictionary_hits_lsorted = natural_sort(dictionary_hits)
dictionary_match_lsorted = natural_sort(dictionary_match)

# create new mapped_exons dir for input
path_to_mapped_exons = './results/A07_mapped_exons/'
path_to_mapped_exons_species_dir = path_to_mapped_exons + SAMPLE + "/"
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
path_to_consensus_dir = "./results/A05_consensus_contigs/" + SAMPLE + "/"
list_consensus_dir = os.listdir(path_to_consensus_dir)
sorted_list_consensus_dir = natural_sort(list_consensus_dir)

create_exon_ffasta(dictionary_hits_lsorted, dictionary_match_lsorted, dictionary_hits, dictionary_match,
                   path_to_mapped_exons_species_dir, path_to_fseq_exons, path_to_fexons_seq,
                   sorted_list_consensus_dir, path_to_consensus_dir, path_to_fmultiple_contigs, path_to_fno_match)

# STEP 5: Creates configuration files for MAFFT
path_to_exons_dir = './results/A07_mapped_exons/' + SAMPLE + '/'
list_exons = os.listdir(path_to_exons_dir)
sorted_list_exons = natural_sort(list_exons)

path_to_exon_env = "./envs/config_exons.yaml"
if not sample_exist(path_to_exon_env):
    append_YAML(path_to_exon_env, sorted_list_exons)

'''
