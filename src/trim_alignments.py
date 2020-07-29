# trim_alignments.py
# This script is to trim the aligned exons for all samples based on the ref exon in ORF.
# The input will be the EXON_NAME.fasta file and the output will be a EXON_NAME.fasta file with the trimmed alignments
# This script can be executed by running $ python trim_alignments.py EXON_NAME
# Made by: Elfy Ly
# Date: 30 June 2020

import sys
import re
import os
from Bio import SeqIO


EXON = sys.argv[1]


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


# Create directory of given path if it doesn't exist
def create_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        print("Directory ", path, " Created ")
    else:
        print("Directory ", path, " already exists")


def create_ffasta(path):
    ftrimmed_exons = open(path, "w+")
    ftrimmed_exons.close()
    print("File: " + path + " has been created.")


def find_ref_start_end(record, ref_start, ref_end):
    for base in record.seq:
        if base == "-":
            ref_start += 1
        else:
            break
    reversed_seq = record.seq[::-1]
    for base in reversed_seq:
        if base == "-":
            ref_end += 1
        else:
            break
    ref_end = len(record.seq) - ref_end
    return ref_start, ref_end


def get_trimmed_alignment(trimmed_alignment, record):
    nbase = 0
    for base in record.seq:
        nbase += 1
        if nbase > ref_start:
            trimmed_alignment += base
        if nbase == ref_end:
            break
    return trimmed_alignment


def append_trimmed_alignments(path_to_ftrimmed_exons, record, trimmed_alignment):
    ftrimmed_exons = open(path_to_ftrimmed_exons, "a+")
    ftrimmed_exons.write(">" + record.id + "\n")
    ftrimmed_exons.write(trimmed_alignment + "\n")
    ftrimmed_exons.close()


def create_YAML(path):
    f = open(path, "w+")
    f.write("exons_NT:\n")
    f.close()
    print(path + " is created")


# Code starts here
path_to_alignments_dir = "results/A11_aligned_exons_ORF/"
# list_alignments = os.listdir(path_to_alignments_dir)
# sorted_list_alignments = natural_sort(list_alignments)

path_to_trimmed_dir = "results/A12_trimmed_alignments/"
create_dir(path_to_trimmed_dir)

# for fexon in sorted_list_alignments:
fexon = EXON + ".fasta"
path_to_fexon = path_to_alignments_dir + fexon
path_to_ftrimmed_exon = path_to_trimmed_dir + fexon
create_ffasta(path_to_ftrimmed_exon)

ref_start = 0
ref_end = 0
for record in SeqIO.parse(path_to_fexon, "fasta"):
    # record.id starting with AT is the ref seq in ORF
    if record.id.startswith("AT"):
        ref_start, ref_end = find_ref_start_end(record, ref_start, ref_end)
    else:
        trimmed_alignment = ""
        trimmed_alignment = get_trimmed_alignment(trimmed_alignment, record)
        append_trimmed_alignments(path_to_ftrimmed_exon, record, trimmed_alignment)

# # create .yaml for macse
# path_to_yaml = "./envs/macse.yaml"
# create_YAML(path_to_yaml)
#
# list_trimmed_exon_ffasta = os.listdir(path_to_trimmed_dir)
# sorted_list_trimmed_ffasta = natural_sort(list_trimmed_exon_ffasta)
# for trimmed_ffasta in sorted_list_alignments:
#     exon_name, fasta = trimmed_ffasta.split(".fasta")
#     macse_NT_output = "results/A13_prot_alignments/" + exon_name + "_NT.fasta"
#
#     fYAML = open(path_to_yaml, "a+")
#     fYAML.write("    " + exon_name + ": " + macse_NT_output + "\n")
#     fYAML.close()
#
# fYAML = open(path_to_yaml, "a+")
# fYAML.write("exons_AA:\n")
# fYAML.close()
#
# for trimmed_ffasta in sorted_list_alignments:
#     exon_name, fasta = trimmed_ffasta.split(".fasta")
#     macse_AA_output = "results/A13_prot_alignments/" + exon_name + "_AA.fasta"
#
#     fYAML = open(path_to_yaml, "a+")
#     fYAML.write("    " + exon_name + ": " + macse_AA_output + "\n")
#     fYAML.close()
#
