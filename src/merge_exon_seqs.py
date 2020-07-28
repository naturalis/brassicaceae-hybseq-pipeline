# merge_exon_seqs.py
# This script is to merge all consensus exons for all samples to one file and adds the original exon sequence
# in ORF for MAFFT.
# It creates a .fasta file for each exon containing al consensus sequences for all samples.
# This script can be executed by running $ python merge_exon_seqs.py
# Made by: Elfy Ly
# Date: 29 June 2020

import os
import re
from Bio import SeqIO


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


def create_fexon(path_to_orf_ref):
    for record in SeqIO.parse(path_to_orf_ref, "fasta"):
        at, exon_name = record.id.split("_")
        path_to_fexon = path_to_all_samples_exons_dir + exon_name + ".fasta"
        exon_file = open(path_to_fexon, "w+")
        exon_file.write(">" + exon_name + "\n")
        exon_file.write(str(record.seq) + "\n")
        exon_file.close()
        print("File: " + path_to_fexon + " has been created.")


def append_seq(path_to_fexon, path_to_fexon_consensus):
    file = open(path_to_fexon, "a+")
    for record in SeqIO.parse(path_to_fexon_consensus, "fasta"):
        file.write(">" + record.id + "\n")
        file.write(str(record.seq) + "\n")
        file.close()


def create_YAML(path):
    f = open(path, "w+")
    f.write("exons:\n")
    f.close()
    print(path + " is created")


# Code starts here
path_to_orf_ref = "./data/exons/ref-at_orf.fasta"

# list to loop over the samples
path_to_consensus_dir = "./results/A09_consensus_exons/"
list_consensus_dir = os.listdir(path_to_consensus_dir)
sorted_consensus_dir = natural_sort(list_consensus_dir)


# creates new dir for output and list to check if file already exists
path_to_all_samples_exons_dir = "./results/A10_all_samples_exons/"
create_dir(path_to_all_samples_exons_dir)
list_exons_dir = os.listdir(path_to_all_samples_exons_dir)
sorted_exons_dir = natural_sort(list_exons_dir)

create_fexon(path_to_orf_ref)

# loops over samples in consensus dir
for sample in sorted_consensus_dir:
    path_to_consensus_sample_dir = path_to_consensus_dir + sample + "/"
    list_consensus_sample_dir = os.listdir(path_to_consensus_sample_dir)
    sorted_consensus_sample_dir = natural_sort(list_consensus_sample_dir)

    for fexon in sorted_consensus_sample_dir:
        path_to_fexon_consensus = path_to_consensus_sample_dir + fexon
        path_to_fexon = path_to_all_samples_exons_dir + fexon

        append_seq(path_to_fexon, path_to_fexon_consensus)

# Creates configuration YAML file for MAFFT in envs exons dir
path_to_fYAML = "./envs/all_exons.yaml"
create_YAML(path_to_fYAML)

for fexon in sorted_exons_dir:
    path_to_fexon = path_to_all_samples_exons_dir + fexon
    exon_name, fasta = fexon.split(".fasta")
    exon_name = exon_name.strip()
    fYAML = open(path_to_fYAML, "a+")
    fYAML.write("    - " + exon_name + "\n")
    fYAML.close()
