# make_exon_consensus.py
# This script is to create a consensus sequence after doing MAFFT assembly
# The input will be the SAMPLE_NAME and the output will be a EXON_NAME.fasta file with the consensus sequence
# This script can be executed by running $ python make_exon_consensus.py [SAMPLE_NAME], for example:
# $ python make_exon_consensus.py SRR8528336
# Made by: Elfy Ly
# Date: 15 June 2020

import sys
import os
from Bio import SeqIO
import re

SAMPLE_NAME = sys.argv[1]


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


def create_empty_dict(path_to_fassembled_exon, ncontigs, cons_dict, seq_length):
    for record in SeqIO.parse(path_to_fassembled_exon, "fasta"):
        if "Contig" in record.id:
            ncontigs += 1
            seq_length = len(str(record.seq))
            for position in range(1, seq_length + 1):
                cons_dict[position] = {}
    return cons_dict, seq_length


def getKeys(cons_dict, position):
    return cons_dict[position].keys()


def append_dict(path_to_fassembled_exon, seq_length, cons_dict):
    for record in SeqIO.parse(path_to_fassembled_exon, "fasta"):
        if "Contig" in record.id:
            for position in range(1, seq_length + 1):
                base = record.seq[position - 1]
                bases = getKeys(cons_dict, position)
                if base not in bases:
                    cons_dict[position][base] = 1
                elif base in bases:
                    cons_dict[position][base] += 1
    return cons_dict


def make_consensus(seq_length, cons_dict):
    consensus_seq = ""
    for position in range(1, seq_length + 1):
        sorted_cons_dict = sorted(cons_dict[position].items(), key=lambda x: x[1])

        nbases = 0
        bases = getKeys(cons_dict, position)
        for base in bases:
            nbases += 1
        if nbases == 1:
            for base in bases:
                consensus_seq += str(base)
        elif nbases == 2:
            bases = getKeys(cons_dict, position)

            base1 = sorted_cons_dict[0][0]
            base2 = sorted_cons_dict[1][0]
            n_bases1 = sorted_cons_dict[0][1]
            n_bases2 = sorted_cons_dict[1][1]

            if '-' in bases:
                for base in bases:
                    if base != "-":
                        consensus_seq += str(base)
            elif n_bases1 == n_bases2:
                consensus_seq += "N"
            elif n_bases1 > n_bases2:
                consensus_seq += base1
            elif n_bases2 > n_bases1:
                consensus_seq += base2

        elif nbases > 2:
            consensus_seq += "N"

    return consensus_seq


def create_fconsensus_exon(path_to_fconsensus_exon, exon_name):
    fconsensus_exon = open(path_to_fconsensus_exon, "w+")
    fconsensus_exon.write(">" + exon_name + "\n")
    fconsensus_exon.write(consensus_seq + "\n")
    fconsensus_exon.close()


# Code starts here
path_to_consensus_exons_dir = "./results/consensus_exons/"
path_to_consensus_exons_species_dir = path_to_consensus_exons_dir + SAMPLE_NAME + "/"
create_dir(path_to_consensus_exons_dir)
create_dir(path_to_consensus_exons_species_dir)

path_to_assembled_exon_dir = "./results/assembled_exons/" + SAMPLE_NAME + "/"
assembled_exon_dir = os.listdir(path_to_assembled_exon_dir)
sorted_assembled_exon_dir = natural_sort(assembled_exon_dir)
consensus_dict = {}
for fassembled_exon in sorted_assembled_exon_dir:
    path_to_fassembled_exon = path_to_assembled_exon_dir + fassembled_exon
    ncontigs = 0
    cons_dict = {}
    seq_length = 0
    cons_dict, seq_length = create_empty_dict(path_to_fassembled_exon, ncontigs, cons_dict, seq_length)
    cons_dict = append_dict(path_to_fassembled_exon, seq_length, cons_dict)

    consensus_seq = make_consensus(seq_length, cons_dict)

    path_to_fconsensus_exon = path_to_consensus_exons_species_dir + fassembled_exon
    exon_name, fasta = fassembled_exon.split(".fasta")
    create_fconsensus_exon(path_to_fconsensus_exon, exon_name)
