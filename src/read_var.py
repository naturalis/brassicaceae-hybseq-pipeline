#!/usr/bin/python
# This script is to read lines in var file after the VARscan step and to make a consensus sequence
# It stores the sequence, number of SNPs and average coverage for every .var file
# It is used by giving the script the .var file as input argument, for example:
# $ python read_var.py results/A04_mapped_contigs/{sample}/var/Contig{nr}_AT_sort.var
# Made by: Elfy Ly
# Date: 28 May 2020

import sys
import os

INPUT_FILE = sys.argv[1]


def create_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        print("Directory " + path + " Created ")
    else:
        print("Directory " + path + " already exists")


# Creates new file for whole consensus.txt
def create_fconsensus(path):
    fconsensus = path + "consensus.txt"
    if not os.path.exists(fconsensus):
        f = open(fconsensus, "w+")
        print("New text file created :" + fconsensus)
        f.close()
    else:
        print(fconsensus + " already exists")


def create_file_empty(path_to_consensus_species_dir):
    fempty = path_to_consensus_species_dir + "empty_contigs.txt"
    if not os.path.exists(fempty):
        f = open(fempty, "w+")
        f.close()
        print("New text file created: empty_contigs.txt")
    else:
        print("Directory " + fempty + " already exists")


# parse var file information for every line and calculates and stores seq, nSNPs and avcov
def calculate_var_stats(seq):
    nsnps = 0
    #avcov = 0
    #nnon_snp = 0
    with open(INPUT_FILE, 'rt') as myfile:
        for myline in myfile:
            if not myline.startswith('Chrom'):
                chrom, position, ref, cons, reads1, reads2, varfreq, strands1, strands2, qual1, qual2, pvalue, \
                mapqual1, mapqual2, reads1plus, reads1minus, reads2plus, reads2minus, varallele = myline.split("\t")

                if cons != 'N':
                    seq += cons
                    #avcov += int(reads2)
                    #nnon_snp += 1
                elif cons == 'N':
                    nsnps += 1
    return seq


# Creates new txt file for consensus sequence
def make_consensus(seq_length, path, contig_nr, seq):
    if seq_length > 0:
        append_fconsensus(path, contig_nr, seq)
        create_fcontig_consensus(path, contig_nr, seq)
    else:
        append_empty_contigs(path_to_consensus_species_dir, contig_nr)
        create_fcontig_consensus(path, contig_nr, seq)
        print("Fasta file is empty: " + "Contig" + contig_nr + ".txt")


def append_fconsensus(path, contig_nr, seq):
    path_to_fconsensus = path + "consensus.txt"
    if os.stat(path_to_fconsensus).st_size == 0:
        f = open(path_to_fconsensus, "a+")
        f.write(">Contig" + contig_nr + "\n" + seq + "\n")
        f.close()
    else:
        contig_check = ">Contig" + contig_nr
        if not string_in_file(path_to_fconsensus, contig_check):
            f = open(path_to_fconsensus, "a+")
            f.write(">Contig" + contig_nr + "\n" + seq + "\n")
            f.close()


def string_in_file(file_name, string_to_search):
    """ Check if any line in the file contains given string """
    with open(file_name, 'r') as read_obj:
        for line in read_obj:
            if string_to_search in line:
                return True
    return False


def create_fcontig_consensus(path, contig_nr, seq):
    f = open(path + "Contig" + contig_nr + ".fasta", "w+")
    f.write(">Contig" + contig_nr + "\n" + seq + "\n")
    f.close()


def append_empty_contigs(path_to_consensus_species_dir, contig_nr):
    f = open(path_to_consensus_species_dir + "empty_contigs.txt", "a+")
    f.write(">Contig" + contig_nr + "\n")
    f.close()


# Code starts here
results_dir, A_dir, sample, var_dir, fvar = INPUT_FILE.split("/")
contig, at, sort = fvar.split("_")
contig, contig_nr = contig.split("Contig")

path_to_consensus_dir = "./results/A05_consensus_contigs/"
create_dir(path_to_consensus_dir)
path_to_consensus_species_dir = path_to_consensus_dir + sample + "/"
create_dir(path_to_consensus_species_dir)
create_fconsensus(path_to_consensus_species_dir)
create_file_empty(path_to_consensus_species_dir)

seq = ""
seq = calculate_var_stats(seq)
seq_length = len(seq)
make_consensus(seq_length, path_to_consensus_species_dir, contig_nr, seq)

