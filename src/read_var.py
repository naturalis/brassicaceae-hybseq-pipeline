#!/usr/bin/python
# This script is to read lines in var file after the VARscan step and to make a consensus sequence
# It loops through every contig in every species_name file and stores the sequence, number of SNPs and avcov of
# every contig
# Made by: Elfy Ly
# Date: 28 May 2020

import os


def count_contigs(path_to_var_dir):
    n_contig_files = 0
    for contig_file in os.listdir(path_to_var_dir):
        n_contig_files += 1
    return n_contig_files


# parse var file information for every line and calculates and stores seq, nSNPs and avcov
def calculate_var_stats(path_to_fvar, seq):
    nsnps = 0
    avcov = 0
    with open(path_to_fvar, 'rt') as myfile:
        for myline in myfile:
            if not myline.startswith('Chrom'):
                chrom, position, ref, cons, reads1, reads2, varfreq, strands1, strands2, qual1, qual2, pvalue, \
                mapqual1, mapqual2, reads1plus, reads1minus, reads2plus, reads2minus, varallele = myline.split("\t")

                if cons != 'N':
                    seq += cons
                elif cons == 'N':
                    nsnps += 1
                    avcov += int(reads2)

        if nsnps > 0:
            avcov /= nsnps

        # print("avcov: " + str(avcov))
        # print("seq: " + seq)
        # print("nsnps: " + str(nsnps))
        return seq


def create_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        print("Directory " + path + " Created ")
    else:
        print("Directory " + path + " already exists")


# Creates new file for whole consensus.txt
def create_fconsensus(path):
    f = open(path + "consensus.txt", "w+")
    print("New text file created :" + path + "consensus.txt")
    f.close()


# Creates new txt file for consensus sequence
def make_consensus(seq_length, path, contig_number, seq, ncontigs_after_VARscan):
    if seq_length > 0:
        append_fconsensus(path, contig_number, seq)
        create_fcontig_consensus(path, contig_number, seq)
        ncontigs_after_VARscan += 1
    else:
        append_empty_contigs(path_to_consensus_species_dir, contig_number)
        create_fcontig_consensus(path, contig_number, seq)
        print("Fasta file is empty: " + "Contig" + str(contig_number) + ".txt")

    return ncontigs_after_VARscan


def append_fconsensus(path, contig_number, seq):
    f = open(path + "consensus.txt", "a+")
    f.write(">Contig" + str(contig_number) + "\n" + seq + "\n")
    f.close()


def create_fcontig_consensus(path_to_consensus_species_dir, contig_number, seq):
    f = open(path_to_consensus_species_dir + "Contig" + str(contig_number) + ".txt", "w+")
    # print("New text file created: " + "Contig" + str(contig_number) + ".txt")
    f.write(">Contig" + str(contig_number) + "\n" + seq + "\n")
    f.close()


def create_file_empty(path_to_consensus_species_dir):
    f = open(path_to_consensus_species_dir + "empty_contigs.txt", "w+")
    f.close()
    print("New text file created: empty_contigs.txt")


def append_empty_contigs(path_to_consensus_species_dir, contig_number):
    f = open(path_to_consensus_species_dir + "empty_contigs.txt", "a+")
    f.write(">Contig" + str(contig_number) + "\n")
    f.close()


# Code starts here
path_to_assembled_exons_dir = "./results/assembled_exons/"

dirs = os.listdir(path_to_assembled_exons_dir)
for species_name in dirs:
    path_to_var_dir = path_to_assembled_exons_dir + species_name + "/var/"
    n_contig_files = count_contigs(path_to_var_dir)

    # Create consensus directory en file for every species if it doesn't exist
    path_to_consensus_dir = "./results/consensus/"
    path_to_consensus_species_dir = path_to_consensus_dir + species_name + "/"

    create_dir(path_to_consensus_dir)
    create_dir(path_to_consensus_species_dir)
    create_fconsensus(path_to_consensus_species_dir)
    create_file_empty(path_to_consensus_species_dir)

    # Loops through all contig files for each species
    ncontigs_after_VARscan = 0
    for contig_number in range(1, n_contig_files + 1):
        fvar_name = "Contig" + str(contig_number) + "_AT_sort.var"
        path_to_fvar = path_to_var_dir + fvar_name

        seq = ""
        seq = calculate_var_stats(path_to_fvar, seq)
        seq_length = len(seq)
        # print("seq length: " + str(seq_length))

        ncontigs_after_VARscan = make_consensus(seq_length, path_to_consensus_species_dir, contig_number, seq,
                                                ncontigs_after_VARscan)

    # Prints number of contig files created after VARscan (some files are empty after VARscan)
    print("contigs after VARscan: " + str(ncontigs_after_VARscan))
