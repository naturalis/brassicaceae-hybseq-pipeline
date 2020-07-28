# trim_prot_alignments.py
# This script is to trim the protein aligned exons for all samples based on the ref exon in ORF.
# The input will be the EXON_NAME_AA.fasta and the output will be a EXON_NAME_AA_trimmed.fasta file with their
# trimmed protein alignments
# This script can be executed by running $ python trim_prot_alignments.py
# Made by: Elfy Ly
# Date: 2 July 2020

import os
import re
import fnmatch
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

ACGT_LENGTH_PERC_THRESH = 0.35
NX_PERC_TRESH = 0.20
DIVERSE_PERC_THRESH = 0.30


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
    ffasta = open(path, "w+")
    ffasta.close()


def replace_to_n(record):
    nt_dict[exon_name][record.id] = {}
    position = 0
    new_seq = ""
    for base in record.seq:
        position += 1
        nt_dict[exon_name][record.id][position] = base
        # change '!' fragment shifts and '-' gaps to 'N' in new sequence
        if base == "-" or base == "!":
            new_seq += "N"
        else:
            new_seq += base
    return new_seq


def count_bases(new_seq):
    nbase_ACGT = 0
    for base in new_seq:
        if base != "N":
            nbase_ACGT += 1
    return nbase_ACGT


# create empty dictionary to count x's per position
def create_empty_nested_dict(new_dict, old_dict):
    for exon in old_dict:
        for sample in old_dict[exon]:
            new_dict[exon] = {}
            for position in old_dict[exon][sample]:
                new_dict[exon][position] = 0
    return new_dict


# fill in the dictionary the number of x's per position
def count_X_dict(new_dict, old_dict):
    for exon in old_dict:
        for sample in old_dict[exon]:
            for position in old_dict[exon][sample]:
                if old_dict[exon][sample][position] == "X":
                    new_dict[exon][position] += 1
    return new_dict


def create_empty_trimmed_dict(aa_dict):
    trimmed_dict = {}
    for exon in aa_dict:
        trimmed_dict[exon] = {}
        for sample in aa_dict[exon]:
            trimmed_dict[exon][sample] = {}
    return trimmed_dict


def calculate_nx_thresh(aa_dict, exon):
    nseq = 0
    for sample in aa_dict[exon]:
        nseq += 1
    nx_thresh = NX_PERC_TRESH * nseq
    return nx_thresh


# initialize the nucleotides positions (regions) for each codon
# 1 aa position consists of 3 nt positions (1 AA - 1 NT, 2 NT, 3 NT)(2e AA - 4 NT, 5 NT, 6 NT) etc.
def initialize_nt_positions(aa_position):
    region = []
    nt1 = ((aa_position-1) * 3)+1   #-1 to intialize the first nucleotides
    nt2 = nt1 + 1
    nt3 = nt2 + 1
    region.append(nt1)
    region.append(nt2)
    region.append(nt3)
    return region


# counts for each exon per position if internal stop codons are present
def count_stop_dict(new_dict, trimmed_aa_x_dict, aa_dict):
    for exon in trimmed_aa_x_dict:
        for sample in trimmed_aa_x_dict[exon]:
            for position in trimmed_aa_x_dict[exon][sample]:
                last_position = len(aa_dict[exon][sample])
                if position != last_position and trimmed_aa_x_dict[exon][sample][position] == "*":
                    new_dict[exon][position] += 1
    return new_dict


def append_trimmed_dict(stop_dict, aa_dict, trimmed_aa_x_dict, trimmed_aa_x_stop_dict,
                    trimmed_nt_x_dict, trimmed_nt_x_stop_dict):
    for exon in stop_dict:
        for sample in aa_dict[exon]:
            for aa_position in stop_dict[exon]:
                if stop_dict[exon][aa_position] == 0:
                    # appends aa to the trimmed aa dictionary if there are no internal stop codons
                    trimmed_aa_x_stop_dict[exon][sample][aa_position] = trimmed_aa_x_dict[exon][sample][aa_position]
                    region = initialize_nt_positions(aa_position)
                    # print(str(aa_position) + " belongs to " + str(region))
                    for nt_position in region:
                        trimmed_nt_x_stop_dict[exon][sample][nt_position] = trimmed_nt_x_dict[exon][sample][nt_position]


def count_nstops(stop_dict):
    total_nstops = 0
    for exon in stop_dict:
        for position in stop_dict[exon]:
            if stop_dict[exon][position] > 0:
                total_nstops += 1
    return total_nstops


def create_empty_diversity_dict(count_aa_diversity, aa_dict, trimmed_aa_x_stop_dict):
    for exon in aa_dict:
        count_aa_diversity[exon] = {}
        for sample in trimmed_aa_x_stop_dict[exon]:
            for aa_position in trimmed_aa_x_stop_dict[exon][sample]:
                count_aa_diversity[exon][aa_position] = {}
    return count_aa_diversity


def count_diversity_dict(count_aa_diversity, trimmed_aa_x_stop_dict):
    for exon in count_aa_diversity:
        for sample in trimmed_aa_x_stop_dict[exon]:
            for aa_position in trimmed_aa_x_stop_dict[exon][sample]:
                aa = trimmed_aa_x_stop_dict[exon][sample][aa_position]
                appended_aa = count_aa_diversity[exon][aa_position].keys()
                if aa in appended_aa:
                    count_aa_diversity[exon][aa_position][aa] += 1
                else:
                    count_aa_diversity[exon][aa_position][aa] = 1
    return count_aa_diversity


def calculate_div_thresh(aa_dict, exon) :
    nseq = 0
    for sample in aa_dict[exon]:
        nseq += 1
    diverse_thresh = nseq * DIVERSE_PERC_THRESH
    return diverse_thresh


def create_empty_final_dict(count_aa_diversity, aa_dict):
    trimmed_dict = {}
    for exon in count_aa_diversity:
        trimmed_dict[exon] = {}
        for sample in aa_dict[exon]:
                trimmed_dict[exon][sample] = ""
    return trimmed_dict


def append_final_trimmed_dict(count_aa_diversity, aa_dict, trimmed_aa_x_stop_diverse_dict, trimmed_aa_x_stop_dict,
                              trimmed_nt_x_stop_diverse_dict, trimmed_nt_x_stop_dict):
    for exon in count_aa_diversity:
        diverse_thresh = calculate_div_thresh(aa_dict, exon)
        for sample in aa_dict[exon]:
            for aa_position in count_aa_diversity[exon]:
                n_most_prevelant = 0
                for aa in count_aa_diversity[exon][aa_position]:
                    if count_aa_diversity[exon][aa_position][aa] > n_most_prevelant:
                        n_most_prevelant = count_aa_diversity[exon][aa_position][aa]
                if n_most_prevelant > diverse_thresh:
                    trimmed_aa_x_stop_diverse_dict[exon][sample] += trimmed_aa_x_stop_dict[exon][sample][aa_position]
                    region = initialize_nt_positions(aa_position)
                    for nt_position in region:
                        trimmed_nt_x_stop_diverse_dict[exon][sample] += trimmed_nt_x_stop_dict[exon][sample][
                            nt_position]


def write_ffasta(final_trimmed_dict, path_to_trimmed_prot_dir, extension):
    for exon in final_trimmed_dict:
        path_to_nt_ffasta = path_to_trimmed_prot_dir + exon + extension
        nt_ffasta = open(path_to_nt_ffasta, "w+")
        for sample in final_trimmed_dict[exon]:
            nt_ffasta.write(">" + sample + "\n")
            nt_seq = final_trimmed_dict[exon][sample]
            nt_ffasta.write(nt_seq + "\n")
        nt_ffasta.close()


# Code starts here
path_to_macse_dir = "./results/A13_prot_alignments/"
list_macse_dir = os.listdir(path_to_macse_dir)
sorted_list_macse_dir = natural_sort(list_macse_dir)

path_to_trimmed_prot_dir = "./results/A14_trimmed_prot/"
create_dir(path_to_trimmed_prot_dir)

pattern = "*_NT.fasta"

max_exons = 0
aa_dict = {}
nt_dict = {}
for file in sorted_list_macse_dir:
    if fnmatch.fnmatch(file, pattern):
        fNT_alignments = file
        name, fasta = fNT_alignments.split(".fasta")
        exon_name, nt = name.split("_")
        max_exons += 1

        path_to_fNT = path_to_macse_dir + fNT_alignments
        aa_dict[exon_name] = {}
        nt_dict[exon_name] = {}
        for record in SeqIO.parse(path_to_fNT, "fasta"):
            # create new sequence to translate gaps and frameshifts to 'N'
            new_seq = replace_to_n(record)
            '''Step 1: Deletes the sample if sequence base length too short (if too many gaps or fragment shifts)
            Final sequences shorter than 35% of unambiguous nucleotide positions based on the reference exon length 
            were removed.'''
            nbase = len(record.seq)
            min_nbases = nbase * ACGT_LENGTH_PERC_THRESH
            nbase_ACGT = count_bases(new_seq)
            if nbase_ACGT > min_nbases:
                # translate dna sequence to protein sequence
                nt_seq = Seq(new_seq, generic_dna)
                protein_seq = nt_seq.translate()

                # assign amino acid per position in aa_dict
                position = 0
                aa_dict[exon_name][record.id] = {}
                for aa in protein_seq:
                    position += 1
                    aa_dict[exon_name][record.id][position] = aa
                protein_seq_length = len(protein_seq)
            else:
                print(record.id + " has too few nucleotide bases: " + str(nbase_ACGT) +
                      ". It's below min_bases: " + str(min_nbases)) + ". This sample has been deleted."

# evaluate protein alignment
'''Step 2: Calculates per exon for every position the number of X's. Positions with > 20% ambiguous amino acids 
(threshold) resulting from unidentified nucleotides (Ns) and gaps (-) were removed'''
x_dict = {}
x_dict = create_empty_nested_dict(x_dict, aa_dict)
x_dict = count_X_dict(x_dict, aa_dict)

# create trimmed dictionary for AA and NT
trimmed_aa_x_dict = create_empty_trimmed_dict(aa_dict)
trimmed_nt_x_dict = create_empty_trimmed_dict(aa_dict)

total_nx = 0
for exon in x_dict:
    nx_thresh = calculate_nx_thresh(aa_dict, exon)
    for sample in aa_dict[exon]:
        for aa_position in x_dict[exon]:
            nx = x_dict[exon][aa_position]
            # checks if number of X's per position are below threshold, if yes: add base to the trimmed dictionary for
            # protein sequences and nucleotide sequences
            if nx <= nx_thresh:
                trimmed_aa_x_dict[exon][sample][aa_position] = aa_dict[exon][sample][aa_position]

                region = initialize_nt_positions(aa_position)
                # print(str(aa_position) + " belongs to " + str(region))
                for nt_position in region:
                    trimmed_nt_x_dict[exon][sample][nt_position] = nt_dict[exon][sample][nt_position]
            else:
                total_nx += 1
                # print("number of X's: " + str(nx) + " on position: " + str(aa_position) + " exceeds the threshold: " +
                #       str(nx_thresh) + " in exon: " + exon)
#
'''Step 3: Checks if internal stop codon (*) indicative of misalignment is present'''
# creates empty dictionary for stop codons
stop_dict = {}
stop_dict = create_empty_nested_dict(stop_dict, trimmed_aa_x_dict)
stop_dict = count_stop_dict(stop_dict, trimmed_aa_x_dict, aa_dict)

# create trimmed AA and NT dictionary for each exon per position without gaps, shifts and internal stop codons
trimmed_aa_x_stop_dict = create_empty_trimmed_dict(aa_dict)
trimmed_nt_x_stop_dict = create_empty_trimmed_dict(aa_dict)
append_trimmed_dict(stop_dict, aa_dict, trimmed_aa_x_dict, trimmed_aa_x_stop_dict,
                    trimmed_nt_x_dict, trimmed_nt_x_stop_dict)

# Counts number of present internal stop codons (*) indicative of misalignment
total_nstops = count_nstops(stop_dict)

'''Step 4: Deletes if a codon position was too diverse (most prevalent amino acid identical for < 30% of the taxa)'''
# create dictionary to count different aa in each position
count_aa_diversity = {}
count_aa_diversity = create_empty_diversity_dict(count_aa_diversity, aa_dict, trimmed_aa_x_stop_dict)
count_aa_diversity = count_diversity_dict(count_aa_diversity, trimmed_aa_x_stop_dict)

trimmed_nt_x_stop_diverse_dict = create_empty_final_dict(count_aa_diversity, aa_dict)
trimmed_aa_x_stop_diverse_dict = create_empty_final_dict(count_aa_diversity, aa_dict)
append_final_trimmed_dict(count_aa_diversity, aa_dict, trimmed_aa_x_stop_diverse_dict, trimmed_aa_x_stop_dict,
                          trimmed_nt_x_stop_diverse_dict, trimmed_nt_x_stop_dict)

'''Step 5: Write fasta files for trimmed AA and NT sequences in directory path_to_trimmed_prot_dir'''
nt = "_NT.fasta"
aa = "_AA.fasta"
write_ffasta(trimmed_nt_x_stop_diverse_dict, path_to_trimmed_prot_dir, nt)
write_ffasta(trimmed_aa_x_stop_diverse_dict, path_to_trimmed_prot_dir, aa)
