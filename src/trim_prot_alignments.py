# trim_alignments.py
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


def write_ffasta(path, record, protein_seq):
    ffasta = open(path, "a+")
    ffasta.write(">" + record.id + "\n")
    ffasta.write(protein_seq + "\n")
    ffasta.close()


# Code starts here
path_to_macse_dir = "./results/13_prot_alignments/"
list_macse_dir = os.listdir(path_to_macse_dir)
sorted_list_macse_dir = natural_sort(list_macse_dir)

path_to_trimmed_prot_dir = "./results/14_trimmed_prot/"
create_dir(path_to_trimmed_prot_dir)

pattern = "*_NT.fasta"
max_exons = 0
exon_dict = {}
nt_dict = {}
for file in sorted_list_macse_dir:
    if fnmatch.fnmatch(file, pattern):
        fNT_alignments = file
        name, fasta = fNT_alignments.split(".fasta")
        exon_name, nt = name.split("_")
        max_exons += 1

        path_to_fNT = path_to_macse_dir + fNT_alignments
        path_to_ftrimmed_prot = path_to_trimmed_prot_dir + fNT_alignments
        create_ffasta(path_to_ftrimmed_prot)

        exon_dict[exon_name] = {}
        nt_dict[exon_name] = {}
        for record in SeqIO.parse(path_to_fNT, "fasta"):
            # create new sequence to translate gaps and frameshifts to 'N'
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

            '''Step 1: Deletes the sample if sequence base length too short (if too many gaps or fragment shifts)
            Final sequences shorter than 35% of unambiguous nucleotide positions based on the reference exon length 
            were removed.'''            ## DIT MOET EIGENLIJK LAATSTE STAP ZIJN, WANT HIJ MOET FINAL SEQUENCE CHECK DOEN?
            # calculate minimum nucleotide length per sequence and delete if below threshold
            nbase = len(record.seq)         ### = REFERENCE EXON LENGTH?? (Is geloof ik gewoon de ORF lengte, dus lengte die elke NT seq heeft in het begin na macse)
            min_nbases = nbase * ACGT_LENGTH_PERC_THRESH
            nbase_ACGT = 0
            for base in new_seq:
                if base != "N":
                    nbase_ACGT += 1

            if nbase_ACGT <= min_nbases:
                print(record.id + " has too few nucleotide bases: " + str(nbase) +
                      ". It's below min_bases: " + str(min_nbases))
                ### MOET NOG GEDAAN WORDEN: verwijder record.id en seq van lijst

            # translate dna sequence to protein sequence
            my_seq = Seq(new_seq, generic_dna)
            protein_seq = my_seq.translate()
            exon_dict[exon_name][record.id] = {}

            position = 0
            for aa in protein_seq:
                position += 1
                exon_dict[exon_name][record.id][position] = aa

                # if aa == "*":
                #     print(fNT_alignments)
                #     print(record.id)
                #     print(protein_seq)

            protein_seq_length = len(protein_seq)

# test
print("originele NT seq:")
print(nt_dict['AT5G05680.1@6'])
print("\n")
print("changed fragments shifts and gaps to 'N' and translated in AA seq:")
print(exon_dict['AT5G05680.1@6'])
print("\n")


# evaluate protein alignment
'''Step 2: Calculates per exon for every position the number of X's. Positions with > 20% ambiguous amino acids 
(threshold) resulting from unidentified nucleotides (Ns) and gaps (-) were removed'''
# create empty dictionary to count x's per position
x_dict = {}
for exon in exon_dict:
    for sample in exon_dict[exon]:
        x_dict[exon] = {}
        for position in exon_dict[exon][sample]:
            x_dict[exon][position] = 0

# fill in the dictionary the number of x's per position
for exon in exon_dict:
    for sample in exon_dict[exon]:
        for position in exon_dict[exon][sample]:
            if exon_dict[exon][sample][position] == "X":
                x_dict[exon][position] += 1

#test
print("dictionary counting the X's of each position of all samples:")
print(x_dict['AT5G05680.1@6'])
print("\n")


# create trimmed dictionary for AA and NT
trimmed_aa_x_dict = {}     # AA
trimmed_nt_x_dict = {}  # NT
n_xs = 0
for exon in x_dict:
    trimmed_aa_x_dict[exon] = {}
    trimmed_nt_x_dict[exon] = {}

    # calculate nx threshold
    nseq = 0
    for sample in exon_dict[exon]:
        nseq += 1
    nx_thresh = NX_PERC_TRESH * nseq

    # add trimmed alignments to new dictionary (aa_x_dict for amino acids and nt_x_dict for nucleotides)
    for sample in exon_dict[exon]:
        trimmed_aa_x_dict[exon][sample] = {}
        trimmed_nt_x_dict[exon][sample] = {}
        for aa_position in x_dict[exon]:
            nx = x_dict[exon][aa_position]
            # checks per position if number of X's are below threshold, if yes: add base to the trimmed dictionary for
            # protein sequences and nucleotide sequences
            if nx <= nx_thresh:
                # write trimmed aa sequence
                trimmed_aa_x_dict[exon][sample][aa_position] = exon_dict[exon][sample][aa_position]

                # initialize the nucleotides positions (regions) for each codon
                region = []
                # 1 aa position consists of 3 nt positions (1 AA - 1 NT, 2 NT, 3 NT)(2e AA - 4 NT, 5 NT, 6 NT) etc.
                nt1 = ((aa_position-1) * 3)+1   #-1 to intialize the first nucleotides
                nt2 = nt1 + 1
                nt3 = nt2 + 1
                region.append(nt1)
                region.append(nt2)
                region.append(nt3)
                # print(str(aa_position) + " belongs to " + str(region))

                # write trimmed nt sequence
                for nt_position in region:
                    # print(nt_position)
                    trimmed_nt_x_dict[exon][sample][nt_position] = nt_dict[exon][sample][nt_position]
            else:
                n_xs += 1
                # print("number of X's: " + str(nx) + " on position: " + str(aa_position) + " exceeds the threshold: " +
                #       str(nx_thresh) + " in exon: " + exon)

# test
print("trimmed dictionary of aa on the positions where amount of x was higher than threshold:")
print(trimmed_aa_x_dict['AT5G05680.1@6'])
print("\n")

print("trimmed dictionary of nucleotide regions of the aa's where X was higher than threshold:")
print(trimmed_nt_x_dict['AT5G05680.1@6'])
print("\n")

# test x's
print("number of x's above threshold: " + str(n_xs))
# print(trimmed_aa_x_dict['AT5G05680.1@6']) # positie 944 verwijderd


'''Step 3: Checks if internal stop codon (*) indicative of misalignment is present'''
# creates empty dictionary for stop codons
stop_dict = {}
for exon in trimmed_aa_x_dict:
    for sample in trimmed_aa_x_dict[exon]:
        stop_dict[exon] = {}
        for position in trimmed_aa_x_dict[exon][sample]:
            stop_dict[exon][position] = 0

# counts for each exon per position if internal stop codons are present
for exon in trimmed_aa_x_dict:
    for sample in trimmed_aa_x_dict[exon]:
        for position in trimmed_aa_x_dict[exon][sample]:
            last_position = len(exon_dict[exon][sample])
            if position != last_position and trimmed_aa_x_dict[exon][sample][position] == "*":
                stop_dict[exon][position] += 1

# test
print("dictionary with counted stop codons for all samples:")
print(stop_dict['AT5G05680.1@6'])
print("\n")


# create trimmed AA and NT dictionary for each exon per position without gaps, shifts and internal stop codons
trimmed_aa_x_stop_dict = {}
trimmed_nt_x_stop_dict = {}
for exon in stop_dict:
    trimmed_aa_x_stop_dict[exon] = {}
    trimmed_nt_x_stop_dict[exon] = {}
    for sample in exon_dict[exon]:
        trimmed_aa_x_stop_dict[exon][sample] = {}
        trimmed_nt_x_stop_dict[exon][sample] = {}
        for aa_position in stop_dict[exon]:
            if stop_dict[exon][aa_position] == 0:
                # appends aa to the trimmed aa dictionary if there are no internal stop codons
                trimmed_aa_x_stop_dict[exon][sample][aa_position] = trimmed_aa_x_dict[exon][sample][aa_position]

                # initialize the nucleotides positions (regions) for each codon (ZELFDE ALS BOVEN FUNCTION)
                region = []
                nt1 = ((aa_position-1) * 3)+1
                nt2 = nt1 + 1
                nt3 = nt2 + 1
                region.append(nt1)
                region.append(nt2)
                region.append(nt3)
                # print(str(aa_position) + " belongs to " + str(region))
                for nt_position in region:
                    trimmed_nt_x_stop_dict[exon][sample][nt_position] = trimmed_nt_x_dict[exon][sample][nt_position]

# test
print("trimmed dictionary of aa on the positions where amount of x was higher than threshold AND stop codons unless"
      "positioned on the last position:")
print(trimmed_aa_x_stop_dict['AT5G05680.1@6'])
print("\n")

print("trimmed dictionary of nucleotide regions of the aa's where X was higher than threshold AND stop codons unless"
      "positioned on the last position:")
print(trimmed_nt_x_stop_dict['AT5G05680.1@6'])
print("\n")


# Counts number of present internal stop codons (*) indicative of misalignment
nstops = 0
for exon in stop_dict:
    for position in stop_dict[exon]:
        if stop_dict[exon][position] > 0:
            nstops += 1
            # print("exon: " + exon + " on position: " + str(position) + " has a stopcodon")
            # print(stop_dict[exon])
print("nr of stop codons: " + str(nstops))


'''Step 4: Deletes if a codon position was too diverse (most prevalent amino acid identical for < 30% of the taxa)'''
# nieuwe lege dict voor trimmed x stop en diverse = (trimmed_aa_x_stop_diverse_dict = {} en trimmed_nt_x_stop_diverse_dict = {})
# nieuwe dict maken voor het tellen van de soorten base in prot alignment
trimmed_aa_x_stop_diverse_dict = {}
trimmed_nt_x_stop_diverse_dict = {}
# Voor elk exon:
for exon in exon_dict:
    nseq = 0
    for sample in exon_dict[exon]:
        # Hoeveel sequenties zijn er
        nseq += 1
    # Bereken threshold (= sequenties * THRESHOLD)
    diverse_thresh = nseq * DIVERSE_PERC_THRESH
    for sample in trimmed_aa_x_stop_dict[exon]:
        # Voor elk positie:
        for aa_position in trimmed_aa_x_stop_dict[exon][sample]:
            # Welke base komt het meest voor (=prevelant_base) EN hoeveel zijn dat er (= n_most_prevelant)

            # Als n_most_prevelant > threshold:
                # toevoegen aan trimmed_aa_x_stop_diverse_dict
                # heel dat region geval
                # toevoegen aan trimmed_nt_x_stop_diverse_dict
            # else:
                # print("in exon [exon] op position [position] komt base [prevalant_base] zo vaak: n_most_prevelant voor
                # en dit is minde dan de threshold (threshold). Deze zijn kolommen zijn verwijderd")

# Deze moet dan pakken van trimmed_aa_x_stop_diverse_dict = {} en trimmed_nt_x_stop_diverse_dict
'''Step 5: Create final dictionary for trimmed AA and NT sequences'''
trimmed_dict = {}
for exon in trimmed_aa_x_stop_dict:
    trimmed_dict[exon] = {}
    for sample in trimmed_aa_x_stop_dict[exon]:
        trimmed_dict[exon][sample] = ""
        for position in trimmed_aa_x_stop_dict[exon][sample]:
            trimmed_dict[exon][sample] += trimmed_aa_x_stop_dict[exon][sample][position]

# test
print("trimmed aa sequences: ")
print(trimmed_dict['AT5G05680.1@6'])
print("\n")

trimmed_nt_dict = {}
for exon in trimmed_nt_x_stop_dict:
    trimmed_nt_dict[exon] = {}
    for sample in trimmed_nt_x_stop_dict[exon]:
        trimmed_nt_dict[exon][sample] = ""
        for position in trimmed_nt_x_stop_dict[exon][sample]:
            trimmed_nt_dict[exon][sample] += trimmed_nt_x_stop_dict[exon][sample][position]

# test
print("trimmed nt sequences: ")
print(trimmed_nt_dict['AT5G05680.1@6'])
print("\n")


'''Step 6: Write fasta files for trimmed AA and NT sequences in directory path_to_trimmed_prot_dir'''

                # write_ffasta(path_to_ftrimmed_prot, record, protein_seq)

# print(max_exons)




