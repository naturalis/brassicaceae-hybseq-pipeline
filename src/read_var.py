#!/usr/bin/python
# This script is to read lines in var file after the VARscan step
# It returns the sequence, number of SNPs and avcov
# Made by: Elfy Ly
# Date: 28 May 2020

import sys
import os

#SPECIES = "SRR8528336"
CONTIGNR = 8

path_to_assembled_exons_dir = "./results/assembled_exons/"

for root_dir, sub_dir, files in os.walk(path_to_assembled_exons_dir):
    for species in sub_dir:
        print(species)

fvar_name = "Contig" + str(CONTIGNR) + "_AT_sort.var"
path_to_fvar = path_to_assembled_exons_dir + species + "/var/" + fvar_name

seq = ""
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

print("avcov: " + str(avcov))
print("seq: " + seq)
print("nsnps: " + str(nsnps))
seq_length = len(seq)
print("seq length: " + str(seq_length))


# Create consensus directory if it doesn't exist
path_to_consensus_dir = "./results/consensus/"
path_to_consensus_species_dir = path_to_consensus_dir + species + "/"

if not os.path.exists(path_to_consensus_dir):
    os.mkdir(path_to_consensus_dir)
    print("Directory " + path_to_consensus_dir + " Created ")
else:
    print("Directory " + path_to_consensus_dir + " already exists")

if not os.path.exists(path_to_consensus_species_dir):
    os.mkdir(path_to_consensus_species_dir)
    print("Directory " + path_to_consensus_species_dir + " Created ")
else:
    print("Directory " + path_to_consensus_species_dir + " already exists")


ncontigs = 0
if seq_length > 0:
    f = open(path_to_consensus_species_dir + "Contig" + str(CONTIGNR) + ".txt", "w+")
    print("New text file created: " + "Contig" + str(CONTIGNR) + ".txt")
    f.write(">Contig" + str(CONTIGNR) + "\n" + seq + "\n")
    f.close()
    ncontigs += 1
else:
    print("Fasta file is not printed: " + "Contig " + str(CONTIGNR))


