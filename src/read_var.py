# This script is to read lines in var file after the VARscan step
# It returns the sequence, number of SNPs and avcov
# Made by: Elfy Ly
# Date: 28 May 2020

SPECIES = "SRR8528336"
CONTIGNR = 8

path_to_assembled_exons_dir = "./results/assembled_exons/"
path_to_fvar = path_to_assembled_exons_dir + SPECIES + "/var/" + "Contig" + str(CONTIGNR) + "_AT_sort.var"

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
                avcov += reads2

if nsnps > 0:
    avcov = avcov / nsnps

print("avcov: " + str(avcov))
print("seq: " + seq)
print("nsnps: " + str(nsnps))
