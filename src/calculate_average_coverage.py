# This script is to calculate the average coverage for every sample
# Execution: python calculate_average_coverage.py [SAMPLENAME]
# Date: 14 August 2020
# Made by: Elfy Ly

import sys
import os
import fnmatch

SAMPLE = sys.argv[1]


path_to_mapped_contigs_dir = "./results/A04_mapped_contigs/"
path_to_mapped_contigs_var_dir = path_to_mapped_contigs_dir + SAMPLE + "/var/"

list_in_var_dir = os.listdir(path_to_mapped_contigs_var_dir)
pattern = "*.var"
ncontigs = 0
coverage_list = []
for file in list_in_var_dir:
    if fnmatch.fnmatch(file, pattern):
        ncontigs += 1
        var_file = file
        path_to_fvar = path_to_mapped_contigs_var_dir + var_file
        # print(var_file)
        nreads = 0
        nsnps = 0
        avcov = 0
        nnon_snp = 0
        with open(path_to_fvar, 'rt') as myfile:
            for myline in myfile:
                if not myline.startswith('Chrom'):
                    chrom, position, ref, cons, reads1, reads2, varfreq, strands1, strands2, qual1, qual2, pvalue, \
                    mapqual1, mapqual2, reads1plus, reads1minus, reads2plus, reads2minus, varallele = myline.split("\t")

                    if cons != 'N':
                        avcov += int(reads2)
                        nnon_snp += 1

                    nreads += 1

        if nnon_snp > 0:
            avcov /= nnon_snp

        # print("average coverage is: " + str(round(avcov)) + " and " + str(nreads) + " reads in file: " + var_file)
        coverage_list.append(round(avcov))

sum_coverage = 0
for contig in coverage_list:
    sum_coverage += contig
avcov = sum_coverage / ncontigs
print("The average coverage of " + SAMPLE + " is " + str(round(avcov)))

path_to_favcov = path_to_mapped_contigs_var_dir + SAMPLE + "_average_coverage.txt"
favcov = open(path_to_favcov, "w+")
favcov.write(str(round(avcov)))
favcov.close()
