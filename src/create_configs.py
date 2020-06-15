# create_configs.py
# This script is to create configuration files for every sample the MAFFT assembly in Snakemake
# It is executed by: $ python create_configs.py [SAMPLE_NAME], for example: $ python create_configs.py SRR8528336
# Made by: Elfy Ly
# Date: 15 June 2020

import sys
import os

SAMPLE_NAME = sys.argv[1]


# Create directory of given path if it doesn't exist
def create_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        print("Directory ", path, " Created ")
    else:
        print("Directory ", path, " already exists")


def create_YAML(path):
    f = open(path, "w+")
    f.write("exons:\n")
    f.close()
    print(path + " is created")


# Code starts here
path_to_MAFFT_configs = "./envs/MAFFT/"
create_dir(path_to_MAFFT_configs)

path_to_fYAML = path_to_MAFFT_configs + SAMPLE_NAME + ".yaml"
create_YAML(path_to_fYAML)

path_to_mapped_exons_dir = "./results/mapped_exons/" + SAMPLE_NAME + "/"
path_to_fseq_exons = path_to_mapped_exons_dir + "stats/sequenced_exons.txt"
with open(path_to_fseq_exons, "rt") as fseq_exons:
    for line in fseq_exons:
        contig_name, exon_name = line.split("\t")
        fYAML = open(path_to_fYAML, "w+")
        fYAML.write("\t" + exon_name + ": " + path_to_mapped_exons_dir + exon_name + ".txt")
        fYAML.close()
fseq_exons.close()


