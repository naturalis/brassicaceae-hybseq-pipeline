# retrieve_exons_sequence_genomes.py
# This script is to retrieve exons from sequenced genomes
# To identify the contigs from the sequenced genomes, each contig has to be retrieved from the sequenced genome first.
# Then, for each sequence query of the sequenced genomes, the query can be BLAT against the database reference.
# In this case, the database reference will be Arabidopsis thaliana.
# Made by: Elfy Ly
# Date: 19 May 2020

import os
from Bio import SeqIO

path_to_sis_database = "/mnt/c/Users/elfyl/PycharmProjects/brassicaceae-hybseq-pipeline-offline/data/sequence_genomes/sis_ref.fna"
path_to_sis_exons_dir = "/mnt/c/Users/elfyl/PycharmProjects/brassicaceae-hybseq-pipeline-offline/data/sequence_genomes/sis_ref"

# Create target Directory if don't exist
if not os.path.exists(path_to_sis_exons_dir):
    os.mkdir(path_to_sis_exons_dir)
    print("Directory ", path_to_sis_exons_dir, " Created ")
else:
    print("Directory ", path_to_sis_exons_dir, " already exists")

# Create new files for every sequence query of the sis database
count_id = 0
for seq_record in SeqIO.parse(path_to_sis_database, "fasta"):
    f = open(path_to_sis_exons_dir + "/" + seq_record.id + ".txt", "w+")
    print("New text file created: " + seq_record.id + ".fa")
    seq_id = seq_record.id
    seq_seq = str(seq_record.seq)
    f.write(">" + seq_id + "\n" + seq_seq)
    f.close()

    count_id += 1

print("Number of sequence records: " + str(count_id))
