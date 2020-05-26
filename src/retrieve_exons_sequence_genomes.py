# retrieve_exons_sequence_genomes.py
# This script is to retrieve exons from sequenced genomes which are also present in the reference genome (A. thaliana).
# To identify the contigs from the sequenced genomes, each contig has to be retrieved from A. thaliana first.
# Then, for each sequence query of A. thaliana, the query can be BLAT against the database reference.
# In this case, the database reference will be S. irio and A. lyrata.
# Made by: Elfy Ly
# Date: 19 May 2020

import os
from Bio import SeqIO

path_to_at_exons_dir = "/mnt/c/Users/elfyl/PycharmProjects/brassicaceae-hybseq-pipeline-offline/results/exons"
path_to_at_dir = "/mnt/c/Users/elfyl/PycharmProjects/brassicaceae-hybseq-pipeline-offline/data/reference_genomes"
path_to_at_reference = path_to_at_dir + "/ref-at.fasta"


# Create exons_AT Directory if don't exist
if not os.path.exists(path_to_at_exons_dir):
    os.mkdir(path_to_at_exons_dir)
    print("Directory ", path_to_at_exons_dir, " Created ")
else:
    print("Directory ", path_to_at_exons_dir, " already exists")

# Create new files for every sequence query of the reference genome A. thaliana
count_id = 0
for seq_record in SeqIO.parse(path_to_at_reference, "fasta"):
    f = open(path_to_at_exons_dir + "/" + seq_record.id + ".txt", "w+")
    print("New text file created: " + seq_record.id + ".fa")
    seq_id = seq_record.id
    seq_seq = str(seq_record.seq)
    f.write(">" + seq_id + "\n" + seq_seq)
    f.close()

    count_id += 1

print("Number of sequence records: " + str(count_id))
