# extract_contigs_YASRA.py
# This script is to extract contigs and to create SAM files with headers for each of them after running YASRA.
# YASRA creates a SAM file without header after mapping.
# Made by: Elfy Ly
# Date: 7 May 2020

import os

path_to_YASRA_fSAM = '/results/alignments/SRR8528336_reads.fq_ref-at.fasta_Mon-May-25-08:50:00-2020/YASRA_related_files/alignments_SRR8528336_reads.fq_ref-at.fasta.sam'
path_to_assembled_exons = './results/assembled_exons'
path_to_species = path_to_assembled_exons + '/SRR8628336'
path_to_txt = path_to_species + '/txt'

# Create assembled_exons directory if it doesn't exist
if not os.path.exists(path_to_assembled_exons):
    os.mkdir(path_to_assembled_exons)
    print("Directory ", path_to_assembled_exons, " Created ")
else:
    print("Directory ", path_to_assembled_exons, " already exists")

# Create species directory if it doesn't exist
if not os.path.exists(path_to_species):
    os.mkdir(path_to_species)
    print("Directory ", path_to_species, " Created ")
else:
    print("Directory ", path_to_species, " already exists")

# Create txt directory if it doesn't exist
if not os.path.exists(path_to_txt):
    os.mkdir(path_to_txt)
    print("Directory ", path_to_txt, " Created ")
else:
    print("Directory ", path_to_txt, " already exists")

ncontigs = 0
nexons = 0
contig_name_temporary = " "

with open(path_to_YASRA_fSAM, 'rt') as myfile:
    for myline in myfile:
        if myline.startswith('>'):
            ncontigs += 1
            qname, flag, contig, pos, mapq, cigar, rnext, pnext, tlen, read, phred  = myline.split("\t")

            contig_number, reference_genome, contig_start, contig_end = contig.split("_")
            contig_length = int(contig_end) - int(contig_start) + 1

            if contig != contig_name_temporary:
                f = open(path_to_assembled_exons + "/" + contig_number + "_" + reference_genome + ".txt", "w+")
                print("New text file created: " + contig_number + "_" + reference_genome)

                f.write("@HD\tVN:1.3\n@SQ\tSN:" + contig + "\tLN:" + str(contig_length) + "\n" + myline)
                f.close()

                nexons +=1
                contig_name_temporary = contig

            elif contig == contig_name_temporary:
                f = open(path_to_assembled_exons + "/" + contig_number + "_" + reference_genome + ".txt", "a+")
                f.write(myline)
                f.close()
#
    f = open(path_to_assembled_exons + "/number_of_contigs_and_exons.txt", "w+")
    f.write("Number of contigs: " + str(ncontigs) + "\nNumber of exons: " + str(nexons))




