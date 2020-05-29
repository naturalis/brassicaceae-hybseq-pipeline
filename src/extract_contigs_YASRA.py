# extract_contigs_YASRA.py
# This script is to extract contigs and to create SAM files with headers for each of them after running YASRA.
# YASRA creates a SAM file without header after mapping.
# Made by: Elfy Ly
# Date: 7 May 2020

import os


# Create directory of given path if it doesn't exist
def create_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        print("Directory ", path, " Created ")
    else:
        print("Directory ", path, " already exists")


def create_new_fSAM(path_to_txt, contig_number, reference_genome, contig, contig_length, myline):
    f = open(path_to_txt + "/" + contig_number + "_" + reference_genome + ".txt", "w+")
    print("New text file created: " + contig_number + "_" + reference_genome)

    f.write("@HD\tVN:1.3\n@SQ\tSN:" + contig + "\tLN:" + str(contig_length) + "\n" + myline)
    f.close()


def write_to_fSAM(path_to_txt, contig_number, reference_genome, myline):
    f = open(path_to_txt + "/" + contig_number + "_" + reference_genome + ".txt", "a+")
    f.write(myline)
    f.close()


def count_save_stats(path_to_txt, ncontigs, nexons):
    f = open(path_to_txt + "/number_of_contigs_and_exons.txt", "w+")
    f.write("Number of contigs: " + str(ncontigs) + "\nNumber of exons: " + str(nexons) + "\n")
    f.close()


# Creating paths for output directory
# Deze paths moeten nog automatisch laten loopen, maar daarvoor moet eerst deze output directory naam gewijzigd worden
path_to_YASRA_fSAM = './results/alignments/SRR8528336_reads.fq_ref-at.fasta_Mon-May-25-08:50:00-2020/' \
                     'YASRA_related_files/alignments_SRR8528336_reads.fq_ref-at.fasta.sam'
path_to_assembled_exons = './results/assembled_exons'
path_to_species = path_to_assembled_exons + '/SRR8528336'
path_to_txt = path_to_species + '/txt'

create_dir(path_to_assembled_exons)
create_dir(path_to_species)
create_dir(path_to_txt)


ncontigs = 0
nexons = 0
contig_name_temporary = " "

with open(path_to_YASRA_fSAM, 'rt') as myfile:
    for myline in myfile:
        if myline.startswith('>'):
            ncontigs += 1
            qname, flag, contig, pos, mapq, cigar, rnext, pnext, tlen, read, phred = myline.split("\t")

            contig_number, reference_genome, contig_start, contig_end = contig.split("_")
            contig_length = int(contig_end) - int(contig_start) + 1

            if contig != contig_name_temporary:
                create_new_fSAM(path_to_txt, contig_number, reference_genome, contig, contig_length, myline)

                nexons +=1
                contig_name_temporary = contig

            elif contig == contig_name_temporary:
                write_to_fSAM(path_to_txt, contig_number, reference_genome, myline)

    count_save_stats(path_to_txt, ncontigs, nexons)
