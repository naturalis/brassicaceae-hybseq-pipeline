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


# creates a list of mapped contigs and their start and end position
def create_mapped_contig_list(path_to_final_assembly, path_to_mapped_contigs):
    f = open(path_to_mapped_contigs, "w+")
    print("New text file created: " + path_to_species + "mapped_contigs.txt")

    with open(path_to_final_assembly, 'rt') as myfile:
        for myline in myfile:
            if myline.startswith('>'):
                contig_name, ref_name, contig_start, contig_end = myline.split('_')
                symbol, contig_name = contig_name.split('>')
                f.write(contig_name + "\t" + contig_start + "\t" + contig_end)
    f.close()


def create_new_fSAM(path_to_txt, contig_number, reference_genome, contig, contig_length, myline):
    f = open(path_to_txt + contig_number + "_" + reference_genome + ".txt", "w+")
    print("New text file created: " + contig_number + "_" + reference_genome)
    f.write("@HD\tVN:1.3\n@SQ\tSN:" + contig + "\tLN:" + str(contig_length) + "\n" + myline)
    f.close()


def write_to_fSAM(path_to_txt, contig_number, reference_genome, myline):
    f = open(path_to_txt + contig_number + "_" + reference_genome + ".txt", "a+")
    f.write(myline)
    f.close()


def count_save_stats(path_to_txt, ncontigs, nexons):
    f = open(path_to_txt + "number_of_contigs_and_exons.txt", "w+")
    f.write("Number of contigs: " + str(ncontigs) + "\nNumber of exons: " + str(nexons) + "\n")
    f.close()


# Code starts here
# Creating paths for output directory
# Deze paths moeten nog automatisch laten loopen, maar daarvoor moet eerst deze output directory naam gewijzigd worden
path_to_YASRA_dir = './results/alignments/SRR8528336_reads.fq_ref-at.fasta_Mon-May-25-08:50:00-2020/YASRA_related_files/'
path_to_YASRA_fSAM = path_to_YASRA_dir + 'alignments_SRR8528336_reads.fq_ref-at.fasta.sam'
path_to_assembled_exons = './results/assembled_exons/'
path_to_species = path_to_assembled_exons + 'SRR8528336/'
path_to_txt = path_to_species + 'txt/'

create_dir(path_to_assembled_exons)
create_dir(path_to_species)
create_dir(path_to_txt)

# Creates assembled contig list and dictionary for start and end position
path_to_final_assembly = path_to_YASRA_dir +'Final_Assembly_SRR8528336_reads.fq_ref-at.fasta'
path_to_mapped_contigs = path_to_species + "mapped_contigs.txt"

create_mapped_contig_list(path_to_final_assembly, path_to_mapped_contigs)

# write new SAM files for every contig after YASRA
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
