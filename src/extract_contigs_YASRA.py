# extract_contigs_YASRA.py
# This script is to extract contigs and to create SAM files with headers for each of them after running YASRA.
# YASRA creates a SAM file without header after alignments.
# This can be executed by running: $ python extract_contigs_YASRA.py [SAMPLE_NAME], for example:
# $ python extract_contigs_YASRA.py SRR8528336
# Made by: Elfy Ly
# Date: 7 May 2020

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


def create_new_fSAM(path_to_sam, contig_number, reference_genome, contig, contig_length, myline):
    f = open(path_to_sam + contig_number + "_" + reference_genome + ".sam", "w+")
    print("New text file created: " + contig_number + "_" + reference_genome)
    f.write("@HD\tVN:1.3\n@SQ\tSN:" + contig + "\tLN:" + str(contig_length) + "\n" + myline)
    f.close()


def write_to_fSAM(path_to_sam, contig_number, reference_genome, myline):
    f = open(path_to_sam + contig_number + "_" + reference_genome + ".sam", "a+")
    f.write(myline)
    f.close()


def count_save_stats(path_to_sam, nreads, ncontigs):
    f = open(path_to_sam + "number_of_reads_and_contigs.txt", "w+")
    f.write("Number of reads: " + str(nreads) + "\nNumber of contigs: " + str(ncontigs) + "\n")
    f.close()


def create_YAML(path):
    f = open(path, "w+")
    f.write("contig_nrs:\n")
    f.close()
    print(path + " is created")


'''Code starts here'''
# Change YASRA output dir name
path_to_mapped_reads_dir = './results/3_mapped_reads/'
path_to_mapped_reads_sample_dir = path_to_mapped_reads_dir + SAMPLE_NAME + "/"
create_dir(path_to_mapped_reads_dir)
create_dir(path_to_mapped_reads_sample_dir)

current_dir = os.path.abspath(os.getcwd())
list_current_dir = os.listdir(current_dir)
for dir in list_current_dir:
    if SAMPLE_NAME in dir:
        old_path_output = "./" + dir + "/"
        os.rename(old_path_output, path_to_mapped_reads_sample_dir)

# Creating paths for output directory
path_to_YASRA_dir = path_to_mapped_reads_sample_dir + "YASRA_related_files/"
path_to_YASRA_fSAM = path_to_YASRA_dir + 'alignments_' + SAMPLE_NAME + '_reads.fq_ref-at.fasta.sam'
path_to_mapped_contigs = './results/4_mapped_contigs/'
path_to_species = path_to_mapped_contigs + SAMPLE_NAME + '/'
path_to_sam = path_to_species + 'sam/'

create_dir(path_to_mapped_contigs)
create_dir(path_to_species)
create_dir(path_to_sam)

# Creates assembled contig list and dictionary for start and end position
path_to_final_assembly = path_to_YASRA_dir + 'Final_Assembly_' + SAMPLE_NAME + '_reads.fq_ref-at.fasta'
path_to_fmapped_contigs = path_to_species + "mapped_contigs.txt"

create_mapped_contig_list(path_to_final_assembly, path_to_fmapped_contigs)

# write new SAM files for every contig after YASRA
nreads = 0
ncontigs = 0
contig_name_temporary = " "
with open(path_to_YASRA_fSAM, 'rt') as myfile:
    for myline in myfile:
        if myline.startswith('>'):
            nreads += 1
            qname, flag, contig, pos, mapq, cigar, rnext, pnext, tlen, read, phred = myline.split("\t")

            contig_number, reference_genome, contig_start, contig_end = contig.split("_")
            contig_length = int(contig_end) - int(contig_start) + 1

            if contig != contig_name_temporary:
                create_new_fSAM(path_to_sam, contig_number, reference_genome, contig, contig_length, myline)

                ncontigs +=1
                contig_name_temporary = contig

            elif contig == contig_name_temporary:
                write_to_fSAM(path_to_sam, contig_number, reference_genome, myline)

    count_save_stats(path_to_sam, nreads, ncontigs)

# creates configuration YAML file in envs contigs dir
path_to_contigs_configs = "./envs/contigs/"
create_dir(path_to_contigs_configs)
path_to_fYAML = path_to_contigs_configs + SAMPLE_NAME + ".yaml"
create_YAML(path_to_fYAML)

for ncontig in range(1, ncontigs + 1):
    path_to_fvar = path_to_species + "var/Contig" + str(ncontig) + "_AT_sort.var"
    fYAML = open(path_to_fYAML, "a+")
    fYAML.write("    " + str(ncontig) + ": " + path_to_fvar + "\n")
    fYAML.close()





