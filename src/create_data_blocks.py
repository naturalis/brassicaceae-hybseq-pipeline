import os
from Bio import SeqIO


path_to_fS1L1M1R = "./envs/S1L1M1R.yaml"
subdata_list = []
fsubdata = open(path_to_fS1L1M1R, "rt")
for line in fsubdata:
    if not line.startswith("exons"):
        symbol, exon_name = line.split("-")
        exon_name = exon_name.strip()
        subdata_list.append(exon_name)

# for exon in subdata_list:
#     #print(subdata_list)
#     print(exon)

path_to_MSA_dir = "./results/A14_trimmed_prot/"
list_in_MSA_dir = os.listdir(path_to_MSA_dir)
pattern = "_NT.fasta"
subdata_dict = {}

for exon in subdata_list:
    ffasta = exon + pattern
    if ffasta in list_in_MSA_dir:
        path_to_ffasta = path_to_MSA_dir + ffasta
        for record in SeqIO.parse(path_to_ffasta, "fasta"):
            seq_length = len(record.seq)
            subdata_dict[exon] = seq_length

start_pos1 = 1
start_pos2 = 2
start_pos3 = 3
end = 0
path_to_fdata_blocks = "./results/A15_S1L1M1R_RAxML/data_blocks.txt"
fdata_blocks = open(path_to_fdata_blocks, "w+")
for exon_name in subdata_list:
    if exon_name in subdata_dict:
        seq_len = subdata_dict[exon_name]
        end += seq_len
        exon_name = exon_name.replace(".", "_")
        exon_name = exon_name.replace("@", "_")
        fdata_blocks.write(exon_name + "_pos1 = " + str(start_pos1) + "-" + str(end) + "\\3;\n")
        fdata_blocks.write(exon_name + "_pos2 = " + str(start_pos2) + "-" + str(end) + "\\3;\n")
        fdata_blocks.write(exon_name + "_pos3 = " + str(start_pos3) + "-" + str(end) + "\\3;\n")
        start_pos1 += seq_len
        start_pos2 += seq_len
        start_pos3 += seq_len


