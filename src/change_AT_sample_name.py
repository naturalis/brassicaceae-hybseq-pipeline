import os
from Bio import SeqIO

path_to_trim_dir = "A14_trimmed_prot"
path_to_new_dir = "A14b_new"



listed_dir = os.listdir(path_to_trim_dir)
for file in listed_dir:
    if file.endswith("_NT.fasta"):
        exon, extension = file.split("_NT.fasta")
        path_to_foriginal = path_to_trim_dir + "/" + file
        path_to_fnew = path_to_new_dir + "/" + file
        with open(path_to_foriginal) as original, open(path_to_fnew, 'w') as corrected:
            records = SeqIO.parse(path_to_foriginal, 'fasta')
            for record in records:
                if record.id == exon:
                    #print(record.id)
                    record.id = "Arabidopsis"
                    record.description = "Thaliana"
                #print(record.id)
                SeqIO.write(record, corrected, 'fasta')



