# extract_hits_psl.py
# This script is to extract the highest hits from psl files after BLAT
# Made by: Elfy Ly
# Date: 22 May 2020

from Bio import SearchIO

# def calculate_sizeMul(blat_hsp):
#     if blat_hsp.query_is_protein:
#         sizeMul = 3
#     else:
#         sizeMul = 1
#     return sizeMul

# moet loopen over contig nummers
path_to_psl = "./results/blat/SRR8528336/contig1_AT.psl"
path_to_fout = "./results/blat/SRR8528336/highest_hits.txt"

qresult = SearchIO.read(path_to_psl, "blat-psl")
# print(qresult)
# print("length: " + str(len(qresult[0])) + "\n")
print(qresult)
n_hits = len(qresult)
print(n_hits)
n_hsp_hits = len(qresult[0])

t_name_highest = ""
q_name_highest = ""
id_pct_highest = 0
max_score = 0
max_hit = 0
for hit in range(n_hits):           # loops over hits
    for n in range(n_hsp_hits):      # loops over hsp
        blat_hsp = qresult[hit][n]
        hit_number = n+1
        # print("hit: " + str(n+1))
        # print(blat_hsp.hit_id)

        # calculate_sizeMul(blat_hsp)

        # print("sizeMul: " + str(sizeMul))
        # print("score: " + str(blat_hsp.score))
        # print("percent ID: " + str(blat_hsp.ident_pct) + "\n")

        if blat_hsp.score > max_score:
            t_name_highest = qresult.id
            q_name_highest = blat_hsp.hit_id
            id_pct_highest = blat_hsp.ident_pct
            max_score = blat_hsp.score
            max_hit = hit_number

    print("tname: " + t_name_highest)
    print("best_exon: " + q_name_highest)
    print("highest hsp pct ID: " + str(id_pct_highest))
    print("max_hsp_hit: " + str(max_hit))
    print("max_hsp_score: " + str(max_score))

f = open(path_to_fout, "w+")
f.write("tName\tqName\tpercent ID\tPSL score\n")
f.write(t_name_highest + "\t" + q_name_highest + "\t" + str(id_pct_highest) + "\t" + str(max_hit) + "\t" +
        str(max_score))
f.close()





