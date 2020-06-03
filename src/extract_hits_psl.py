# extract_hits_psl.py
# This script is to extract the highest hits from psl files after BLAT
# Made by: Elfy Ly
# Date: 22 May 2020

from Bio import SearchIO

MIN_HITS = 5

path_to_psl = "/mnt/c/Users/elfyl/PycharmProjects/brassicaceae-hybseq-pipeline-offline/results/blat/at_sis/at_sis_KE154134.1.psl"          # 45 matches
path_to_psl_2 = "/mnt/c/Users/elfyl/PycharmProjects/brassicaceae-hybseq-pipeline-offline/results/blat/at_sis/at_sis_ASZH01000019.1.psl"    # should be empty
path_to_psl_3 = "/mnt/c/Users/elfyl/PycharmProjects/brassicaceae-hybseq-pipeline-offline/results/blat/at_sis/at_sis_KE154096.1.psl"        # 29 matches
path_to_psl_4 = "/results/blat/SRR8528336/contig1_AT.psl"

qresult = SearchIO.read(path_to_psl_3, "blat-psl")
# print(qresult)
# print("length: " + str(len(qresult[0])) + "\n")
print(qresult)
n_hits = len(qresult)
print(n_hits)
n_hsp_hits = len(qresult[0])

highest_id_pct = 0
max_score = 0
max_hit = 0
best_exon = ""
for n in range(n_hsp_hits):
    blat_hsp = qresult[0][n]      # first hit, loop over hsp
    hit_number = n+1
    # print("hit: " + str(n+1))
    # print(blat_hsp.hit_id)
    if blat_hsp.query_is_protein:
        sizeMul = 3
    else:
        sizeMul = 1
    # print("sizeMul: " + str(sizeMul))
    # print("score: " + str(blat_hsp.score))
    # print("percent ID: " + str(blat_hsp.ident_pct) + "\n")
    if blat_hsp.score > max_score:
        highest_id_pct = blat_hsp.ident_pct
        max_score = blat_hsp.score
        max_hit = hit_number
        best_exon = blat_hsp.hit_id

print("max_score: " + str(max_score))
print("max_hit: " + str(max_hit))
print("highest pct ID: " + str(highest_id_pct))
print("best_exon: " + best_exon)



