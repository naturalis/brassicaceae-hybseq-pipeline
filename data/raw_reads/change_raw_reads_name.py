# This script is to change the name of raw reads samples in /data/raw_reads if the files are ending with .fastq.gz
# This script is specifically made for files sequenced by Naturalis, for example:
# S0275_Brassi-hy01_61237-5_GACCAAGTTAATATAAGGCG_L001_R1_001_BHMWGWDRXX.filt.fastq.gz
# Made by: Elfy Ly
# Date: 25 August 2020

import os
from os import listdir
from os.path import isfile, join

current_dir = os.path.abspath(os.getcwd())
onlyfiles = [f for f in listdir(current_dir) if isfile(join(current_dir, f))]

test = []
for file in onlyfiles:
    if file.endswith(".fastq.gz"):
        sample_name, brassi, number, primer, l, direction, nr, extention = file.split("_")
        # sample_name, direction = file.split("_")
        new_name = ""
        if "1" in direction:
            new_name = sample_name + "_1.fastq.gz"
        elif "2" in direction:
            new_name = sample_name + "_2.fastq.gz"
        os.rename(file, new_name)
