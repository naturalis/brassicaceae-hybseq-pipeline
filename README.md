# brassicaceae-hybseq-pipeline

This bioinformatics pipeline is used for the phylogenetic reconstruction of Brassicaceae using hybrid sequencing data.

For Windows 10: 
 - make sure miniconda3 is installed for Linux64 platform.
 - bash commands can be written in Anaconda prompt (miniconda3) or another program such as Ubuntu. 

---

1) Download the 2 FastQ (FTP) (fastq.gz) files of SRR8528336 in BioProject PRJNA518905 from:
https://www.ebi.ac.uk/ena/data/view/SRX5331770
and put the files in the path: ~/data/raw_reads/SRR8528336/
so that: ./data/raw_reads/SRR8528336/SRR8528336_1.fastq.gz and ./data/raw_reads/SRR8528336/SRR8528336_2.fastq.gz 
These are paired end reads: forward (1) and reverse (2)

2) Install and download the necessary software packages for this project written in brassicaceae-hybseq-pipeline.yaml. 
Create the environment by running the following command in the main folder:
`$ conda env create -f ./envs/brassicaceae-hybseq-pipeline.yaml`

3) Go to the ./src folder

4) Install alignreads folder directory by running: 
`$ git clone https://github.com/zachary-foster/alignreads`

5) Install the software packages by running the command:
`$ python ./alignreads/install.py ./src/installed_alignreads`
choose the recommended versions:
YASRA-2.33.tar.gz
lastz-1.03.02.tar.gz
mummer 3.23

6) Go back to the main folder

7) To only use 'alignreads' instead of calling the python file, run the command:
`$ export PATH="$PATH:./src/installed_alignreads/alignreads"`
This is necessary for the use of Snakemake

8) In the main folder where the Snakefile is, run the Snakefile with the command:
`$ snakemake`


