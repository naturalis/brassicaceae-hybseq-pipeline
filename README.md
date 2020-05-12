# brassicaceae-hybseq-pipeline

This bioinformatics pipeline is used for the phylogenetic reconstruction of Brassicaceae using hybrid sequencing data.

For Windows 10: 
 - make sure miniconda3 is installed for Linux64 platform:
`$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

 - bash commands can be written in Anaconda prompt (miniconda3) or another program such as Ubuntu. 

---

1) Download the 2 FastQ (FTP) (fastq.gz) files of SRR8528336 in BioProject PRJNA518905 from:
https://www.ebi.ac.uk/ena/data/view/SRX5331770 by:
a) Create the path by:
`$ mkdir data/raw_reads/SRR8528336/`
b) Go to this path by:
`$ cd data/raw_reads/SRR8528336/`
c) Download the files by running the command:
`$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR852/006/SRR8528336/SRR8528336_1.fastq.gz`
`$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR852/006/SRR8528336/SRR8528336_2.fastq.gz` 
d) At the end you should get:
./data/raw_reads/SRR8528336/SRR8528336_1.fastq.gz and ./data/raw_reads/SRR8528336/SRR8528336_2.fastq.gz 
These are paired end reads: forward (1) and reverse (2)
e) Go back to the main folder.

2) Download the sequence genome Arabidopsis lyrata subsp. lyrata (GCA_000004255.1) from:
https://www.ncbi.nlm.nih.gov/assembly/GCF_000004255.2
a) By creating a new path:
`$ data/sequence_genomes`
b) Go to this path:
`$ cd data/sequence_genomes/download`
c) Download list of all GenBank files:
`$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt`
d) Get the download link of the correct genome by:
`$ grep -E 'GCF_000004255.2' assembly_summary_genbank.txt | cut -f 20 > ftp_folder.txt`
e) Create folder and script to save the links:
`awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}' ftp_folder.txt > download_fna_files.sh`
f) Download the sequence genome by:
`$ source download_fna_files.sh`
g) Unpack the .gz file by:
`$ gunzip GCA_000004255.1_v.1.0_genomic.fna.gz`
h) Change the .fna file to arl_ref.fa file and paste this file back to the sequence_genomes directory by:
`$ cp GCA_000004255.1_v.1.0_genomic.fna ../arl_ref.fa`

3) Install and download the necessary software packages for this project written in brassicaceae-hybseq-pipeline.yaml. 
Create the environment by running the following command in the main folder:
`$ conda env create -f ./envs/brassicaceae-hybseq-pipeline.yaml`

4) Go to the ./src folder

5) Install alignreads folder directory by running: 
`$ git clone https://github.com/zachary-foster/alignreads`

6) Install the software packages by running the command:
`$ python ./alignreads/install.py ./src/installed_alignreads`
choose the recommended versions:
YASRA-2.33.tar.gz
lastz-1.03.02.tar.gz
mummer 3.23

7) Go back to the main folder

8) To only use 'alignreads' instead of calling the python file, run the command:
`$ export PATH="$PATH:./src/installed_alignreads/alignreads"`
This is necessary for the use of Snakemake

9) In the main folder where the Snakefile is, run the Snakefile with the command:
`$ snakemake`


