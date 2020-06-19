# brassicaceae-hybseq-pipeline
This is obtained by running the command:
$ git clone https://github.com/naturalis/brassicaceae-hybseq-pipeline.git

This bioinformatics pipeline is used for the phylogenetic reconstruction of Brassicaceae using hybrid sequencing data.

For Windows 10: 
 - make sure miniconda3 is installed for Linux64 platform:
`$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

 - bash commands can be written in Anaconda prompt (miniconda3) or another program such as Ubuntu. 

---

## Get raw reads
Download the 2 FastQ (FTP) (fastq.gz) files of the samples you like from https://www.ebi.ac.uk/ena/browse.
For example: SRR8528336 in BioProject PRJNA518905 from https://www.ebi.ac.uk/ena/data/view/SRR8528336 by:
1) Go to the data dir by:
`$ cd data/`
2) Create the raw_reads directories by:
`$ mkdir raw_reads`
3) Go to this path by:
`$ cd raw_reads`
4) Get the .fastq.gz files for the forward (1) and reverse (2) paired end reads by copying the link adress of FASTQ files (FTP). 
Download the files by running the command:
`$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR852/006/SRR8528336/SRR8528336_1.fastq.gz`
`$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR852/006/SRR8528336/SRR8528336_2.fastq.gz` 
5) Unpack these files by:
`$ gunzip SRR8528336_1.fastq.gz` 
`$ gunzip SRR8528336_2.fastq.gz`
6) Go back to the main folder
`$ cd ~/brassicaceae-hybseq-pipeline/`

## Get sequenced genome
Download the sequence genome Arabidopsis lyrata subsp. lyrata (GCA_000004255.1) from:
https://www.ncbi.nlm.nih.gov/assembly/GCF_000004255.2
1) By creating a new path:
`$ data/sequence_genomes`
2) Go to this path:
`$ cd data/sequence_genomes/download`
3) Download list of all GenBank files:
`$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt`
4) Get the download link of the correct genome and save this in ftp_folder.txt by:
`$ grep -E 'GCF_000004255.2' assembly_summary_genbank.txt | cut -f 20 > ftp_folder.txt`
5) Create folder and script to save the links:
`awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}' ftp_folder.txt > download_fna_files.sh`
6) Download the sequence genome by:
`$ source download_fna_files.sh`
7) Unpack the .gz file by: (is not necessary)
`$ gunzip GCA_000004255.1_v.1.0_genomic.fna.gz`
8) Change the .fna file to arl_ref.fa file and move this file back to the sequence_genomes directory by:
`$ mv GCA_000004255.1_v.1.0_genomic.fna ../arl_ref.fa`

## Create and activate environment
1) Install and download the necessary software packages for this project written in brassicaceae-hybseq-pipeline.yaml. 
Create the environment by running the following command in the main folder:
`$ conda env create -f ./envs/brassicaceae-hybseq-pipeline.yaml`
2) Activate this environment by:
`$ conda activate brassicaceae-hybseq-pipeline`

## Download alignreads.py
(alignreads.py is currently installed in /usr/local/src/alignreads)
1) Go to the ./src folder
2) Install alignreads folder directory by running: 
`$ git clone https://github.com/zachary-foster/alignreads`
3) Install the software packages by running the command:
`$ python ./alignreads/install.py ./installed_alignreads`
4) choose the recommended versions:
(8) lastz-1.03.02.tar.gz
(1) YASRA-2.33.tar.gz
(1) mummer 3.23
5) Go back to the main folder
6) To only use 'alignreads' instead of calling the python file, run the command:
`$ export PATH="$PATH:./src/installed_alignreads/alignreads"`
This is necessary for the use of Snakemake
(at this moment: $ export PATH="$PATH:~/usr/local/src/alignreads/alignreads)


## Run snakemake
1) In the main folder where the Snakefile is, run the Snakefile with the command:
`$ snakemake`


