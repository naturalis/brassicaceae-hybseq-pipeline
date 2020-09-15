# brassicaceae-hybseq-pipeline
This is obtained by running the command:
`$ git clone https://github.com/naturalis/brassicaceae-hybseq-pipeline.git`

Go to this directory by running the command:
`$ cd brassicaceae-hybseq-pipeline`

This bioinformatics pipeline is used for the phylogenetic reconstruction of Brassicaceae using hybrid sequencing data.

For Windows 10: 
 - make sure miniconda3 is installed for Linux64 platform:
`$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

 - bash commands can be written in Anaconda prompt (miniconda3) or another program such as Ubuntu. 

---
## Create and activate environment
This conda environment contains all tools and versions necessary for the execution of this pipeline
1) Install and download the necessary software packages for this project written in brassicaceae-hybseq-pipeline.yaml. 
Create the environment by running the following command in the main folder:
- `$ conda env create -f ./envs/brassicaceae-hybseq-pipeline.yaml`
2) Activate this environment by:
- `$ conda activate brassicaceae-hybseq-pipeline`

## Get raw reads 
##### From BioProject:
Download the 2 FastQ (FTP) (fastq.gz) files of the samples you like from https://www.ebi.ac.uk/ena/browse.
For example: SRR8528336 in BioProject PRJNA518905 from https://www.ebi.ac.uk/ena/data/view/SRR8528336 by:

1) Go to the data/raw_reads dir by:
`$ cd data/raw_reads`
2) Get the .fastq.gz files for the forward (1) and reverse (2) paired end reads by copying the link adress of FASTQ files (FTP). 
Download the files by running the command:
- `$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR852/006/SRR8528336/SRR8528336_1.fastq.gz`
- `$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR852/006/SRR8528336/SRR8528336_2.fastq.gz` 
3) Unpack these files by:
- `$ gunzip SRR8528336_1.fastq.gz` 
- `$ gunzip SRR8528336_2.fastq.gz`
4) Go back to the main folder by: `$ cd ~/brassicaceae-hybseq-pipeline/`

##### From elsewhere:
1) Go to the data/raw_reads dir by:
`$ cd data/raw_reads`
2) Download the .fastq.gz files in the path: brassicaceae-hybseq-pipeline/data/raw_reads/
3) Make sure the name of the files look like: SAMPLECODE_1.fastq.gz and SAMPLECODE_2.fastq.gz.
- If sequenced samples are sequenced by Naturalis, change_raw_reads_name.py can be ran to change all names at once:
- `$ python change_raw_reads_name.py`
4) Unpack these files one by one by:
- `$ gunzip [SAMPLECODE]_1.fastq.gz` 
- `$ gunzip [SAMPLECODE]_2.fastq.gz`
Or all at once:
- `$ gunzip *.fastq.gz`
5) Go back to the main folder by: `$ cd ~/brassicaceae-hybseq-pipeline/`

## Download alignreads.py
##### In Naturalis high-mem:
1) alignreads.py is installed in /usr/local/src/alignreads
2) To only use 'alignreads' instead of calling the python file (necessary for the execution of Snakemake), run the command:
- `$ export PATH="$PATH:~/usr/local/src/alignreads/alignreads"`

##### If alignreads.py has not been installed yet:
1) Go to the ./src folder
2) Download alignreads folder directory by running: 
- `$ git clone https://github.com/zachary-foster/alignreads`
3) Install the software packages by running the command:
- `$ python ./alignreads/install.py ./installed_alignreads`
4) choose the recommended versions:
(8) lastz-1.03.02.tar.gz
(1) YASRA-2.33.tar.gz
(1) mummer 3.23
5) Go back to the main folder: `$ cd ..`
6) To only use 'alignreads' instead of calling the python file (necessary for the execution of Snakemake), run the command:
- `$ export PATH="$PATH:./src/installed_alignreads/alignreads"`

## Download phyutilities v2.7.1:
For the concatenation step in Snakemake_A5_RAxML
1) Go to the ./src folder
2) Donwload phyutilities by running:
- `$ wget https://github.com/blackrim/phyutility/releases/download/v2.7.1/phyutility_2.7.1.tar.gz`
3) Unzip the .gz file by running:
- `$ gunzip phyutility_2.7.1.tar.gz` 
4) Unpack the .tar file by:
- `tar -xvf phyutility_2.7.1.tar`
5) Go back to the main folder: `$ cd ..`


## Run snakemake
1) Open the snakefile by running:
- `$ nano Snakefile_brassicaceae-hybseq-pipeline`
2) Adjust in the Snakefiles (A1, A2 and A3) for the samples of interest (can be one or multiple) by changing the sample names of 'SAMPLES'.
For example: SAMPLES = "S0603 S0560".split() to SAMPLES = "SAMPLES = "SRR8528336".split()
3) Adjust in Snakefile_A1_Mapping the ORIGIN variable to "naturalis", "donovan" or "nikolov".
For example: ORIGIN = "naturalis" to ORIGIN = "donovan"
(This is for the count of raw reads in their FASTQ file. Counting the number of raw reads is based on the sequence identifier which differs per sequencer: Naturalis @A00, Nikolov @SRR and donovan @M01)
4) To run the parallel pipeline, in the main folder where the Snakefile is, run the Snakefile with the commands in the order:
- `$ snakemake --snakefile Snakefile_A1_Mapping -F part1`
- `$ snakemake --snakefile Snakefile_A1_Mapping -F part2`
- `$ snakemake --snakefile Snakefile_A2_Contigs -F part3`
- `$ snakemake --snakefile Snakefile_A2_Contigs -F part4`
- `$ snakemake --snakefile Snakefile_A3_Exons -F part5`
- `$ snakemake --snakefile Snakefile_A3_Exons -F part6`
5) To run the single pipeline, make sure the path ./results/A09_consensus_exons only contains the samples of interest and that at least one sample covers all exons to make Snakefile_A4_MSA work.
- `$ snakemake --snakefile Snakefile_A3_Exons -F merge_exon_seqs`
- `$ snakemake --snakefile Snakefile_A4_MSA -F all`
- `$ snakemake --snakefile Snakefile_A5_RAxML -F phyutilities`

## Build phylogenetic tree
6) The concatenated MSA output can be found in the path: results/A15_S1L1M1R_RAxML/S1L1M1R_msa.aln
This nexus file should be converted to a phy file (I used Geneious).
7) To select the best model using PartitionFinder, datablocks (how the genes are splitted in all 3 codon positions) should be given as input. 
This can be calculated by running: 
- `$ python src/create_data_blocks.py`
The output will be in the path: results/A15_S1L1M1R_RAxML/datablocks.txt. This data should be copied to the partition_finder.cfg configuration file (follow the PartitionFinder manual). The output looks like this:

>AT5G23110_1_6_pos1 = 1-219\3;
>AT5G23110_1_6_pos2 = 2-219\3;
>AT5G23110_1_6_pos3 = 3-219\3;
>AT5G23110_1_5_pos1 = 220-480\3;
>AT5G23110_1_5_pos2 = 221-480\3;
>...
>AT4G20090_1_1_pos3 = 527529-529296\3;
>AT2G40840_1_12_pos1 = 529297-529539\3;
>AT2G40840_1_12_pos2 = 529298-529539\3;
>AT2G40840_1_12_pos3 = 529299-529539\3;
8) Run PartitionFinder for model selection (follow PF2 manual)
9) Run RAxML for phylogenetic tree reconstruction (follow RAxML manual)


## Get sequenced genome
(Only necessary if pipeline for sequenced genomes is finished - currently work in progress).
Download the sequence genome Arabidopsis lyrata subsp. lyrata (GCA_000004255.1) from:
https://www.ncbi.nlm.nih.gov/assembly/GCF_000004255.2
1) By creating a new path:
`$ data/sequence_genomes`
2) Go to this path:
`$ cd data/sequence_genomes/download`
3) Download list of all GenBank files:
- `$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt`
4) Get the download link of the correct genome and save this in ftp_folder.txt by:
- `$ grep -E 'GCF_000004255.2' assembly_summary_genbank.txt | cut -f 20 > ftp_folder.txt`
5) Create folder and script to save the links:
- `awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}' ftp_folder.txt > download_fna_files.sh`
6) Download the sequence genome by:
`$ source download_fna_files.sh`
7) Unpack the .gz file by: 
`$ gunzip GCA_000004255.1_v.1.0_genomic.fna.gz`
8) Change the .fna file to arl_ref.fa file and move this file back to the sequence_genomes directory by:
- `$ mv GCA_000004255.1_v.1.0_genomic.fna ../arl_ref.fa`

