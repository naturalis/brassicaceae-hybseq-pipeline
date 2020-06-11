# variables for every species
SAMPLES = ["SRR8528336", "SRR8528337"]

#forward or reverse | paired or unpaired
FRPU = ["forward_trim_paired", "forward_trim_unpaired", "reverse_trim_paired", "reverse_trim_unpaired"]

# variables within contigs: should be range(1, ALL CONTIG FILES IN DIR (FOR EVERY SPECIES))
CONTIGS_IN_DIR = 2113   #1905 for SRR8528337
CONTIG_NRS = range(1, CONTIGS_IN_DIR+1)   #["1", "2", "3"]

# all variables within snakemake
# raw_reads_forward_samples = expand("data/raw_reads/{samples}_1.fastq.gz", samples = SAMPLES)
# raw_reads_reverse_samples = expand("data/raw_reads/{samples}_2.fastq.gz", samples = SAMPLES)
# raw_reads_samples = expand("data/raw_reads/{samples}_count_reads.txt", samples = SAMPLES)
# deduplication_variables = expand("results/deduplicated_reads/SRR8528336/SRR8528336_{frpu}_dedupl.fq", frpu = FRPU)
# fsam_variables = expand("results/mapped_contigs/SRR8528336/sam/Contig{nr}_AT.sam", nr = CONTIG_NRS)
# fbam_variables = expand("results/mapped_contigs/SRR8528336/bam/Contig{nr}_AT.bam", nr = CONTIG_NRS)
# sorted_fbam_variables = expand("results/mapped_contigs/SRR8528336/sorted_bam/Contig{nr}_AT_sort.bam", nr = CONTIG_NRS)
# pileup_variables = expand("results/mapped_contigs/SRR8528336/pileup/Contig{nr}_AT_sort.pileup", nr = CONTIG_NRS)
# var_variables = expand("results/mapped_contigs/SRR8528336/var/Contig{nr}_AT_sort.var", nr = CONTIG_NRS)
blat_variables = expand("results/blat/SRR8528336/contig{nr}_AT.psl", nr = CONTIG_NRS)


rule all:
    input:
        # raw_reads_samples
        # raw_reads_forward_samples, raw_reads_reverse_samples, deduplication_variables, fsam_variables, fbam_variables,
        #  sorted_fbam_variables, pileup_variables, var_variables
        blat_variables

rule gunzip:
    input:
        "data/raw_reads/SRR8528336_1.fastq.gz",
        "data/raw_reads/SRR8528336_2.fastq.gz"
    shell:
        "gunzip {input}"

rule count_raw_reads:
    input:
        forward="data/raw_reads/{samples}_1.fastq",
        reverse="data/raw_reads/{samples}_2.fastq"
    output:
        "data/raw_reads/{samples}_count_reads.txt"
    shell:
        "echo $(cat {input.forward} | grep 'HISEQ' | wc -l) + $(cat {input.reverse} | grep 'HISEQ' | wc -l) | "
        "bc > {output}"

# preprocessing raw reads before alignment
rule trimming:
    input:
        "data/raw_reads/SRR8528336_1.fastq",
        "data/raw_reads/SRR8528336_2.fastq"
    output:
        expand("results/trimmed_reads/SRR8528336/SRR8528336_{FRPU}.fq", FRPU = FRPU)
    shell:
        "trimmomatic PE -phred33 {input} {output} "
        "ILLUMINACLIP:trimmomatic_adapter/TruSeq3-PE-2.fa:2:30:10 "
        "LEADING:20 "
        "TRAILING:20 "
        "SLIDINGWINDOW:5:20 "
        "MINLEN:36"

rule count_reads_trimming:
    input:
        expand("results/trimmed_reads/SRR8528336/SRR8528336_{FRPU}.fq", FRPU = FRPU)
    output:
        "results/trimmed_reads/SRR8528336/SRR8528336_count_reads.txt"
    shell:
        "echo $(cat {input} | wc -l)/4|bc >> {output}"

rule deduplication:
    input:
        "results/trimmed_reads/SRR8528336/SRR8528336_{frpu}.fq"
    output:
        "results/deduplicated_reads/SRR8528336/SRR8528336_{frpu}_dedupl.fq"
    shell:
        "fastx_collapser -v -i {input} -o {output}"

rule combine:
    input:
        expand("results/deduplicated_reads/SRR8528336/SRR8528336_{FRPU}_dedupl.fq", FRPU = FRPU)
    output:
        "results/deduplicated_reads/SRR8528336/SRR8528336_reads.fq"
    shell:
         "cat {input} > {output}"

rule count_reads_deduplication:
    input:
        "results/deduplicated_reads/SRR8528336/SRR8528336_reads.fq"
    output:
        "results/deduplicated_reads/SRR8528336/SRR8528336_count_reads.txt"
    shell:
        "grep '>' {input} | wc -l > {output}"

# reference mapping and de novo using YASRA/alignreads.py
# make sure to: $ export PATH="$PATH:~/usr/local/src/alignreads/alignreads"
rule alignreads:
    input:
        "results/deduplicated_reads/SRR8528336/SRR8528336_reads.fq",
        "data/reference_genomes/ref-at.fasta"
    output:
        "results/align"
    shell:
        "alignreads {input} "
        "--single-step "
        "--read-type solexa "
        "--read-orientation linear "
        "--percent-identity medium "
        "--depth-position-masking 5- "
        "--proportion-base-filter 0.7-"
        "--output-directory {output}"

# extract contigs from created SAM file after YASRA and create per contig new SAM files with headers
rule extract_contigs:
    input:
        "src/extract_contigs_YASRA.py"
    shell:
        "python3 {input}"

rule convert_to_fSAM:
    input:
        "results/mapped_contigs/SRR8528336/txt/Contig{nr}_AT.txt"
    output:
        temp("results/mapped_contigs/SRR8528336/sam/Contig{nr}_AT.sam")
    shell:
        "cp {input} {output}"

rule convert_to_fBAM:
    input:
        "results/mapped_contigs/SRR8528336/sam/Contig{nr}_AT.sam"
    output:
        temp("results/mapped_contigs/SRR8528336/bam/Contig{nr}_AT.bam")
    shell:
        "samtools view -bS {input} > {output}"

rule sort_fBAM:
    input:
        "results/mapped_contigs/SRR8528336/bam/Contig{nr}_AT.bam"
    output:
        temp("results/mapped_contigs/SRR8528336/sorted_bam/Contig{nr}_AT_sort.bam")
    shell:
        "samtools sort -m5G {input} -o {output}"

rule convert_to_fpileup:
    input:
        "results/mapped_contigs/SRR8528336/sorted_bam/Contig{nr}_AT_sort.bam"
    output:
        temp("results/mapped_contigs/SRR8528336/pileup/Contig{nr}_AT_sort.pileup")
    shell:
        "samtools mpileup -B {input} > {output}"

rule SNP_calling:
    input:
        "results/mapped_contigs/SRR8528336/pileup/Contig{nr}_AT_sort.pileup"
    output:
        "results/mapped_contigs/SRR8528336/var/Contig{nr}_AT_sort.var"
    shell:
        "varscan pileup2cns {input} "
        "--min-freq-for-hom 0.6 "
        "--min-coverage 5 "
        "--min-var-freq 0.6 "
        "--p-value 0.1 "
        "--min-reads2 5 "
        "> {output}"

rule make_consensus:
    input:
        "src/read_var.py"
    shell:
        "python3 {input}"

rule BLAT_assembled:
    input:
        "data/exons/exons_AT.fasta",
        "results/consensus/SRR8528336/Contig{nr}.txt"
    output:
        "results/blat/SRR8528336/contig{nr}_AT.psl"
    shell:
        "blat "
        "-t=dnax "
        "-q=dnax "
        "-stepSize=5 "
        "-repMatch=2253 "
        "-minScore=0 "
        "-minIdentity=0 "
        "{input} {output}"

rule extract_hits_psl:
    input:
        "src/extract_hits_psl.py"
    shell:
        "python3 {input}"


# SEQUENCE GENOMES
# extract reads from sequence genomes
rule extract_reads:
    input:
        "src/retrieve_exons_sequence_genomes.py"
    shell:
        "python {input}"

# identify contigs from sequence genomes
rule blat_sequence_genomes:
    input:
        "data/reference_genomes/ref-at.fasta",
        "data/sequence_genomes/sis_ref/KE154134.1.txt"
    output:
        "results/blat/at_sis_KE154134.1.psl"
    shell:
         "blat "
         "-t=dnax "
         "-q=dnax "
         "-stepSize=5 "
         "-repMatch=2253 "
         "-minScore=20 "
         "-minIdentity=0 "
         "{input} {output}"
