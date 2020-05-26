# SAMPLES = ["SRR8528337"] # variables for every species

#forward or reverse | paired or unpaired
FRPU = ["forward_trim_paired", "forward_trim_unpaired", "reverse_trim_paired", "reverse_trim_unpaired"]

# variables within contigs: should be range(1, ALL CONTIG FILES IN DIR (FOR EVERY SPECIES))
CONTIG_NRS = range(1, 2114)   #["1", "2", "3"]

# all variables within snakemake
deduplication_variables = expand("results/deduplicated_reads/SRR8528337/SRR8528337_{frpu}_dedupl.fq", frpu = FRPU)
fsam_variables = expand("results/exons/sam/Contig{nr}_AT.sam", nr = CONTIG_NRS)
fbam_variables = expand("results/exons/bam/Contig{nr}_AT.bam", nr = CONTIG_NRS)
sorted_fbam_variables = expand("results/exons/sorted_bam/Contig{nr}_AT_sort.bam", nr = CONTIG_NRS)
pileup_variables = expand("results/exons/pileup/Contig{nr}_AT_sort.pileup", nr = CONTIG_NRS)
var_variables = expand("results/exons/var/Contig{nr}_AT_sort.var", nr = CONTIG_NRS)


rule all:
    input:
        deduplication_variables, fsam_variables, fbam_variables, sorted_fbam_variables, pileup_variables, var_variables


# preprocessing raw reads before alignment
rule trimming:
    input:
        "data/raw_reads/SRR8528337_1.fastq.gz",
        "data/raw_reads/SRR8528337_2.fastq.gz"
    output:
        expand("results/trimmed_reads/SRR8528337/SRR8528337_{FRPU}.fq", FRPU = FRPU)
    shell:
        "trimmomatic PE -phred33 {input} {output} "
        "ILLUMINACLIP:trimmomatic_adapter/TruSeq3-PE-2.fa:2:30:10 "
        "LEADING:20 "
        "TRAILING:20 "
        "SLIDINGWINDOW:5:20 "
        "MINLEN:36"

rule count_reads_trimming:
    input:
        expand("results/trimmed_reads/SRR8528337/SRR8528337_{FRPU}.fq", FRPU = FRPU)
    output:
        "results/trimmed_reads/SRR8528337/SRR8528337_count_reads.txt"
    shell:
        "echo $(cat {input} | wc -l)/4|bc >> {output}"

rule deduplication:
    input:
        "results/trimmed_reads/SRR8528337/SRR8528337_{frpu}.fq"
    output:
        "results/deduplicated_reads/SRR8528337/SRR8528337_{frpu}_dedupl.fq"
    shell:
        "fastx_collapser -v -i {input} -o {output}"

rule combine:
    input:
        expand("results/deduplicated_reads/SRR8528337/SRR8528337_{FRPU}_dedupl.fq", FRPU = FRPU)
    output:
        "results/deduplicated_reads/SRR8528337/SRR8528337_reads.fq"
    shell:
         "cat {input} > {output}"

rule count_reads_deduplication:
    input:
        "results/deduplicated_reads/SRR8528337/SRR8528337_reads.fq"
    output:
        "results/deduplicated_reads/SRR8528337/SRR8528337_count_reads.txt"
    shell:
        "grep '>' {input} | wc -l > {output}"

# reference mapping and de novo using YASRA/alignreads.py
# make sure to: $ export PATH="$PATH:~/usr/local/src/alignreads/alignreads"
# this rule to be changed into multiple rules within the alignreads.py
rule alignreads:
    input:
        "results/deduplicated_reads/SRR8528337/SRR8528337_reads.fq",
        "data/reference_genomes/ref-at.fasta"
        #"src/installed_alignreads/alignreads/YASRA-2.33/test_data/454.fa",
        #"src/installed_alignreads/alignreads/YASRA-2.33/test_data/rhino_template.fa"
    shell:
        "alignreads {input} "
        "--single-step "
        "--read-type solexa "
        "--read-orientation linear "
        "--percent-identity medium "
        "--depth-position-masking 5- "
        "--proportion-base-filter 0.7-"

# extract contigs from created SAM file after YASRA and create per contig new SAM files with headers
rule extract_contigs:
    input:
        "src/extract_contigs_YASRA.py"
    shell:
        "python3 {input}"

rule convert_to_fSAM:
    input:
        "results/exons/txt/Contig{nr}_AT.txt"
    output:
        temp("results/exons/sam/Contig{nr}_AT.sam")
    shell:
        "cp {input} {output}"

rule convert_to_fBAM:
    input:
        "results/exons/sam/Contig{nr}_AT.sam"
    output:
        temp("results/exons/bam/Contig{nr}_AT.bam")
    shell:
        "samtools view -bS {input} > {output}"

rule sort_fBAM:
    input:
        "results/exons/bam/Contig{nr}_AT.bam"
    output:
        temp("results/exons/sorted_bam/Contig{nr}_AT_sort.bam")
    shell:
        "samtools sort -m5G {input} -o {output}"

rule convert_to_fpileup:
    input:
        "results/exons/sorted_bam/Contig{nr}_AT_sort.bam"
    output:
        temp("results/exons/pileup/Contig{nr}_AT_sort.pileup")
    shell:
        "samtools mpileup -B {input} > {output}"

rule SNP_calling:
    input:
        "results/exons/pileup/Contig{nr}_AT_sort.pileup"
    output:
        "results/exons/var/Contig{nr}_AT_sort.var"
    shell:
        "varscan pileup2cns {input} "
        "--min-freq-for-hom 0.6 "
        "--min-coverage 5 "
        "--min-var-freq 0.6 "
        "--p-value 0.1 "
        "--min-reads2 5 "
        "> {output}"


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
