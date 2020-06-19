# variables for every species
SAMPLES = ["SRR8528336", "SRR8528337", "SRR8528338"]
SAMPLE_NAME = "SRR8528336"

#forward or reverse | paired or unpaired
FRPU = ["forward_trim_paired", "forward_trim_unpaired", "reverse_trim_paired", "reverse_trim_unpaired"]

# variables within contigs: should be range(1, ALL CONTIG FILES IN DIR (FOR EVERY SPECIES))
CONTIGS_MAPPED = 2146 #1905 for SRR8528337   #2113 for SRR8528336
CONTIG_NRS = range(1, CONTIGS_MAPPED + 1)   #["1", "2", "3"]

# configfile: "./envs/MAFFT/" + SAMPLE_NAME + ".yaml"

# all variables within snakemake
# raw_reads_samples = expand("data/raw_reads/{samples}_count_reads.txt", samples = SAMPLES)
deduplication_variables = expand("results/2_deduplicated_reads/"+ SAMPLE_NAME + "/"+ SAMPLE_NAME + "_{frpu}_dedupl.fq", frpu = FRPU)
var_variables = expand("results/4_mapped_contigs/"+ SAMPLE_NAME + "/var/Contig{nr}_AT_sort.var", nr = CONTIG_NRS)
blat_variables = expand("results/6_identified_contigs_blat/"+ SAMPLE_NAME + "/contig{nr}_AT.psl", nr = CONTIG_NRS)
# mafft_variables = expand("results/8_assembled_exons/"+ SAMPLE_NAME + "/{exon}.fasta", exon=config["exons"])


rule all:
    input:
        # deduplication_variables
        var_variables
        # blat_variables
        # mafft_variables

rule count_raw_reads:
    input:
        forward="data/raw_reads/"+ SAMPLE_NAME + "_1.fastq",
        reverse="data/raw_reads/"+ SAMPLE_NAME + "_2.fastq"
    output:
        "data/raw_reads/"+ SAMPLE_NAME + "_count_reads.txt"
    shell:
        "echo $(cat {input.forward} | grep '@SRR' | wc -l) + $(cat {input.reverse} | grep '@SRR' | wc -l) | "
        "bc > {output}"

# preprocessing raw reads before alignment
rule trimming:
    input:
        "data/raw_reads/"+ SAMPLE_NAME + "_1.fastq",
        "data/raw_reads/"+ SAMPLE_NAME + "_2.fastq"
    output:
        expand("results/1_trimmed_reads/"+ SAMPLE_NAME + "/"+ SAMPLE_NAME + "_{FRPU}.fq", FRPU = FRPU)
    shell:
        "trimmomatic PE -phred33 {input} {output} "
        "ILLUMINACLIP:trimmomatic_adapter/TruSeq3-PE-2.fa:2:30:10 "
        "LEADING:20 "
        "TRAILING:20 "
        "SLIDINGWINDOW:5:20 "
        "MINLEN:36"

rule count_reads_trimming:
    input:
        expand("results/1_trimmed_reads/"+ SAMPLE_NAME + "/"+ SAMPLE_NAME + "_{FRPU}.fq", FRPU = FRPU)
    output:
        "results/1_trimmed_reads/"+ SAMPLE_NAME + "/"+ SAMPLE_NAME + "_count_reads.txt"
    shell:
        "echo $(cat {input} | wc -l)/4|bc >> {output}"

rule deduplication:
    input:
        "results/1_trimmed_reads/"+ SAMPLE_NAME + "/"+ SAMPLE_NAME + "_{frpu}.fq"
    output:
        "results/2_deduplicated_reads/"+ SAMPLE_NAME + "/"+ SAMPLE_NAME + "_{frpu}_dedupl.fq"
    shell:
        "fastx_collapser -v -i {input} -o {output}"

rule combine:
    input:
        expand("results/2_deduplicated_reads/"+ SAMPLE_NAME + "/"+ SAMPLE_NAME + "_{FRPU}_dedupl.fq", FRPU = FRPU)
    output:
        "results/2_deduplicated_reads/"+ SAMPLE_NAME + "/"+ SAMPLE_NAME + "_reads.fq"
    shell:
         "cat {input} > {output}"

rule count_reads_deduplication:
    input:
        "results/2_deduplicated_reads/"+ SAMPLE_NAME + "/"+ SAMPLE_NAME + "_reads.fq"
    output:
        "results/2_deduplicated_reads/"+ SAMPLE_NAME + "/"+ SAMPLE_NAME + "_count_reads.txt"
    shell:
        "grep '>' {input} | wc -l > {output}"

# reference mapping and de novo using YASRA/alignreads.py
# make sure to: $ export PATH="$PATH:~/usr/local/src/alignreads/alignreads"
rule alignreads:
    input:
        "results/2_deduplicated_reads/"+ SAMPLE_NAME + "/"+ SAMPLE_NAME + "_reads.fq",
        "data/reference_genomes/ref-at.fasta"
    output:
        "results/3_mapped_reads/" + SAMPLE_NAME
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
        "python3 {input} " + SAMPLE_NAME

rule convert_to_fBAM:
    input:
        "results/4_mapped_contigs/"+ SAMPLE_NAME + "/sam/Contig{nr}_AT.sam"
    output:
        temp("results/4_mapped_contigs/"+ SAMPLE_NAME + "/bam/Contig{nr}_AT.bam")
    shell:
        "samtools view -bS {input} > {output}"

rule sort_fBAM:
    input:
        "results/4_mapped_contigs/"+ SAMPLE_NAME + "/bam/Contig{nr}_AT.bam"
    output:
        temp("results/4_mapped_contigs/"+ SAMPLE_NAME + "/sorted_bam/Contig{nr}_AT_sort.bam")
    shell:
        "samtools sort -m5G {input} -o {output}"

rule convert_to_fpileup:
    input:
        "results/4_mapped_contigs/"+ SAMPLE_NAME + "/sorted_bam/Contig{nr}_AT_sort.bam"
    output:
        "results/4_mapped_contigs/"+ SAMPLE_NAME + "/pileup/Contig{nr}_AT_sort.pileup"
    shell:
        "samtools mpileup -B {input} > {output}"

rule SNP_calling:
    input:
        "results/4_mapped_contigs/"+ SAMPLE_NAME + "/pileup/Contig{nr}_AT_sort.pileup"
    output:
        "results/4_mapped_contigs/"+ SAMPLE_NAME + "/var/Contig{nr}_AT_sort.var"
    shell:
        "varscan pileup2cns {input} "
        "--min-freq-for-hom 0.6 "
        "--min-coverage 5 "
        "--min-var-freq 0.6 "
        "--p-value 0.1 "
        "--min-reads2 5 "
        "> {output}"

rule make_contig_consensus:
    input:
        "src/read_var.py"
    shell:
        "python3 {input} " + SAMPLE_NAME

rule BLAT_assembled:
    input:
        "data/exons/exons_AT.fasta",
        "results/5_consensus_contigs/"+ SAMPLE_NAME + "/Contig{nr}.txt"
    output:
        "results/6_identified_contigs_blat/"+ SAMPLE_NAME + "/contig{nr}_AT.psl"
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
        "python3 {input} " + SAMPLE_NAME

rule MAFFT_assembly:
    input:
        "results/7_mapped_exons/"+ SAMPLE_NAME + "/{exon}.fasta"
    output:
        "results/8_assembled_exons/"+ SAMPLE_NAME + "/{exon}.fasta"
    shell:
        "mafft "
        "--maxiterate 1000 "
        "--oldgenafpair "
        "{input} > {output}"

rule make_exon_consensus:
    input:
        "src/make_exon_consensus.py"
    shell:
        "python3 {input} " + SAMPLE_NAME


''' SEQUENCE GENOMES '''
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
