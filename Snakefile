# mapping against reference using minimap2 instead of YASRA
rule minimap2:
    input:
        "data/reference_genomes/ref-at.fasta data/raw_reads/SRR8528336_1.fastq.gz",
        "data/raw_reads/SRR8528336_2.fastq.gz"
    output:
        "results/alignments/minimap2/SRR8528336_ref-at.sam"
    shell:
        "minimap2 -ax sr {input} > {output}"

rule convert_to_fBAM:
    input:
        "results/alignments/minimap2/SRR8528336_ref-at.sam"
    output:
        "results/alignments/minimap2/SRR8528336_ref-at.bam"
    shell:
        "samtools view -bS {input} > {output}"

rule sort_fBAM:
    input:
        "results/alignments/minimap2/SRR8528336_ref-at.bam"
    output:
        "results/alignments/minimap2/SRR8528336_ref-at_sort.bam"
    shell:
        "samtools sort -m5G {input} -o {output}"

rule convert_to_fpileup:
    input:
        "results/alignments/minimap2/SRR8528336_ref-at_sort.bam"
    output:
        "results/alignments/minimap2/SRR8528336_ref-at_sort.pileup"
    shell:
        "samtools mpileup -B {input} > {output}"

rule SNP_calling:
    input:
        "results/alignments/minimap2/SRR8528336_ref-at_sort.pileup"
    output:
        "results/alignments/minimap2/SRR8528336_ref-at_sort.var"
    shell:
        "varscan pileup2cns {input} "
        "--min-freq-for-hom 0.6 "
        "--min-coverage 5 "
        "--min-var-freq 0.6 "
        "--p-value 0.1 "
        "--min-reads2 5 "
        "> {output}"


# extract exons from sequence genomes
rule blat_sequence_genomes:
    input:
        "data/reference_genomes/ref-at.fasta",
        "data/sequence_genomes/arl_ref.fa"
    output:
        "results/blat/at_arl.psl"
    shell:
         "blat "
         "-t=dnax "
         "-q=dnax "
         "-stepSize=5 "
         "-repMatch=2253 "
         "-minScore=20 "
         "-minIdentity=0 "
         "{input} {output}"
