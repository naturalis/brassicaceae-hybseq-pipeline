SAMPLES = ["SRR8528336"]
FRPU = ["forward_trim_paired","forward_trim_unpaired","reverse_trim_paired","reverse_trim_unpaired"]
#forward or reverse | paired or unpaired

rule trimming:
    input:
        "data/raw_reads/SRR8528336_1.fastq.gz",
        "data/raw_reads/SRR8528336_2.fastq.gz"
    output:
        temp(expand("trimmed_reads/SRR8528336/SRR8528336_{FRPU}.fq", FRPU = FRPU))
    shell:
        "trimmomatic PE -phred33 {input} {output} ILLUMINACLIP:trimmomatic_adapter/TruSeq3-PE-2.fa:2:30:10 LEADING:20 "
        "TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36"

rule count_reads_trimming:
    input:
        expand("trimmed_reads/SRR8528336/SRR8528336_{FRPU}.fq", FRPU = FRPU)
    output:
        "trimmed_reads/SRR8528336/SRR8528336_count_reads.txt"
    shell:
        "echo $(cat {input} | wc -l)/4|bc >> {output}"

rule deduplication_1:
    input:
        "trimmed_reads/SRR8528336/SRR8528336_forward_trim_paired.fq"
    output:
        "deduplicated_reads/SRR8528336/SRR8528336_forward_trim_paired_dedupl.fq"
    shell:
        "fastx_collapser -v -i {input} -o {output}"

rule deduplication_2:
    input:
        "trimmed_reads/SRR8528336/SRR8528336_forward_trim_unpaired.fq"
    output:
        "deduplicated_reads/SRR8528336/SRR8528336_forward_trim_unpaired_dedupl.fq"
    shell:
        "fastx_collapser -v -i {input} -o {output}"

rule deduplication_3:
    input:
        "trimmed_reads/SRR8528336/SRR8528336_reverse_trim_paired.fq"
    output:
        "deduplicated_reads/SRR8528336/SRR8528336_reverse_trim_paired_dedupl.fq"
    shell:
        "fastx_collapser -v -i {input} -o {output}"

rule deduplication_4:
    input:
        "trimmed_reads/SRR8528336/SRR8528336_reverse_trim_unpaired.fq"
    output:
        "deduplicated_reads/SRR8528336/SRR8528336_reverse_trim_unpaired_dedupl.fq"
    shell:
        "fastx_collapser -v -i {input} -o {output}"

rule combine:
    input:
        temp(expand("deduplicated_reads/SRR8528336/SRR8528336_{FRPU}_dedupl.fq", FRPU = FRPU))
    output:
        "deduplicated_reads/SRR8528336/SRR8528336_reads.fq"
    shell:
         "cat {input} > {output}"

rule count_reads_deduplication:
    input:
        "deduplicated_reads/SRR8528336/SRR8528336_reads.fq"
    output:
        "deduplicated_reads/SRR8528336/SRR8528336_count_reads.txt"
    shell:
        "grep '>' {input} | wc -l > {output}"
