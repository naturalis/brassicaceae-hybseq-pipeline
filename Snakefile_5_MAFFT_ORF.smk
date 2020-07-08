# variables for every species
SAMPLES = "SRR8528338 SRR8528339 SRR8528340".split()
SAMPLE_NAME = "SRR8528336"

configfile: "./envs/exons.yaml"

exon_variables = expand("results/11_aligned_exons_ORF/{exon}.fasta", exon=config["exons"])

rule all:
    input:
        exon_variables

rule MAFFT_assembly:
    input:
        "results/10_all_samples_exons/{exon}.fasta"
    output:
        "results/11_aligned_exons_ORF/{exon}.fasta"
    shell:
        "mafft "
        "--maxiterate 1000 "
        "--oldgenafpair "
        "{input} > {output}"

rule trim_alignments:
    input:
         "src/trim_alignments.py"
    output:
          "results/12_trimmed_alignments/"
    shell:
         "python {input}"


