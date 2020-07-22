# variables for every species
SAMPLES = ["SRR8528336", "SRR8528337", "SRR8528338", "SRR8528339", "SRR8528340", "SRR8528341"]
SAMPLE_NAME = "SRR8528337"

configfile: "./envs/MAFFT/" + SAMPLE_NAME + ".yaml"

mafft_variables = expand("results/8_aligned_exons/"+ SAMPLE_NAME + "/{exon}.fasta", exon=config["exons"])

rule all:
    input:
        mafft_variables

rule MAFFT_assembly:
    input:
        "results/7_mapped_exons/"+ SAMPLE_NAME + "/{exon}.fasta"
    output:
        "results/8_aligned_exons/"+ SAMPLE_NAME + "/{exon}.fasta"
    shell:
        "mafft "
        "--maxiterate 1000 "
        "--oldgenafpair "
        "{input} > {output}"

rule make_exon_consensus:
    input:
        "src/make_exon_consensus.py"
    shell:
        "python {input} " + SAMPLE_NAME

rule merge_exon_seqs:
    input:
         "src/merge_exon_seqs.py"
    output:
        "results/10_all_samples_exons"
    shell:
         "python {input}"
