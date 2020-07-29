# variables for every species
SAMPLES = "SRR8528338 SRR8528339 SRR8528340".split()

configfile: "./envs/macse.yaml"

macse_AA_variables = expand("results/13_prot_alignments/{exon}_AA.fasta", exon=config["exons_AA"])
masce_NT_variables = expand("results/13_prot_alignments/{exon}_NT.fasta", exon=config["exons_NT"])


rule all:
    input:
        macse_AA_variables,
        masce_NT_variables

rule macse:
    input:
         macse = "src/macse_v1.2.jar",
         fasta = "results/12_trimmed_alignments/{exon}.fasta"
    output:
         NT = "results/13_prot_alignments/{exon}_NT.fasta",
         AA = "results/13_prot_alignments/{exon}_AA.fasta"
    shell:
         "java -jar {input.macse} -prog alignSequences -seq {input.fasta} -out_NT {output.NT} -out_AA {output.AA}"
