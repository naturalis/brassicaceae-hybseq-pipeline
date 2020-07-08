#!/bin/bash

# Script that codes the pipeline for SNP calling with SAMtools/VarScan
# Example: makeconsensus.sh <sam filename> 
# Requires set PICARD_PATH in the executable path. Also, make sure you have samtools installed and added to your PATH variable.

export OMP_NUM_THREADS=1
thread=1

# Options for variant calling (same as the Zingiberales people except the minimum coverage that I cranked up to 5):
# minimum coverage at a position
mincoverage=5
# minimum frequency of observed allele
minvarfreq=0.6
# p-value
pvalue=0.1
# minimum number of reads supporting position
minreads=5

samfile=$(basename "$1")
pref="${samfile%.*}"
outdir=`pwd`

echo "The sam file is: $samfile"
echo "Output will be on: $outdir"

# sort the SAM file 
echo "Sort the SAM file and make BAM file ..."
samtools view -bS ${samfile} > ${pref}.bam
samtools sort -m5G ${pref}.bam -o ${pref}_sort.bam 

# make a pileup file
echo "Make a pileup file ..."
samtools mpileup -B ${pref}_sort.bam > ${pref}.pileup

# Variant SNPs and INDEL calling && filtering
java -Xmx5G -jar $VarScan_PATH/VarScan.v2.4.1.jar pileup2cns ${pref}.pileup --min-freq-for-hom $minvarfreq --min-coverage $mincoverage --min-var-freq $minvarfreq --p-value $pvalue --min-reads2 $minreads > ${pref}.var

# clean the *.bam files
rm ${pref}.bam ${pref}_sort.bam ${pref}.pileup
