#!/bin/bash

# Script that codes the pipeline for short-read assembly using Trimmomatic/Fastx_collapse/Alignreads 
# Example: runhmap <reference> <read1> <read2> <file prefix> 
# Requires set TRIM_PATH, and fastx-tools and alignreads in the executable path.

# Common path to the reads database
rdata_path=/Users/lukenikolov/Desktop/raw_data/bastet.ccg.uni-koeln.de/downloads/lnikolov 

thread=2
export OMP_NUM_THREADS=$thread
ref=$1
read1=$2
read2=$3
pref=$4
outdir=$(pwd)

# non-default options for alignreads
depthmasking=5
propbfilter=0.7
percentidentity=medium

echo "The reference is: $ref"
echo "Reads are stored in: $read1 $read2"
echo "The files prefix is: $pref"
echo "We work on: $thread threads"
echo "Output will be on: $outdir"

# trim the reads
echo "Trimming reads ..."
java -jar $TRIM_PATH/trimmomatic-0.36.jar PE -phred33 "$rdata_path/$read1" "$rdata_path/$read2" "${pref}_forward_trim_paired.fq" "${pref}_forward_trim_unpaired.fq" "${pref}_rev_trim_paired.fq" "${pref}_rev_trim_unpaired.fq" LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36

# remove duplicate reads 
fastx_collapser -v -i "${pref}_forward_trim_paired.fq" -o "${pref}_forward_trim_paired_dedupl.fq"
fastx_collapser -v -i "${pref}_rev_trim_paired.fq" -o "${pref}_rev_trim_paired_dedupl.fq" 
fastx_collapser -v -i "${pref}_forward_trim_unpaired.fq" -o "${pref}_forward_trim_unpaired_dedupl.fq"
fastx_collapser -v -i "${pref}_rev_trim_unpaired.fq" -o "${pref}_rev_trim_unpaired_dedupl.fq"

# remove intemediate files
rm "${pref}_forward_trim_paired.fq" "${pref}_rev_trim_paired.fq" "${pref}_forward_trim_unpaired.fq" "${pref}_rev_trim_unpaired.fq"  

# combine all reads in one file
cat "${pref}_forward_trim_paired_dedupl.fq" "${pref}_rev_trim_paired_dedupl.fq" "${pref}_forward_trim_unpaired_dedupl.fq" "${pref}_rev_trim_unpaired_dedupl.fq" > "${pref}_reads.fq"

# remove separate files
rm "${pref}_forward_trim_paired_dedupl.fq" "${pref}_rev_trim_paired_dedupl.fq" "${pref}_forward_trim_unpaired_dedupl.fq" "${pref}_rev_trim_unpaired_dedupl.fq"

# align using alignreads
echo "Alignment ..."
echo "Using alignreads"
echo "Non-default parameters:"
echo "Percent identity=${percentidentity}"
echo "Depth-masking=${depthmasking}-"
echo "Proportion-base-filter=${propbfilter}-"

alignreads "${pref}_reads.fq" "$ref" --single-step --read-type solexa --read-orientation linear --percent-identity ${percentidentity} --depth-position-masking "${depthmasking}-" --proportion-base-filter "${propbfilter}-" 

# remove the trimmed and de-duplicated reads
rm "${pref}_reads.fq"

