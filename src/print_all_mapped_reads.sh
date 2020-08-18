dir=results/A03_mapped_reads/                                                                                       file=/contigs_formed.txt
file=/mapped_reads.txt

for d in $dir*;
do
    echo "$d $(cat $d$file)"
done
