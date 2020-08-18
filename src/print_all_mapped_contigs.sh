dir=results/A03_mapped_reads/
file=/contigs_formed.txt

for d in $dir*;
do
    echo "$d $(cat $d$file)"
done
