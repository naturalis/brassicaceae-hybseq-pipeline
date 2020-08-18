dir=results/A06_identified_contigs_blat/
file=/stats/identified_contigs.txt

for d in $dir*;
do
    echo "$d $(cat $d$file)"
done
