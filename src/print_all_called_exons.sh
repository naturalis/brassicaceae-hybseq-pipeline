dir=results/A09_consensus_exons/
file=/called_exons.txt

for d in $dir*;
do
    echo "$d $(cat $d$file)"
done
