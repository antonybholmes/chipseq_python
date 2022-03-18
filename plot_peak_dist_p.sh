# make the tss plots
f1=$1
n1=$2
c1=$3
f2=$4
n2=$5
c2=$6
out=$7

d=`pwd`

echo $out

Rscript /ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/plot_peak_dist_p.r ${f1} ${n1} ${c1} ${f2} ${n2} ${c2} ${out}
