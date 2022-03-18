file=$1

out=`echo ${file} | sed -r 's/\.bed/.gff/'`

echo ${file} ${out}

cat ${file} | sed 1d | awk '{print $1"\tpeak_"NR"\tpeak\t"$2"\t"$3"\t"$4"\t+\t.\tpeak_id \"peak_"NR"\""}' > ${out}
