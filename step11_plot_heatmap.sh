cores=4

matches=$1
p=$2

matches=`echo ${matches} | sed -r 's/\,/ /g'`


pwd=`pwd`

for match in `echo ${matches}`
do
    matrix=`find ${pwd} | grep -v superenhancer | grep -v hubs | grep -v rose | grep -P '.matrix.gz$' | grep ${match} | head -1`

    dir=`dirname ${matrix}`

    echo ${dir} ${matrix}	

    cd ${dir}

    cmd=(qsub -N heatmap -v matrix=${matrix},cores=${cores} -pe smp ${cores} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_plot_heatmap.sh)
    echo "${cmd[@]}"
    eval ${cmd[@]}
        
    cd ${pwd}
done
