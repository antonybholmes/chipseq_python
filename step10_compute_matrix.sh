cores=4

matches=$1
p=$2

matches=`echo ${matches} | sed -r 's/\,/ /g'`


pwd=`pwd`

for match in `echo ${matches}`
do
    bw=`find ${pwd} | grep -v superenhancer | grep -v hubs | grep -v rose | grep -P '.bw$' | grep ${match} | head -1`
    regions=`find ${pwd} | grep chipseeqer | grep -v superenhancer | grep -v hubs | grep -v rose | grep -P '.bed$' | grep -P "p${p}" | grep ${match} | head -1`
    
    dir=`dirname ${bw}`

    echo ${dir} ${f}	

    cd ${dir}

    cmd=(qsub -N compmat -v bw=${bw},regions=${regions},cores=${cores} -pe smp ${cores} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_compute_matrix.sh)
    echo "${cmd[@]}"
    eval ${cmd[@]}
        
    cd ${pwd}
done
