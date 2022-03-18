cores=4

matches=$1

matches=`echo ${matches} | sed -r 's/\,/ /g'`


pwd=`pwd`

for match in `echo ${matches}`
do
    for f in `find . | grep -v superenhancer | grep -v hubs | grep -v rose | grep -P '(rmdup|merged).bam$' | grep ${match}`
    do
        dir=`dirname ${f}`

        echo ${dir} ${f}	

        cd ${dir}

        cmd=(qsub -N bigwig -v bam=${f},cores=${cores} -pe smp ${cores} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_bam2bigwig.sh)
        echo "${cmd[@]}"
        eval ${cmd[@]}
        
        cd ${pwd}

        #break
    done
done
