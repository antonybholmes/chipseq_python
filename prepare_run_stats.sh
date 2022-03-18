
matches=$1

pwd=`pwd`

matches=`echo ${matches} | sed -r 's/\,/ /g'`

for match in `echo ${matches}`
do
    echo Match on ${match}
    
    for f in `find . | grep -P '\.bam$' | grep -v rmdup | grep -v merged | grep -v rose | grep -v superenhancers | grep ${match}`
    do
        echo Found bam ${f}
        
        dir=`dirname ${f}`
        file=`basename ${f}`
      
        cd ${dir}

	cmd="qsub -v bam=${file} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_stats.sh"
        echo ${cmd}
        ${cmd}
      
        cd ${pwd}
      
        #break
    done
done
