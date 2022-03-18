POWERS="100 1000 10000 100000 1000000 10000000" # 100000000" # 1000000000"

matches=$1
mode=count

pwd=`pwd`

matches=`echo ${matches} | sed -r 's/\,/ /g'`

for match in `echo ${matches}`
do
    echo Match on ${match}
    
    for f in `find . | grep -P '(rmdup|merged).bam$' | grep -v rose | grep -v superenhancers | grep ${match}`
    do
        echo Found bam ${f}
        
        dir=`dirname ${f}`
        file=`basename ${f}`
      
        cd ${dir}
      
        for power in `echo ${POWERS}`
        do
            cmd=(qsub -N bincount -v bam=${file},power=${power},mode=${mode} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_bincounts.sh)
            echo "${cmd[@]}"
            eval ${cmd[@]}
            
            #break
        done
      
        cd ${pwd}
      
        #break
    done
done
