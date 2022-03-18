matches=$1

matches=`echo ${matches} | sed -r 's/\,/ /g'`

pwd=`pwd`

for match in `echo ${matches}`
do
    echo ${match}
    
    for d in `find ${pwd} -mindepth 1 -maxdepth 1 -type d | grep ${match}`
    do
        echo ${d}
          
        cd ${d}

        if [ ! -d "fastqc" ]
        then
          mkdir fastqc
        fi

        for f in `find fastq | grep -P 'fastq.gz$'`
        do
            #f2=`echo ${f} | sed -r 's/R1/R2/'`
            echo ${f}
            dir=`dir ${f}`
            name=`basename ${dir}`
          

            cmd=(qsub -N fastqc -v name=${name},f1=${f} /ifs/scratch/cancer/Lab_RDF/ngs/rna_seq/scripts/qsub_fastqc.sh)
            echo "${cmd[@]}"
            eval ${cmd[@]}
          
            #break
        done
        
        cd ${pwd}
    done
done
