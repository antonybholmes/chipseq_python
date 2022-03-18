matches=$1
genome=hg19

matches=`echo ${matches} | sed -r 's/\,/ /g'`


pwd=`pwd`

for match in `echo ${matches}`
do
    for f in `find . | grep -v superenhancer | grep -v hubs | grep -P 'bed$' | grep ${match}`
    do
        dir=`dirname ${f}`

        echo ${dir} ${f}	

        cd ${dir}

        qsub -v bed=${f},genome=${genome} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_bigbed.sh

        cd ${pwd}

        #break
    done
done
