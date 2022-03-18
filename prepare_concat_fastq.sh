#!/bin/bash
#
#
# Concat fastq files
#


# a pattern
matches=$1


pwd=`pwd`

c=0

matches=`echo ${matches} | sed -r 's/\,/ /g'`

for match in `echo ${matches}`
do
    for d in `find . -maxdepth 1 -type d | grep -P ${match}`
    do
        cd ${d}
        
        echo ${d}
        cmd=(qsub -N concat /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_concat_fastq.sh)
        
        echo "${cmd[@]}"
        eval ${cmd[@]}
        
        cd ${pwd}
    done
done	
	
