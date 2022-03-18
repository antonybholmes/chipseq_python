#!/bin/bash
#
#
# Run multiple alignments based on a directory pattern match
#

BATCH_SIZE=30

# a pattern
matches=$1

# len=36 #100

# Two mismatches use -2
mismatches=$2 #-10 #-5

genome=$3

mode=$4

cores=$5

if [[ -z ${mode} ]]
then
  mode=single
fi

if [[ -z ${cores} ]]
then
  cores=8 #4
fi

pwd=`pwd`

c=0

matches=`echo ${matches} | sed -r 's/\,/ /g'`

for match in `echo ${matches}`
do
    for d in `find . -maxdepth 1 -type d | grep -P ${match}`
    do
        # Do not process if the bam file already exists
        #if [[ `find ${d} | grep bam | wc -l` -gt 0 ]]
        #then
        #	continue
        #fi
        
        #if [[ ${mode} == "paired" ]]
        #then
            #if [[ `find ${d} | grep merged.bam | wc -l` -gt 0 ]]
            #then
                # If merged already exists, skip processing
            #    continue
            #fi
        #fi
        
        cd ${d}
        
        # Old bowtie method
        #qsub -v match=${match},mismatches=${mismatches},genome=${genome} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/submit_align_bowtie2.sh
        
        echo ${d}
        cmd=(qsub -N hisat2 -v match=${match},mismatches=${mismatches},genome=${genome},mode=${mode},cores=${cores} -pe smp ${cores} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_hisat2.sh)
        
        echo "${cmd[@]}"
        eval ${cmd[@]}
        
        cd ${pwd}
        
        c=$((c+1))
        
        if [[ ${c} -eq ${BATCH_SIZE} ]]
        then
            break
        fi
    done
done	
	
