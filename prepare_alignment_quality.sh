#
# Run multiple alignments based on a directory pattern match
#

# a pattern
matches=$1
genome=$2

pwd=`pwd`

matches=`echo ${matches} | sed -r 's/\,/ /g'`

for match in `echo ${matches}`
do
    for d in `find . -maxdepth 1 -type d | grep -P ${match}`
    do
        if [[ ! -d ${d} ]]
        then
            continue
        fi
        
        # Do not process if the reads file cannot be found
        #if [[ `ls ${d}/*.fastq.gz | wc -l` -eq 0 ]]
        #then
        #	continue
        #fi
        
        cd ${d}
        
        #fastq=`ls *.fastq.gz` #`ls *.fastq.gz | head -1`
        #fastq=`echo ${fastq} | sed -r 's/\s/,/g'`
      
        # rm submit_align_bowtie2_gz.sh.*
        
        cmd=(qsub -v match=${match},genome=${genome} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_alignment_quality.sh)
        echo "${d} ${cmd[@]}"
        eval ${cmd[@]}
            
        cd ${pwd}
    done
done
	
	
