#!/bin/bash -l
#$ -l mem=8G,time=2::
#$ -cwd
#$ -S /bin/bash

# merge reads on a per chromosome basis



# export SAMTOOLS=/ifs/scratch/cancer/Lab_RDF/ngs/tools/samtools-1.8/bin/samtools
export SCRIPT_DIR=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts

echo
echo "==== Start "`date`" ===="
echo
echo "dir: `pwd`"
echo "os: `cat /etc/redhat-release`"
echo


dir=`dirname ${bam}`

if [[ ! -f ${bam} ]]
then
    tmp=${bam}.tmp.1
    tmp2=${bam}.tmp.2

    ls ${dir}/*merged*tmp | sort -n -k1 -t'.' > ${dir}/merge_files.txt

    echo "Merging ${bam}..."
    cmd=(samtools merge -f -b ${dir}/merge_files.txt ${tmp}) #*merged.chr.*tmp
    echo "${cmd[@]}"
    eval ${cmd[@]}

    # Re-sort by position otherwise indexing can fail
    echo "Sorting ${tmp} into ${bam}..."
    cmd=(samtools sort -@ ${cores} -o ${tmp2} ${tmp})
    echo "${cmd[@]}"
    eval ${cmd[@]}

    rm ${tmp}
    mv ${tmp2} ${bam}
    
    # Remove temp files
    # rm ${tmp}
    rm ${dir}/*merged*tmp
fi

if [[ -f ${bam} && ! -f ${bam}.bai ]]
then
    echo "Indexing ${bam}..."
    samtools index ${bam}
fi




# Remove markdup as merge is better
# find ${dir} | grep markdup.bam | xargs rm


${SCRIPT_DIR}/alignment_quality.sh ${lib} ${genome} hisat2


date
echo "======== Finished ========"
