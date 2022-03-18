# merge reads on a per chromosome basis

#!/bin/bash -l
#$ -l mem=8G,time=1::
#$ -cwd
#$ -S /bin/bash


export DIR=/ifs/scratch/cancer/Lab_RDF/ngs
export SCRIPT_DIR=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts
export HISAT2_INDEXES=/ifs/scratch/cancer/Lab_RDF/ngs/references/hisat2/${genome}/
export HISAT2=/ifs/scratch/cancer/Lab_RDF/ngs/tools/hisat2-2.0.5/hisat2
export SAMTOOLS=/ifs/scratch/cancer/Lab_RDF/ngs/tools/samtools-1.8/bin/samtools

echo ======== Start ========
date

dir=`dirname ${bam}`

oid=`printf "%02d" ${id}`

# create a temp output file for this chr
tmp=${dir}/${oid}.merged.${chr}.tmp #`mktemp ${bam}.merge.${chr}.tmp.XXXXXXXX`

cmd=(python 
    /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/python/merge_pairs/merge_pairs_chr.py
    --samtools=${SAMTOOLS}
    --out=${tmp}
    ${bam}
    ${chr})
    
echo "${cmd[@]}"
eval ${cmd[@]}

date
echo ======== Finished ========
