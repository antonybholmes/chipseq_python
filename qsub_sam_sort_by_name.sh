#!/bin/bash -l
#$ -l mem=2G,time=1::
#$ -pe smp 8
#$ -cwd
#$ -S /bin/bash


export DIR=/ifs/scratch/cancer/Lab_RDF/abh2138
export SCRIPT_DIR=/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts
export HISAT2_INDEXES=/ifs/scratch/cancer/Lab_RDF/abh2138/references/hisat2/${genome}/
export HISAT2=/ifs/scratch/cancer/Lab_RDF/abh2138/tools/hisat2-2.0.5/hisat2
export SAMTOOLS=/ifs/scratch/cancer/Lab_RDF/abh2138/tools/samtools-1.8/bin/samtools


echo ======== Start ========
date

out_bam=`echo ${bam} | sed -r 's/sorted\.//'`
cmd="${SAMTOOLS} sort -@ 8 -n -o ${out_bam} ${bam}"
echo ${cmd}
${cmd}

date
echo ======== Finished ========
