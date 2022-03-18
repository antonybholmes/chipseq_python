#!/bin/bash -l
#$ -l mem=1G,time=1::
#$ -cwd
#$ -S /bin/bash

echo ======== Start ========
date

export SCRIPT_DIR=/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts
export HISAT2_INDEXES=/ifs/scratch/cancer/Lab_RDF/abh2138/references/hisat2/${genome}/
export HISAT2=/ifs/scratch/cancer/Lab_RDF/abh2138/tools/hisat2-2.0.5/hisat2
export SAMTOOLS=/ifs/scratch/cancer/Lab_RDF/abh2138/tools/samtools-1.8/bin/samtools

pwd=`pwd`
lib=`basename ${pwd}`

echo Library ${lib}


echo Alignment quality...
# create alignment quality file
. ${SCRIPT_DIR}/alignment_quality.sh ${lib} ${genome} hisat2

date
echo ======== Finished ========
