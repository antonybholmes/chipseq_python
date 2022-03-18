#!/bin/bash -l
#$ -l mem=8G,time=1::
#$ -cwd
#$ -S /bin/bash

echo ======== Start ========
date

echo ${bam}

python3 /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/python/bincounts.py \
${bam} \
hg19 \
${power} \
${mode}

date
echo ======== Finished ========
