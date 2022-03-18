#!/bin/bash -l
#$ -l mem=4G,time=4::
#$ -cwd
#$ -S /bin/bash

echo ======== Start ========
date

out=${bam}.stats.txt

echo ${bam} ${out}
samtools stats ${bam} > ${out}

date
echo ======== Finished ========
