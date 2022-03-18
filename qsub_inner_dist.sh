#!/bin/bash -l
#$ -l mem=4G,time=2::
#$ -cwd
#$ -S /bin/bash

echo ======== Start ========
date

/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/inner_dist.sh ${bam}

date
echo ======== Finished ========
