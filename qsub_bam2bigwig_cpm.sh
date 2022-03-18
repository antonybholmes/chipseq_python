#!/bin/bash
#$ -l mem=8G,time=2::
#$ -cwd
#$ -S /bin/bash

echo "======== Start ========"
date

bam=`basename ${bam}`

bw=${bam}.cpm.bw

cmd=(bamCoverage 
    --bam ${bam}
    --outFileFormat bigwig
    --outFileName ${bw}
    --binSize 10
    --normalizeUsing CPM
    --ignoreForNormalization chrX chrM
    --centerReads
    --numberOfProcessors ${cores})


echo "${cmd[@]}"
eval ${cmd[@]}

date
echo "======== Finished ========"
