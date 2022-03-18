#!/bin/bash
#$ -l mem=8G,time=2::
#$ -cwd
#$ -S /bin/bash

#export PATH=/ifs/scratch/cancer/Lab_RDF/ngs/apps/python/miniconda3/bin:$PATH

#base=`basename ${bam}`
#bw=`echo ${base} | sed -r 's/\.bam/\.bw/'`
#echo ${base} ${bw}
#bamCoverage --bam ${base} -o ${bw} -of bigwig -p4



# Size from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html

echo "======== Start ========"
date

bam=`basename ${bam}`

bw=${bam}.bw

cmd=(bamCoverage 
    --bam ${bam}
    --outFileFormat bigwig
    --outFileName ${bw}
    --binSize 10
    --normalizeUsing RPGC
    --effectiveGenomeSize 2864785220
    --ignoreForNormalization chrX chrM
    --centerReads
    --numberOfProcessors ${cores})
    
echo "${cmd[@]}"
eval ${cmd[@]}

date
echo "======== Finished ========"

