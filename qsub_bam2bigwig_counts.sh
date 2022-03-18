#!/bin/bash
#$ -l mem=8G,time=2::
#$ -cwd
#$ -S /bin/bash

export PATH=/ifs/scratch/cancer/Lab_RDF/ngs/apps/python/miniconda3/bin:$PATH


# Size from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html

echo "======== Start ========"
date

bam=`basename ${bam}`

bw=${bam}.counts.bw

cmd=(bamCoverage 
    --bam ${bam}
    --outFileFormat bigwig
    --outFileName ${bw}
    --binSize 10
    --centerReads
    --normalizeUsing None
    --numberOfProcessors ${cores})
    
echo "${cmd[@]}"
eval ${cmd[@]}

date
echo "======== Finished ========"

