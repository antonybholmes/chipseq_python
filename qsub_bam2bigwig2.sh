#!/bin/bash
#$ -l mem=4G,time=1::
#$ -cwd
#$ -S /bin/bash

# Size from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html

bw=`echo ${bam} |sed -r 's/\.bam/.bw/'`

cmd=(bamCoverage 
    --bam a.bam 
    -o ${bw}
    --binSize 10
    --normalizeUsing RPGC
    --effectiveGenomeSize 2864785220
    --ignoreForNormalization chrX chrM
    --extendReads
    --centerReads
    --numberOfProcessors ${cores})
