#!/bin/bash
#$ -l mem=4G,time=1::
#$ -cwd
#$ -S /bin/bash

# Size from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html

bw=`echo ${bam} |sed -r 's/\.bam/.cpm.bw/'`

cmd=(bamCoverage 
    --bam a.bam 
    -o ${bw}
    --binSize 10
    --normalizeUsing CPM
    --ignoreForNormalization chrX chrM
    --extendReads
    --centerReads
    --numberOfProcessors ${cores})
