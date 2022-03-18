#!/bin/bash
#$ -l mem=8G,time=1::
#$ -cwd
#$ -S /bin/bash

# Size from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html

echo "======== Start ========"
date

f=`basename ${bw}`

matrix=${f}.${name}.matrix.gz
tsv=${f}.${name}.matrix.tsv
bed=${f}.${name}.sorted.bed

bp=5000

cmd=(computeMatrix 
    reference-point
    --referencePoint center
    --upstream ${bp}
    --downstream ${bp}
    --binSize 100
    --regionsFileName ${regions}
    --scoreFileName ${bw}
    --outFileName ${matrix}
    --outFileNameMatrix ${tsv}
    --numberOfProcessors ${cores}
    --verbose)
    
echo "${cmd[@]}"
eval ${cmd[@]}


out=${matrix}.png


cmd=(plotHeatmap 
    --matrixFile ${matrix}
    --outFileName ${out} 
    --colorMap Reds 
    --sortRegions no 
    --yMax 20 
    --zMax 20 
    --dpi 300)
    
echo "${cmd[@]}"
eval ${cmd[@]}


date
echo "======== Finished ========"

