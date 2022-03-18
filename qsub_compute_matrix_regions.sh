#!/bin/bash
#$ -l mem=8G,time=1::
#$ -cwd
#$ -S /bin/bash

# Size from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html

echo "======== Start ========"
date

f=`basename ${bw}`

matrix=${f}.${name}.matrix.gz
tsv=${f}.${name}.regions.tsv
bed=${f}.${name}.sorted.bed

bp=5000
binsize=20 # 50 bins

cmd=(computeMatrix 
    scale-regions
    --regionsFileName ${regions}
    --scoreFileName ${bw}
    --outFileName ${matrix}
    --outFileNameMatrix ${tsv}
    --outFileSortedRegions ${bed}
    --numberOfProcessors ${cores}
    --binSize ${binsize}
    --verbose)
    
echo "${cmd[@]}"
eval ${cmd[@]}


#out=${matrix}.png


#cmd=(plotHeatmap 
    #--matrixFile ${matrix}
    #--outFileSortedRegions ${bed}
    #--outFileName ${out} 
    #--colorMap Reds 
    #--yMax 20 
    #--zMax 20 
    #--dpi 300)
    
#echo "${cmd[@]}"
#eval ${cmd[@]}


date
echo "======== Finished ========"

