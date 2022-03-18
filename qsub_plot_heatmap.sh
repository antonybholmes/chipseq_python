#!/bin/bash
#$ -l mem=8G,time=2::
#$ -cwd
#$ -S /bin/bash

# Size from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html

echo "======== Start ========"
date

${matrix}=`basename ${matrix}`

out=${matrix}.png


cmd=(plotHeatmap -m ${matrix} -out ${out} --dpi 600)
    
echo "${cmd[@]}"
eval ${cmd[@]}

date
echo "======== Finished ========"

