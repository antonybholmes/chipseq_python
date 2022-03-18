#!/bin/bash
#$ -l mem=1G,time=:10:
#$ -cwd
#$ -S /bin/bash

export BEDSORT=/ifs/scratch/cancer/Lab_RDF/ngs/tools/ucsc/bedSort
export BIGBED=/ifs/scratch/cancer/Lab_RDF/ngs/tools/ucsc/bedToBigBed


# make a temp dir to run rose in

echo ======== Start ========
date

base=`basename ${bed}`

# for bed in *.bed
# do
    bb=`echo ${base} | sed -r 's/bed/bb/'`
    chr_size=/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/${genome}.chrom.sizes

    # remove header and keep just the chr,start, end
    cat ${base} | grep -P '^chr' | cut -f1,2,3 > ${base}.tmp

    # sort on hg18
    cmd="${BEDSORT} ${base}.tmp ${base}.tmp.sorted"
    echo ${cmd}
    ${cmd}

    # convert to big bed
    cmd="${BIGBED} ${base}.tmp.sorted ${chr_size} ${bb}"
    echo ${cmd}
    ${cmd}

    rm ${base}.tmp
    rm ${base}.tmp.sorted
# done

date
echo ======== Finished ========

