#!/bin/bash
#$ -l mem=2G,time=1::
#$ -cwd
#$ -S /bin/bash

# Extract reads for each chr since chipseeqer requires them to be
# separated

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)


pwd=`pwd`
      
for bam in `find . | grep genome | grep ${genome} | grep -v bai | grep -v rose | sort`
do
    file=`basename ${bam}`
    dir=`dirname ${bam}`

    echo ${bam}

    cd ${dir}

    # Get rid of job files
    rm submit_extract_chr_from_bam.sh.e*
    rm submit_extract_chr_from_bam.sh.o*

    rm reads.chr*

    for chr in ${CHROMS[*]}
    do
        qsub -N extract_${chr} -v file=${file},chr=${chr} /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/qsub_chr_reads.sh
    done

    cd ${pwd}
done

# Once the extraction is done, run chipseeqer
qsub -N merge_all -hold_jid "extract_chr*" -v chip_dir=${d},input_dir=${input_dir},read_length=${read_length},p=${p},genome=${genome} /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/qsub_chipseeqer.sh
