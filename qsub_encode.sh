#!/bin/bash
#$ -l mem=4G,time=1::
#$ -t 1-25
#$ -cwd

#$ -S /bin/bash

# ##############$ -t 1-25

PYTHONPATH=/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/python


echo ======== Start ========
date

if [ ${genome} == "hg19" ]
then
  if [ "${SGE_TASK_ID}" == 25 ]
  then
    chr=chrM
  elif [ "${SGE_TASK_ID}" == 24 ]
  then
    chr=chrY
  elif [ "${SGE_TASK_ID}" == 23 ]
  then
    chr=chrX
  else
    chr=chr${SGE_TASK_ID}
  fi
else
  if [ "${SGE_TASK_ID}" == 22 ]
  then
    chr=chrM
  elif [ "${SGE_TASK_ID}" == 21 ]
  then
    chr=chrY
  elif [ "${SGE_TASK_ID}" == 20 ]
  then
    chr=chrX
  else
    chr=chr${SGE_TASK_ID}
  fi
fi

if [ ${genome} == "hg19" ]
then
  chr_size_file=/ifs/scratch/cancer/Lab_RDF/abh2138/references/ucsc/hg19_chromosome_sizes.txt
else
  chr_size_file=/ifs/scratch/cancer/Lab_RDF/abh2138/references/ucsc/mm10_chromosome_sizes.txt
fi

file=reads.${chr}

echo ${chr} ${out}
	
#python ${PYTHONPATH}/encode_reads_16bit.py ${chr_size_file} ${file} ${chr} ${read_length} ${window} > ${out}
python ${PYTHONPATH}/encode_reads.py ${chr_size_file} ${file} ${chr} ${read_length} ${window}


echo end
date
