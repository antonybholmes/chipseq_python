#!/bin/bash
#$ -l mem=8G,time=24::
#$ -cwd
#$ -S /bin/bash

# export PATH=/ifs/scratch/cancer/Lab_RDF/abh2138/miniconda2/bin:${PATH}
#export PATH=/ifs/scratch/cancer/Lab_RDF/ngs/apps/python/miniconda2/bin:${PATH}

conda activate python27

export PATH=${PATH}:/ifs/scratch/cancer/Lab_RDF/ngs/tools/weblogo/weblogo
export PATH=${PATH}:/ifs/scratch/cancer/Lab_RDF/ngs/tools/blat
export PATH=${PATH}:/ifs/scratch/cancer/Lab_RDF/ngs/tools/homer/bin
export PATH=${PATH}:/ifs/scratch/cancer/Lab_RDF/ngs/tools/samtools-1.8/bin #0.1.19
export PATH=${PATH}:/nfs/apps/R/3.1.2/:/nfs/apps/R/3.1.2/bin

export SAMPLE_DIR=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/data/samples
export ROSE_PATH=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/tools/ROSE
export PATH=${PATH}:${ROSE_PATH}
export PYTHONPATH=${PYTHONPATH}:${ROSE_PATH}

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/nfs/apps/boost/1.5.4/lib/

# make a temp dir to run rose in

echo
echo "==== Start "`date`" ===="
echo
echo "dir: `pwd`"
echo "os: `cat /etc/redhat-release`"
echo

out=rose_s${stitch}_t${tss}_p${p}

mkdir -p ${out}

cmd=(python ROSE_main.py -g ${genome} -i ${gff_name} -r ${chip_name} -o ${out} -c ${input_name} -s ${stitch} -t ${tss})
echo "${cmd[@]}"
eval ${cmd[@]}

echo "==== End "`date`" ===="
echo

