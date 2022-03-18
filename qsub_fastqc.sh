#!/bin/bash
#$ -l mem=4G,time=:30:
#$ -cwd
#$ -S /bin/bash

# export PATH=$PATH:/nfs/apps/java/1.8.0_91/bin/
# export JAVA_HOME=$JAVA_HOME:/nfs/apps/java/1.8.0_91/bin/

# module load java

echo
echo "==== Start "`date`" ===="
echo
echo "dir: `pwd`"
echo "os: `cat /etc/redhat-release`"
echo

sample=`basename ${name} | sed -r 's/\..+//'`

d=fastqc/${sample}

mkdir -p ${d}

/ifs/scratch/cancer/Lab_RDF/ngs/tools/FastQC/fastqc --quiet --outdir ${d} ${f1}

echo "==== End "`date`" ===="
echo

