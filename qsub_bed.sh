#!/bin/bash
#$ -l mem=1G,time=1::
#$ -cwd
#$ -S /bin/bash


echo
echo "==== Start "`date`" ===="
echo
echo "dir: `pwd`"
echo "os: `cat /etc/redhat-release`"
echo

/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/create_peak_bed_file.sh ${genome}

echo "==== End "`date`" ===="
echo
