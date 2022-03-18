#!/bin/bash
#$ -l mem=1G,time=1::
#$ -cwd
#$ -S /bin/bash

# export SAMTOOLS=/ifs/scratch/cancer/Lab_RDF/ngs/tools/samtools-1.8/bin/samtools

out=reads.${chr}

echo
echo "==== Start "`date`" ===="
echo
echo "dir: `pwd`"
echo "os: `cat /etc/redhat-release`"
echo

if [[ -f "${file}" ]]
then
	echo "Extract reads for ${chr}: ${file} -> ${out}"
	samtools view -F 4 ${file} ${chr} > ${out}
	# ${SAMTOOLS} view -b -F 4 -o ${out} ${file} ${chr}
fi

echo "==== End "`date`" ===="
echo
