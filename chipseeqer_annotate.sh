#!/bin/bash
#$ -l mem=8G,time=2::
#$ -cwd
#$ -S /bin/bash

export CHIPSEEQERDIR=/ifs/scratch/cancer/Lab_RDF/abh2138/chip-seq/ChIPseeqer/dist
export PATH=$PATH:$CHIPSEEQERDIR:$CHIPSEEQERDIR/SCRIPTS
export PERL5LIB=$PERL5LIB:$CHIPSEEQERDIR:$CHIPSEEQERDIR/SCRIPTS

echo start
date

for file in `ls | grep -P 'TF_targets.+txt$' | grep -P -v 'refSeq'`
do
	echo ${file}
	echo ChIPseeqerAnnotate --genome=hg18 --peakfile=${file}
	ChIPseeqerAnnotate --genome=hg18 --peakfile=${file}
done

echo end
date
