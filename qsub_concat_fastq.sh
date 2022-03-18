# merge reads on a per chromosome basis

#!/bin/bash -l
#$ -l mem=1G,time=1::
#$ -cwd
#$ -S /bin/bash

echo ======== Start ========
date


cat `ls *_1.fastq.gz | grep -v concat | sort` > concat_1.fastq.gz
cat `ls *_2.fastq.gz | grep -v concat | sort` > concat_2.fastq.gz

date
echo ======== Finished ========
