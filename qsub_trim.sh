#!/bin/bash
#$ -l mem=1G,time=4::
#$ -cwd
#$ -S /bin/bash

trimmed=`echo ${fastq} | sed -r 's/fastq.gz/trim.noadapt.fastq/'`
gz=${trimmed}.gz
cores=1

echo
echo "==== Start "`date`" ===="
echo
echo "dir: `pwd`"
echo "os: `cat /etc/redhat-release`"
echo


#cmd="cutadapt \
#-q 20,20 \
#--cores=${cores} \
#-o ${trimmed} ${fastq} \
#-a AGATCGGAAGAG \
#--minimum-length=20"

echo "Trimming ${fastq1}..."
t1=`date -u '+%s'`

# trim 1 bp to ensure pairs don't overlap, trim 2 bp from ends in
# case still contaminated around adapter
cmd=(trim_galore --cores ${cores} --three_prime_clip_R1 2 ${fastq1})

echo "${cmd[@]}"
eval ${cmd[@]}

t2=`date -u '+%s'`
mins=$(((${t2}-${t1})/60))
echo "`date '+%b %d %H:%M:%S'` - Finished in ${mins} mins."
echo
#fi

echo "==== End "`date`" ===="
echo

