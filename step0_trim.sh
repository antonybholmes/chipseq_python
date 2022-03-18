matches=$1
matches=`echo ${matches} | sed -r 's/\,/ /g'`

cores=4

pwd=`pwd`

for match in `echo ${matches}`
do
    echo ${match}
    
    for d in `find ${pwd} -mindepth 1 -maxdepth 1 -type d | grep ${match}`
    do
        echo ${d}/fastq
        cd ${d}/fastq

		f1=`find . | grep ${match} | grep -P 'fastq.gz$' | grep merge | grep -v fq | grep -v trash | grep -v trim | grep -P '_1'`

		if [ ! -f "${f1}" ]
		then
			echo Looking for non merged
			f1=`find . | grep ${match} | grep -P 'fastq.gz$' | grep -v fq | grep -v trash | grep -v trim | grep -P '_1'`
		fi
		
		if [ ! -f "${f1}" ]
		then
			echo Looking for non merged R1/R2
			f1=`find . | grep ${match} | grep -P 'fastq.gz$' | grep -v fq | grep -v trash | grep -v trim | grep -P 'R1'`
		fi
		
		
		echo Found ${f1}...
		
		# -l centos7=1
		
		cmd=(qsub -N trim_g -v fastq1=${f1},cores=${cores} -pe smp ${cores} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_trim.sh)
		echo "${cmd[@]}"
		eval ${cmd[@]}
		
		cd ${pwd}
    done
done
