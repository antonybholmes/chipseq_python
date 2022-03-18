# Extract reads for each chr since chipseeqer requires them to be
# separated

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)
#CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

BATCH_SIZE=100

# a pattern
matches=$1
genome=$2
mode=$3

if [[ -z ${mode} ]]
then
  mode=single
fi

if [[ ${mode} == "paired.v2" ]]
then
  pattern=markdup.bam
elif [[ ${mode} == "v1" ]]
then
  pattern=sorted.rmdup.bam
else
  pattern=merged.bam
fi

pwd=`pwd`

c=0

matches=`echo ${matches} | sed -r 's/\,/ /g'`


for match in `echo ${matches}`
do
    echo ${match} `pwd`
    
    for d in `find . -maxdepth 2 -type d | grep -P ${match} | sort`
    do
        echo ${d}
        
        # if target files have not been made, lets make some read
        # files
        
        cd ${d} #/reads
        
        cwd=`pwd`
      
        for bam in `find . | grep genome | grep ${genome} | grep -P "${pattern}$" | grep -v bai | grep -v rose | sort`
        do
            file=`basename ${bam}`
            dir=`dirname ${bam}`
        
            echo ${bam}
        
            cd ${dir}
        
            # Get rid of job files
            # rm submit_extract_chr_from_bam.sh.e*
            # rm submit_extract_chr_from_bam.sh.o*

            # rm reads.chr*

            for chr in ${CHROMS[*]}
            do
                #echo ${chr}
                qsub -N ex_${chr} -v file=${file},chr=${chr} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_chr_reads.sh
            done
        
            cd ${cwd}
        done
      
        cd ${pwd}
        
        c=$((c+1))
        
        if [[ ${c} -eq ${BATCH_SIZE} ]]
        then
            break
        fi
    done
done
