#export SAMTOOLS=/ifs/scratch/cancer/Lab_RDF/ngs/tools/samtools-1.8/bin/samtools

sample=$1
genome=$2
mapper=$3

bam_file=`find . | grep -P "${genome}\.bam$" | grep ${mapper} | grep -v rose | head -1`
dir=`dirname ${bam_file}`

echo "Alignment quality ${sample} ${dir} ${genome} ${bam_file}..."

#
# Count the number of reads from the fastq file
#

merged=`find . | grep -P 'merged\.bam$' | wc -l`


file1=`find . | grep -P '(fq|fastq).gz' | grep R1 | sort | tail -1`
file2=`find . | grep -P '(fq|fastq).gz' | grep R2 | sort | tail -1`

if [[ -f ${file2} ]]
then
    mode=paired
else
    mode=single
fi

echo "Mode ${mode}"

f=${dir}/reads_count_${genome}.tsv
    
if [[ ! -f "${f}" ]]
then
    if [[ ${mode} == "paired" ]]
    then
        echo "Counting paired end reads in ${file1}, writing to ${f}..."
        
        files=${file1} # ${file2}"
        reads=`zcat ${file1} | grep -P '@' | wc -l`
        #reads2=`zcat ${file2} | grep -P '@' | wc -l`
      
        #echo ${reads} > ${dir}/reads_count_1_${genome}.tsv
        #echo ${reads2} > ${dir}/reads_count_2_${genome}.tsv
    else
		echo "Counting single end reads in ${file1}, writing to ${f}..."
        reads=`zcat ${file1} | grep -P '@' | wc -l`
        #reads2="n/a"
    fi
    
    echo ${reads} > ${f}.tmp
    mv ${f}.tmp ${f}
fi

reads=`cat ${f}`

#f=${dir}/trimmed_reads_count_${genome}.tsv

#if [[ ! -f ${f} ]]
#then
    #if [[ ${mode} == "paired" ]]
    #then
        #trimmed1=`ls *.trimmed.fastq.gz | grep R1 | tail -1`
        #trimmed2=`ls *.trimmed.fastq.gz | grep R2 | tail -1`
    #else
        #trimmed1=`ls *.trimmed.fastq.gz | tail -1`
    #fi
    
    #if [[ -f ${trimmed1} ]]
    #then
        #echo "Counting trimmed reads in ${trimmed1}..."
        #trimmed_reads=`zcat ${trimmed1} | grep -P '@' | wc -l`
    #else
        #trimmed_reads=${reads}
    #fi
    
    #echo ${trimmed_reads} > ${f}
#fi

#trimmed_reads=`cat ${f}`

#
# Count reads from the BAM file
#



#reads=`samtools view -c ${bam_file}`
#echo reads ${reads}

if [[ -f ${bam_file} ]]
then
    echo "Counting mapped reads in ${bam_file}, writing to ${f}..."

    f=${dir}/mapped_reads_count_${genome}.tsv
    
    if [[ ! -f ${f} ]]
    then
        if [[ ${mode} == "paired" ]]
        then
            # just count one of the pairs
            mapped_reads=`samtools view -f 64 -c ${bam_file}`
        else
            mapped_reads=`samtools view -F 4 -c ${bam_file}`
        fi
        
        echo ${mapped_reads} > ${f}.tmp
        mv ${f}.tmp ${f}
    fi
    
    mapped_reads=`cat ${f}`
else
    sorted_file=`find . | grep -P 'sorted\.bam$' | grep ${genome} | grep ${mapper} | grep -v rose | head -1`
        
    if [[ -f ${sorted_file} ]]
    then
        mapped_reads=`samtools view -F 4 -c ${sorted_file}`
    else
        mapped_reads=0
    fi
fi

# count the unique mapped reads

f=${dir}/unique_reads_count_${genome}.tsv
    
if [[ ! -f ${f} ]]
then    
    if [[ ${mode} == "paired" ]]
    then
        merged_file=`find . | grep -P 'merged\.bam$' | grep ${genome} | grep ${mapper} | grep -v rose | head -1`
        
        echo "Counting uniquely mapped paired end reads in ${merged_file}, writing to ${f}..."

        # Merged reads are no longer paired so we only need to count
        # those that mapped and not worry abour -f 3 and whether
        # samples are correctly paired
        unique_reads=`samtools view -F 4 -c ${merged_file}`
        echo ${unique_reads} > ${f}
    else
        echo "Counting uniquely mapped reads..."
        unique_file=`find . | grep -P 'markdup\.bam$' | grep ${genome} | grep ${mapper} | grep -v rose | head -1`
        unique_reads=`samtools view -F 4 -c ${unique_file}`
    fi
    
    echo ${unique_reads} > ${f}.tmp
    mv ${f}.tmp ${f}
fi

unique_reads=`cat ${f}`


duplicate_reads=$((mapped_reads-unique_reads))
echo duplicate_reads ${duplicate_reads} ${mapped_reads} ${unique_reads}

f=${dir}/duplicate_reads_count_${genome}.tsv
if [[ ! -f ${f} ]]
then
    echo unique_reads ${duplicate_reads}
    echo ${duplicate_reads} > ${f}.tmp
    mv ${f}.tmp ${f}
fi

p_mapped_reads=`awk -v r=${reads} -v d=${mapped_reads} 'BEGIN {printf("%0.2f", d / r * 100)}'`
p_dup_reads=`awk -v r=${reads} -v d=${duplicate_reads} 'BEGIN {printf("%0.2f", d / r * 100)}'`
p_unique_reads=`awk -v r=${reads} -v u=${unique_reads} 'BEGIN {printf("%0.2f", u / r * 100)}'`

out=${dir}/alignment_quality_${genome}.tsv

echo -e "Name\tMode\tReads\tMapped Reads\t% Mapped Reads\tDuplicate Reads\t% Duplicate Reads\tUnique Reads\t% Unique Reads" > ${out}
echo -e ${sample}"\t"${mode}"\t"${reads}"\t"${mapped_reads}"\t"${p_mapped_reads}"\t"${duplicate_reads}"\t"${p_dup_reads}"\t"${unique_reads}"\t"${p_unique_reads} >> ${out}

#sys.stdout.write("\t".join([sample, "{:,}".format(total_reads), "{:,}".format(total_mapped_reads), "{:,}".format(non_duplicate_count), "{:.2f}".format(p)]) + "\n")  
#sys.stdout.write("\t".join([sample, "{:,}".format(total_reads), "{:,}".format(total_mapped_reads), "{:,}".format(duplicate_reads), "{:.2f}".format(p_dup), "{:,}".format(unique_reads), "{:.2f}".format(unique_reads_p)]) + "\n")  



#f=`ls | grep -P 'submit_align_bowtie.+\.e\d+' | sort | head -1`
#python /ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/alignment_quality.py ${sample} ${f} ${unique_reads} > alignment_quality.tsv
