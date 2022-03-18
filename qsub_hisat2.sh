#!/bin/bash -l
#$ -l mem=16G,time=16::
#$ -cwd
#$ -S /bin/bash

export DIR=/ifs/scratch/cancer/Lab_RDF/ngs
export SCRIPT_DIR=${DIR}/chip_seq/scripts
export HISAT2_INDEXES=${DIR}/genomes/indexes/hisat2/${genome}/
#export HISAT2=${DIR}/tools/hisat2-2.0.5/hisat2
#export SAMTOOLS=${DIR}/tools/samtools-1.8/bin/samtools


#CHROMS=(chr5 chr6) 
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)


echo
echo "==== Start "`date`" ===="
echo
echo "dir: `pwd`"
echo "os: `cat /etc/redhat-release`"
echo

echo "Mode: ${mode}"

if [[ ${mode} == "paired" ]]
then
	# tail means that if the trimmed versions exist, we use those in
	# preference to the regular fastq
	file1=`ls fastq/*.fastq.gz | grep R1 | tail -1`
	file2=`ls fastq/*.fastq.gz | grep R2 | tail -1`
	files=${file1},${file2}
	echo "Found pairs ${file1} and ${file2}..."
else
	file1=`ls fastq/*.gz | sort | tail -1`
	#file1=`ls fastq/*.fastq.gz | tail -1`
	echo "Found ${file1}..."
fi

# If concat of reads exists, use that in preference to other files
#concat=`ls fastq/concat*.fastq.gz | sort | head -1`

#if [[ -f "${concat}" ]]
#then
    #file1=${concat}
    #file2=`ls fastq/concat*.fastq.gz | sort | tail -1`
    #files=${file1},${file2}
#fi


#output=`echo ${fastq} | sed -r 's/ .+//' | sed -r 's/\.fastq\.gz//'`

pwd=`pwd`
lib=`basename ${pwd}`

echo "Library ${lib}"

# create sam file name

# make a directory to store our results
# output_dir=../${genome}/align_2_mismatches/reads/${lib}
# output_dir=../${genome}/align_2_mismatches/reads
fastq_dir=${pwd}/fastq
genome_dir=${pwd}/genome/${genome}



output_dir=${genome_dir}/hisat2
echo Making directory ${output_dir}...
mkdir --parents ${output_dir}

bam=${lib}_${genome}.bam
sorted_bam=`echo ${bam} | sed -r 's/\.bam/.sorted.bam/'`
name_sort_bam=`echo ${bam} | sed -r 's/\.bam/.name_sorted.bam/'`
fixmate_bam=`echo ${bam} | sed -r 's/\.bam/.fixmate.bam/'`
sorted_fixmate_bam=`echo ${bam} | sed -r 's/\.bam/.fixmate_sorted.bam/'`
markdup_bam=`echo ${sorted_fixmate_bam} | sed -r 's/\.bam/.markdup.bam/'`
unique_bam=`echo ${markdup_bam} | sed -r 's/\.bam/.unique.bam/'`
merge_bam=`echo ${unique_bam} | sed -r 's/\.bam/.merged.bam/'`

bam=1.${bam}
sorted_bam=2.${sorted_bam}
name_sort_bam=3.${name_sort_bam}
fixmate_bam=4.${fixmate_bam}
sorted_fixmate_bam=5.${sorted_fixmate_bam}
markdup_bam=6.${markdup_bam}
unique_bam=7.${unique_bam}
merge_bam=8.${merge_bam}

bam=${output_dir}/${bam}
sorted_bam=${output_dir}/${sorted_bam}
name_sort_bam=${output_dir}/${name_sort_bam}
fixmate_bam=${output_dir}/${fixmate_bam}
sorted_fixmate_bam=${output_dir}/${sorted_fixmate_bam}
markdup_bam=${output_dir}/${markdup_bam}
unique_bam=${output_dir}/${unique_bam}
merge_bam=${output_dir}/${merge_bam}


trim5=0 #10
trim3=0 #10

# Allow multiple, then find unique
k=5 #1

# mismatches is always a negative number
mismatches="-"`echo ${mismatches} | sed -r 's/^-//'`

# mismatches should be given as -1,-2,-3 etc

t1=`echo ${file1} | sed -r 's/.gz/.trimmed.gz/'`
#t1_tmp=`mktemp ${t1}.tmp.XXXXXXXX`


if [[ ${mode} == "paired" ]]
then
  # paired end mode
  
  t2=`echo ${file2} | sed -r 's/.gz/.trimmed.gz/'`
  #t2_tmp=`mktemp ${t2}.tmp.XXXXXXXX`
  
  echo "Paired end mode: ${file1} ${file2}"
  
  trim_cmd=(cutadapt \
    -q 10,10 \
    --cores=${cores} \
    --pair-filter=any \
    --minimum-length 50 \
    -o ${t1_tmp} \
    -p ${t2_tmp} \
    ${file1} \
    ${file2})
  
  # paired
  cmd=(hisat2 \
  --threads ${cores} \
  --add-chrname \
  --score-min C,${mismatches} \
  --rdg 100,100 \
  --rfg 100,100 \
  --mp 1,1 \
  --no-softclip \
  --no-mixed \
  --no-discordant \
  --trim5 ${trim5} \
  --trim3 ${trim3} \
  --ignore-quals \
  --no-spliced-alignment \
  -q \
  -x ${genome} \
  -k ${k} \
  -1 ${file1} \
  -2 ${file2}) # ${t1} ${t2}
  echo test "${cmd[@]}"
else
  #single
  
  echo "Single end mode: ${file1}"
  
  trim_cmd=(cutadapt \
    -q 10,10 \
    --cores=${cores} \
    --minimum-length 50 \
    -o {t1_tmp} \
    ${file1})
  
  cmd=(hisat2 \
  --threads ${cores} \
  --add-chrname \
  --score-min C,${mismatches} \
  --rdg 100,100 \
  --rfg 100,100 \
  --mp 1,1 \
  --no-softclip \
  --trim5 ${trim5} \
  --trim3 ${trim3} \
  --no-spliced-alignment \
  --ignore-quals \
  -q \
  -x ${genome} \
  -k ${k} \
  -U ${file1})
fi

# Trim reads if trimmed files do not exist
#if [[ ! -f ${t1} ]]
#then
    #echo "${trim_cmd[@]}"
    #eval ${trim_cmd[@]}
    
    #mv ${t1_tmp} ${t1}
    #mv ${t2_tmp} ${t2}
#fi


#exit

# Align trimmed
if [[ ! -f "${bam}" ]]
then
    tmp=${bam}.tmp
    echo "Aligning..."
    echo "${cmd[@]}" ${tmp}
    
    eval ${cmd[@]} | grep -v -P '\tchrUn' | samtools view -F 4 -Sb - > ${tmp}
    mv ${tmp} ${bam}
fi

if [[ -f "${bam}" && ! -f "${sorted_bam}" ]]
then
    echo "Sorting: ${bam} -> ${sorted_bam}..."
    tmp=${sorted_bam}.tmp
    samtools sort -@ ${cores} -o ${tmp} ${bam}
    mv ${tmp} ${sorted_bam}
fi

if [[ -f "${sorted_bam}" && ! -f "${sorted_bam}.bai" ]]
then
    echo "Indexing ${sorted_bam}..."
    cmd=(samtools index ${sorted_bam})
    echo "${cmd[@]}"
    eval ${cmd[@]}
fi


# Convert to sorted bam

echo
echo "Checking for ${name_sort_bam}..."

if [[ ! -f "${name_sort_bam}" && ! -f "${markdup_bam}" && ! -f "${merge_bam}" ]]
then
    echo "Sorting by name: ${bam} -> ${name_sort_bam}..."
    tmp=${name_sort_bam}.tmp
    samtools sort -n -@ ${cores} -o ${tmp} ${bam}
    mv ${tmp} ${name_sort_bam}
fi

echo
echo "Checking for ${fixmate_bam}..."

if [[ ! -f "${fixmate_bam}" && ! -f "${markdup_bam}" && ! -f "${merge_bam}" ]]
then
    echo "Running fixmate: ${name_sort_bam} -> ${fixmate_bam}..."
    tmp=${fixmate_bam}.tmp
    samtools fixmate -m ${name_sort_bam} ${tmp}
    mv ${tmp} ${fixmate_bam}
    
    # Remove to save space
    #rm ${name_sort_bam}
fi

# Convert to sorted bam
echo
echo "Checking for ${sorted_fixmate_bam}..."

if [[ ! -f "${sorted_fixmate_bam}" && ! -f "${markdup_bam}" && ! -f "${merge_bam}" ]]
then
    echo "Sorting: ${fixmate_bam} -> ${sorted_fixmate_bam}..."
    tmp=${sorted_fixmate_bam}.tmp
    samtools sort -@ ${cores} -o ${tmp} ${fixmate_bam}
    mv ${tmp} ${sorted_fixmate_bam}
    
    # Free space
    #rm ${fixmate_bam}
fi




# remove pcr duplicates
# rmdup_bam=`echo ${sorted_fixmate_bam} | sed -r 's/\.bam/.rmdup.bam/'`
# echo Removing duplicates ${sorted_fixmate_bam} ${rmdup_bam}
# samtools rmdup -s ${sorted_fixmate_bam} ${rmdup_bam}

echo
echo "Checking for ${markdup_bam}..."

if [[ ! -f "${markdup_bam}" && ! -f "${merge_bam}" ]]
then
    echo "Removing duplicates: ${sorted_fixmate_bam} -> ${markdup_bam}..."
    tmp=${markdup_bam}.tmp
    samtools markdup -r -s ${sorted_fixmate_bam} ${tmp}
    mv ${tmp} ${markdup_bam}

    # Free space
    #rm ${sorted_fixmate_bam}
fi

if [[ -f "${markdup_bam}" && ! -f "${markdup_bam}.bai" ]]
then
    echo "Indexing ${markdup_bam}..."
    samtools index ${markdup_bam}
fi

echo
echo "Checking for ${unique_bam}..."

if [[ ! -f "${unique_bam}" && ! -f "${merge_bam}" ]]
then
    echo "Finding unique reads in ${markdup_bam}..."
    
    tmp=${unique_bam}.tmp
    
    # match on the header (@) or tag
    if [[ ${mode} == "paired" ]]
    then
        # Keep proper pairs with unique alignments
        samtools view -h -f 3 ${markdup_bam} | grep -P '(^@|NH:i:1)' | samtools view -Sb - > ${tmp}
    else
        # single
        samtools view -h -F 4 ${markdup_bam} | grep -P '(^@|NH:i:1)' | samtools view -Sb - > ${tmp}
    fi
    
    mv ${tmp} ${unique_bam}
    
    # Free space
    #rm ${markdup_bam}
fi

if [[ -f "${unique_bam}" && ! -f "${unique_bam}.bai" ]]
then
    echo "Indexing ${unique_bam}..."
    cmd=(samtools index ${unique_bam})
    echo "${cmd[@]}"
    eval ${cmd[@]}
fi


#if [[ ! -f unique_reads.txt ]]
#then
    #echo Counting unique reads...
    ## for large files cannot sort or get unique
    #samtools view -F 4 -c ${unique_bam} > unique_reads.txt
#fi

if [[ ${mode} == "single" && -f "${unique_bam}" && -f "${unique_bam}.bai" ]]
then
	echo "Linking single unique file to merge..."
	# use symlink to make merge file from unique to save time
	ln -s ${unique_bam} ${merge_bam}
	ln -s ${unique_bam}.bai ${merge_bam}.bai
fi

echo
echo "Checking for ${merge_bam}..."

if [[ ${mode} == "paired" && ! -f "${merge_bam}" ]]
then
    echo "Merging pairs in ${merge_bam}..."
    
    #tmp=`mktemp ${merge_bam}.tmp.XXXXXXXX`
    #cmd="python3 /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/python/merge_pairs/merge_pairs.py \
    #--samtools=samtools \
    #--out=${tmp} \
    #${markdup_bam}"
    #echo ${cmd}
    #${cmd}
    #mv ${tmp} ${merge_bam}
    
    #
    # Merge reads in each chr separately (read pairs should not span
    # multiple chromosomes so no need to address this case) so that
    # each pair becomes one read (with missing gap inserted or overlap
    # accounted for so that it is not doubled).
    #
    
    id=1
    
    for chr in ${CHROMS[*]}
    do
        echo ${chr}
        
        # tmp=`mktemp ${merge_bam}.tmp.XXXXXXXX`
        # qsub -N merge_${chr} -v bam=${unique_bam},chr=${chr},id=${id},cores=${cores} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_merge_chr.sh
        cmd=(qsub -N ${chr}_merge -v bam=${unique_bam},chr=${chr},id=${id},cores=${cores} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_merge_chr.sh)
        echo "${cmd[@]}"
        eval ${cmd[@]}
        
        id=$((id+1))
    done
    
    # once all jobs complete, merge the chromosomes into one
    cmd=(qsub -N merge_all -hold_jid "chr*_merge" -v bam=${merge_bam},lib=${lib},genome=${genome},cores=${cores} -pe smp ${cores} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_merge_chrs.sh)
    echo "${cmd[@]}"
    
    eval ${cmd[@]}
    #qsub -N merge_all -hold_jid "chr*_merge" -v bam=${merge_bam},lib=${lib},genome=${genome},cores=${cores} -pe smp ${cores} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_merge_chrs.sh
else
    echo "Merge found, running alignment quality..."
    
    if [[ ! -f "${merge_bam}.bai" ]]
    then
        echo "Indexing ${merge_bam}..."
        cmd=(samtools index ${merge_bam})
        # echo "${cmd[@]}"
        eval ${cmd[@]}
    fi

    # in the single end read case, or when the merge file has been, 
    # previously generate test the alignment quality
    
    cmd=(${SCRIPT_DIR}/alignment_quality.sh ${lib} ${genome} hisat2)
    echo "${cmd[@]}"
    eval ${cmd[@]}
fi


echo "==== End "`date`" ===="
echo
