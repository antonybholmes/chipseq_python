

#export PATH=/ifs/scratch/cancer/Lab_RDF/ngs/apps/python/miniconda2/bin:$PATH


export PATH=${PATH}:/ifs/scratch/cancer/Lab_RDF/ngs/tools/weblogo/weblogo
export PATH=${PATH}:/ifs/scratch/cancer/Lab_RDF/ngs/tools/blat
export PATH=${PATH}:/ifs/scratch/cancer/Lab_RDF/ngs/tools/homer/bin
export PATH=${PATH}:/ifs/scratch/cancer/Lab_RDF/ngs/tools/samtools-1.8/bin
export PATH=${PATH}:/nfs/apps/R/3.1.2/:/nfs/apps/R/3.1.2/bin

export SAMPLE_DIR=/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/data/samples
export ROSE_PATH=/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/tools/ROSE
export PATH=${PATH}:${ROSE_PATH}
export PYTHONPATH=${PYTHONPATH}:${ROSE_PATH}

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/nfs/apps/boost/1.5.4/lib/

input=$1
matches=$2
p=$3
genome=$4
mode=$5

if [[ -z ${mode} ]]
then
  mode=single
fi

stitch=12500
tss=2000
# genome=hg19



search_dir=`pwd` #/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/data/samples/
depth=7

echo "find ${search_dir} -maxdepth ${depth} | grep -P ${input} | grep -v rose | grep -v superenhancers | grep -P 'merged.bam$' | head -1"

if [[ ${mode} == "single" ]]
then
    input_bam=`find ${search_dir} -maxdepth ${depth} | grep -P ${input} | grep -v rose | grep -v superenhancers | grep -P '(rmdup|markdup).bam$' | head -1`
else
    input_bam=`find ${search_dir} -maxdepth ${depth} | grep -P ${input} | grep -v rose | grep -v superenhancers | grep -P 'merged.bam$' | head -1` 
fi

echo ${input_bam} ${search_dir} ${input}


matches=`echo ${matches} | sed -r 's/\,/ /g'`

pwd=`pwd`

for match in `echo ${matches}`
do
    if [[ ${mode} == "single" ]]
    then
        chip_bam=`find ${search_dir} -maxdepth ${depth} | grep -P ${match} | grep -v rose | grep -v superenhancers | grep -P '(rmdup|markdup).bam$' | head -1`
    else
        chip_bam=`find ${search_dir} -maxdepth ${depth} | grep -P ${match} | grep -v rose | grep -v superenhancers | grep -P 'merged.bam$' | head -1`  
    fi
	
	echo "find ${search_dir} -maxdepth ${depth} | grep -P ${match} | grep -v rose | grep -P "p${p}" | grep -P '.gff$' | head -1"
    peaks_gff=`find ${search_dir} -maxdepth ${depth} | grep -P ${match} | grep -v rose | grep -P "p${p}" | grep -P '.gff$' | head -1`

    echo ${chip_bam}
    echo ${input_bam}
    echo ${peaks_gff}

    bam_dir=`dirname ${chip_bam}`

    echo "Changing to bam dir ${bam_dir}"
    
    cd ${bam_dir}

    #rm submit_run_rose*

    # make a temp dir to run rose in
    rose_dir=rose #_s${stitch}_t${tss}_p${p}

    #rm -rf ${rose_dir}
 
    # FOR TESTING
    #cd ${pwd}
    #continue
 
    chip_name=`echo ${chip_bam} | sed -r 's/^.+\///'`
    gff_name=`echo ${peaks_gff} | sed -r 's/^.+\///'`
    input_name=`echo ${input_bam} | sed -r 's/^.+\///'`
    
    mkdir -p ${rose_dir}/annotation
        
    #
    # Use symlinks to run scripts
    #
    
    # copy the rose files to this dir
    #cp ${ROSE_PATH}/*.py ${rose_dir}
    #cp ${ROSE_PATH}/*.R ${rose_dir}
    
    for f in `find ${ROSE_PATH} | grep -P '(py|R)$'`
    do
        echo ${f}
        name=`basename ${f}`
        ln -s ${f} ${rose_dir}/${name}
    done
    
    
    
    # copy annotation files
    for f in `find ${ROSE_PATH}/annotation | grep ucsc`
    do
        echo ${f}
        name=`basename ${f}`
        ln -s ${f} ${rose_dir}/annotation/${name}
    done
    

	if [[ -f "${chip_bam}" ]]
	then
		echo "chip ${chip_name}"
		ln -s ${chip_bam} ${rose_dir}/chip.bam
		ln -s ${chip_bam}.bai ${rose_dir}/chip.bam.bai
	fi
	
	if [[ -f "${input_bam}" ]]
	then
		echo "input ${input_name}"
		ln -s ${input_bam} ${rose_dir}/input.bam
		ln -s ${input_bam}.bai ${rose_dir}/input.bam.bai
	fi
	
	if [[ -f "${peaks_gff}" ]]
	then
		echo "gff ${gff_name}"
		ln -s ${peaks_gff} ${rose_dir}/${gff_name}
    fi
    
    echo "Changing to rose dir ${rose_dir}"
    cd ${rose_dir}
    


    cmd=(qsub \
    -N rose \
    -v stitch=${stitch},tss=${tss},chip_name=chip.bam,input_name=input.bam,gff_name=${gff_name},genome=${genome},p=${p} \
    /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_rose.sh)
    
    echo "${cmd[@]}"
    eval ${cmd[@]}
    
    cd ${pwd}

done
