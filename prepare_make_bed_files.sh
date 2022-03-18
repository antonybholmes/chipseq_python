# run multiple samples against their input

BATCH_SIZE=20

matches=$1
genome=$2

SAMPLE_DIR=`pwd` #/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/data/samples/hg19/david

pwd=`pwd`

# Allow multiple matches to be parsed out
matches=`echo ${matches} | sed -r 's/\,/ /g'`

c=0

for match in `echo ${matches}`
do
    echo ${match}
    
    for d in `find ${SAMPLE_DIR} -maxdepth 1 -type d | grep ${match} | grep -v -i input`
    do
        #if [[ `find ${d} | grep TF_targets | wc -l` -eq 0 ]]
        #then
        
        # No target files, so lets make some
        
        echo dir ${d}
        
        sample_id=`basename ${d}` # echo ${d} | sed -r 's/\.\///' | sed -r 's/\/.+//'`
		
        analysis_dir=${SAMPLE_DIR}/${sample_id}/chipseeqer

        cd ${analysis_dir}
		
		for d2 in `find ${analysis_dir} -maxdepth 1 -type d`
		do
			cd ${d2}
			echo "here ${d2}"
			qsub -v genome=${genome} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_bed.sh
			cd ${analysis_dir}
        done
        
        cd ${pwd}
        
        c=$((c+1))
        
        #if [[ ${c} -eq ${BATCH_SIZE} ]]
        #then
        #    break
        #fi
    done
done
