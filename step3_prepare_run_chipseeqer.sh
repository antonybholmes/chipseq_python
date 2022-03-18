# run multiple samples against their input

BATCH_SIZE=20

input=$1
# a pattern
matches=$2
read_length=$3
pvalues=$4
genome=$5

SAMPLE_DIR=`pwd` #/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/data/samples/hg19/david

input_dir=`find . -maxdepth 1 -type d | grep ${input} | head -1`
input_id=`basename ${input_dir}` # | sed -r 's/\.\///' | sed -r 's/\/.+//'`
input_dir=`find ${SAMPLE_DIR} -maxdepth 1 -type d | grep ${input} | head -1`

pwd=`pwd`

# Allow multiple matches to be parsed out
matches=`echo ${matches} | sed -r 's/\,/ /g'`

# turn comma separated list into space separated list
pvalues=`echo ${pvalues} | sed -r 's/\,/ /g'`

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

        analysis_dir=${sample_id}/chipseeqer/${input_id}
        mkdir -p ${analysis_dir}

        cd ${analysis_dir}
      
        #rm submit_run_chipseeqer.sh.*

        for p in `echo ${pvalues}`
        do
            echo ${sample_id} ${p} ${input_dir} ${d} ${analysis_dir}
        
            qsub -v chip_dir=${d},input_dir=${input_dir},read_length=${read_length},p=${p},genome=${genome} /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/qsub_chipseeqer.sh
        done
        #fi
        
        cd ${pwd}
        
        c=$((c+1))
        
        #if [[ ${c} -eq ${BATCH_SIZE} ]]
        #then
        #    break
        #fi
    done
done
