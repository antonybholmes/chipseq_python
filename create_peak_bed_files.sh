export SCRIPTDIR=/ifs/scratch/cancer/Lab_RDF/ngs/scripts/python
export PYTHONPATH=/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/python

genome=$1

echo Making bed files...

pwd=`pwd`

for f in `find . | grep -P 'TF_targets.*txt' | grep -v submission`
do
    dir=`dirname ${f}`
    file=`basename ${f}`
    id=`echo ${file} | sed -r 's/TF_targets.//' | sed -r 's/\.txt//'`
    
    echo "Going into ${dir}..."
    
    cd ${dir}
    
    bed=Peaks_${id}_${genome}.bed #graph

    echo Make bed $bed

    python3 ${PYTHONPATH}/create_peak_bed_file.py ${file} ${id} ${genome} > ${bed}
    
    md5=${bed}.md5sum
    md5sum ${bed} > ${md5}
    
    /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/bed2gff.sh ${bed}
    
    cd ${pwd}
    
    #break
done
