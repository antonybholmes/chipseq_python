export SCRIPTDIR=/ifs/scratch/cancer/Lab_RDF/ngs/scripts/python
export PYTHONPATH=/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/python

genome=$1
pwd=`pwd`

echo "Making bed files in ${pwd} for ${genome}..."

for file in `find . | grep -P 'TF_targets.+(txt|tsv)' | grep -v submission`
do
    dir=`dirname ${file}`
    name=`basename ${file}`
    id=`echo $name | sed -r 's/TF_targets.//' | sed -r 's/\.txt//' | sed -r 's/\.tsv//'`
    bed=Peaks_${id}_${genome}.bed #graph

    echo "Make bed $bed"

    python ${PYTHONPATH}/create_peak_bed_file.py ${file} ${id} ${genome} > ${bed}
  
    /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/bed2gff.sh ${bed}
done
