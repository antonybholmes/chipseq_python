#file=$1

# a pattern
matches=$1

read_length=$2 #36 #101 #$2

genome=$3

pwd=`pwd`

matches=`echo ${matches} | sed -r 's/\,/ /g'`

for match in `echo ${matches}`
do
    for f in `find . | grep ${match} | grep -P 'reads.chr1$'`
    do
      d=`dirname ${f}`
        
        echo ${d} ${read_length} ${genome}
      
      cd ${d}
      
      /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/create_meta.sh ${genome}
      
      #exit

      rm qsub_encode.sh.*
      
      qsub -v read_length=${read_length},genome=${genome},window=10 /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/qsub_encode.sh
      qsub -v read_length=${read_length},genome=${genome},window=100 /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/qsub_encode.sh
      qsub -v read_length=${read_length},genome=${genome},window=1000 /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/qsub_encode.sh
      qsub -v read_length=${read_length},genome=${genome},window=10000 /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/qsub_encode.sh
      qsub -v read_length=${read_length},genome=${genome},window=100000 /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/qsub_encode.sh
      qsub -v read_length=${read_length},genome=${genome},window=1000000 /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/qsub_encode.sh
      qsub -v read_length=${read_length},genome=${genome},window=10000000 /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/qsub_encode.sh
      qsub -v read_length=${read_length},genome=${genome},window=100000000 /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/qsub_encode.sh
      qsub -v read_length=${read_length},genome=${genome},window=1000000000 /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/qsub_encode.sh
      
      cd ${pwd}
    done
done
