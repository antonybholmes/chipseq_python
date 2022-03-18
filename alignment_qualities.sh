# a pattern
match=$1
genome=$2

pwd=`pwd`

for d in `find samples -maxdepth 3 | grep -P ${match}`
do
	if [[ ! -d ${d} ]]
	then
		continue
	fi
	
  sample=`echo ${d} | sed -r 's/^.+\///'`
  
  echo ${d} ${sample}
  
  cd ${d}
  
  /ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/alignment_quality.sh ${sample} ${genome}
  
  cd ${pwd}
done
