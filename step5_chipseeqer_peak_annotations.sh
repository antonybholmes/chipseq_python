#export PATH=/ifs/scratch/cancer/Lab_RDF/ngs/apps/python/miniconda3/bin:$PATH

# a pattern
matches=$1
genome=$2

pwd=`pwd`
#/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/data/samples

#matches=`echo ${matches} | sed -r 's/\,/ /g'`

/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/chipseeqer_peak_annotation.sh ${matches} ${genome}


exit(0)

for match in `echo ${matches}`
do
    echo Searching for ${match}
    
    for d in `find ${pwd} -maxdepth 4 -type d | grep -P ${match} | grep -v chipseeqer | grep -v -i input | sort`
    do
      # skip if the folder does not contain a TF_targets file
      files=`find ${d} -type f | grep TF_targets | wc -l`
      
      if [[ ${files} -eq "0" ]]
      then
        continue
      fi
      
        
      sample=`echo ${d} | sed -r 's/^.+\///'`
      
      echo "########"
      echo "sample: ${sample}"
      echo "dir: ${d}"
      echo "########"
      
      cd ${d}/chipseeqer
      
      pwd2=`pwd`

      pwd      

      #for d2 in `ls`
      #do
      #  if [[ ! -d ${d2} ]]
      #  then
      #    continue
      #  fi
        
      #  cd ${d2}
        
      /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/chipseeqer_peak_annotation.sh ${genome}
        
      #  cd ${pwd2}
      #done
      
      cd ${pwd}
    done
done
