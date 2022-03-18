matches=$1
genome=$2
matches=`echo ${matches} | sed -r 's/\,/ /g'`

for match in `echo ${matches}`
do
    #echo find . -maxdepth 3  grep -P ${match}  head -1
    f=`find . | grep -P ${match} | grep alignment_quality_${genome}.txt | head -1`
    head -1 ${f}
    break
done



for match in `echo ${matches}`
do
    for f in `find . | grep -P ${match} | grep alignment_quality_${genome}.txt | sort`
    do
        sample=`echo ${d} | sed -r 's/^.+\///'`
      
        cat ${f} | tail -1
    done
done
