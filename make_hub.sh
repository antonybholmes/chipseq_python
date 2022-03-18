name=$1
genome=$2
matches=$3

pwd=`pwd`

dir=hub/${name}

# Make hub exist
mkdir -p ${dir}
mkdir -p ${dir}/${genome}

matches=`echo ${matches} | sed -r 's/\,/ /g'`

# copy files to make hub
for match in `echo ${matches}`
do
    for f in `find ${pwd} | grep ${match} | grep -P '(bb|bw)$' | grep -v hub`
    do
        echo Copying ${f} to ${dir}/${genome}
        #cp ${f} ${dir}/${genome}
		ln -s ${f} ${dir}/${genome}/
    done
done

#for match in `echo ${matches}`
#do
#    for f in `find . | grep ${match} | grep -P 'bb$' | grep -v hubs`
#    do
#        echo Copying ${f} to ${dir}/${genome}
#        cp ${f} ${dir}/${genome}
#    done
#done

# make hub file
hub=${dir}/hub.txt
echo hub ${name} > ${hub}
echo shortLabel ${name} >> ${hub}
echo longLabel ${name} >> ${hub}
echo genomesFile genomes.txt >> ${hub}
echo email abh2138@cumc.columbia.edu >> ${hub}
echo descriptionUrl ${name}.html >> ${hub}


# make genome file

echo -e "genome ${genome}" > ${dir}/genomes.txt
echo -e "trackDb hg19/trackDb.txt" >> ${dir}/genomes.txt

# make trackdb file
trackdb=${dir}/${genome}/trackDb.txt

if [[ -f ${trackdb} ]]
then
    rm ${trackdb}
fi

for f in ${dir}/${genome}/*.bw
do
  echo ${f}
  
  base=`basename ${f}`
  n=`echo ${base} | sed -r 's/\..+//' | sed -r 's/Peaks_//'`
  
  #cat template.track | sed -r "s/insert_name/${name}/g" | sed -r "s/insert_file/${f}/g" > ${n}.track
  
  
  echo track ${n} >> ${trackdb}
  echo bigDataUrl ${base} >> ${trackdb}
  echo shortLabel ${n} >> ${trackdb}
  echo longLabel ${n} >> ${trackdb}
  echo type bigWig >> ${trackdb}
  echo >> ${trackdb}
done

for f in ${dir}/${genome}/*.bb
do
  echo ${f}
  
  base=`basename ${f}`
  n=`echo ${base} | sed -r 's/\..+//'`
  
  echo track ${n} >> ${trackdb}
  echo bigDataUrl ${base} >> ${trackdb}
  echo shortLabel ${n} >> ${trackdb}
  echo longLabel ${n} >> ${trackdb}
  echo type bigBed >> ${trackdb}
  echo >> ${trackdb}
done

aws s3 sync ${dir}/ s3://files.rdf-lab.org/ucsc/hubs/chipseq/human/hg19/bpm/${name}/
