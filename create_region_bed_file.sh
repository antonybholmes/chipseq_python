function find_bed_file {
	file=$1
	
	dir=`echo $file | sed -r 's/\/[^\/]+$//'`
  echo $dir >&2
  p=`echo $file | grep -ohP 'p\d+'`
  echo $p >&2
	ret=`find ${dir} | grep -P '\.bed' | grep ${p} | sed -r 's/TF_targets/Peaks/' | head -1`
	
	echo $ret #$p $dir
}
	

file=$1
id=$2
f1=$3
f2=$4

echo Create Region BED File...

echo "f1 " $f1

# find the bed version of each file
bed_file_1=`find_bed_file $f1` #`echo ${f1} | sed -r 's/TF_targets/Peaks/' | sed -r 's/\.txt/.bed/'`
bed_file_2=`find_bed_file $f2` #`echo ${f2} | sed -r 's/TF_targets/Peaks/' | sed -r 's/\.txt/.bed/'`

echo "bed1 " $f1 $bed_file_1
echo "bed2 " $f2 $bed_file_2

#id=`echo ${file} | sed -r 's/^.*Peaks_//' | sed -r 's/_Overlaps\.txt//'`
out=`echo ${file} | sed -r 's/\.txt/_Overlaps.bed/'`

echo $file $id ${bed_file_1} ${bed_file_2} ${out}

python /ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/python/create_region_bed_file.py \
${file} \
${id} \
${bed_file_1} \
${bed_file_2} > ${out}
