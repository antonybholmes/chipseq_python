#!/bin/bash
#$ -l mem=10G,time=8::
#$ -cwd
#$ -S /bin/bash


CHIPSEEQERDIR=/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/tools/ChIPseeqer/dist
export PATH=$PATH:$CHIPSEEQERDIR:$CHIPSEEQERDIR/SCRIPTS


echo
echo "==== Start "`date`" ===="
echo
echo "dir: `pwd`"
echo "os: `cat /etc/redhat-release`"
echo
echo "dir: ${input_dir}"

#lane=`echo ${chip_dir} | sed -r 's/^.+\///' | sed -r 's/\/$//'`
lane=`echo ${chip_dir} | sed -r 's/^.+\///'` # | sed -r 's/^\..*\///' | sed -r 's/\/.+//'`

# Remove trailing and leading forward slashes to leave just the
# file name. Remove the beginning of the file name to leave the name
# as Input_SeqId to keeps names short
#input_lane=`echo ${input_dir} | sed -r 's/^.+\///' | sed -r 's/\/$//' | sed -r 's/^.+Input/Input/'`


input_lane=`echo ${input_dir} | sed -r 's/^.+\///'`

#length=101

#echo ChIPseeqer.bin -chipdir ${chip_dir}/reads -inputdir ${input_dir}/reads -t 5 -fold_t 2 -readlen ${length} -minlen ${length} -format bowtiesam -chrdata /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/ChIPseeqer/data/hg19.chrdata -outfile TF_targets_${lane}_vs_${input_lane}_p5.txt
#ChIPseeqer.bin -chipdir ${chip_dir}/reads -inputdir ${input_dir}/reads -t 5 -fold_t 2 -readlen ${length} -minlen ${length} -format bowtiesam -chrdata /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/ChIPseeqer/data/hg19.chrdata -outfile TF_targets_${lane}_vs_${input_lane}_p5.txt

#genome=hg19
fold=2
#read_length=101
min_length=100
mindist=100 #6000 #100
mapper=hisat2


reads_dir=genome/${genome}/${mapper} #reads_${genome}

out=TF_targets_${lane}_vs_${input_lane}_p${p}_md${mindist}.tsv

cmd=(ChIPseeqer.bin \
-chipdir ${chip_dir}/${reads_dir} \
-inputdir ${input_dir}/${reads_dir} \
-t ${p} \
-fold_t ${fold} \
-readlen ${read_length} \
-minlen ${min_length} \
-mindist ${mindist} \
-format bowtiesam \
-chrdata /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/tools/ChIPseeqer/data/${genome}.chrdata \
-outfile ${out})

echo "${cmd[@]}"
eval ${cmd[@]}

#ChIPseeqer.bin -chipdir ${chip_dir}/${reads_dir} -inputdir ${input_dir}/${reads_dir} -t ${p} -fold_t ${fold} -readlen ${read_length} -minlen ${min_length} -mindist ${mindist} -format bowtiesam -chrdata /ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/ChIPseeqer/data/${genome}.chrdata -outfile TF_targets_${lane}_vs_${input_lane}_p${p}_md${mindist}.txt

/ifs/scratch/cancer/Lab_RDF/ngs/chip_seq/scripts/create_peak_bed_file.sh ${genome}


echo "==== End "`date`" ===="
echo
