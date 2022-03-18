f1=$1
g1=$2
c1=$3

python_script_dir=/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/python

echo Region TSS...

# first find the tss

tss_bin_size=0.1

count=`cat ${f1} | sed 1d | wc -l`

bins=`echo ${f1} | sed -r "s/\.txt/_TSS_Bins_${tss_bin_size}_kb_n${count}.txt/"`
log_bins=`echo ${f1} | sed -r "s/\.txt/_TSS_Log10_Bins_${tss_bin_size}_n${count}.txt/"`

out=`echo ${f1} | sed -r "s/\.txt/_TSS_n${count}.txt/"`

# try to shrink the group by removing misc stuff so it fits better on the label
g1=`echo ${g1} | sed -r 's/CB[0-9]_//g'`


# cut the closest tss column out

echo $f1 ${bins} ${log_bins}

python /ifs/scratch/cancer/Lab_RDF/abh2138/scripts/python/cut_col_by_name.py \
${f1} \
"Region TSS Closest Distance" \
"Best P-Value (ChIPseeqer)" \
"Number Of Overlapping" \
> a

#cat a | sed 1d | awk '{p=($2+$3)/2;print $1"\t"p"\t"$4}' > b
cp a b

head -1 a | cut -f1,2,3 > h
cat b | awk '$3 == 2' > c
cat h c > ${out}
#rm a b h

# bin the file

python ${python_script_dir}/tss_bin.py \
${out} \
${tss_bin_size} \
> ${bins}
python ${python_script_dir}/xlsx_tss.py ${bins}
python ${python_script_dir}/xlsx_tss.py ${bins} locked

python ${python_script_dir}/tss_log10_bin.py \
${out} \
${tss_bin_size} \
> ${log_bins}
python ${python_script_dir}/xlsx_tss_log10.py ${log_bins}
python ${python_script_dir}/xlsx_tss_log10.py ${log_bins} locked


# plot the ss
python ${python_script_dir}/plot_region_tss.py ${out} ${g1} ${c1}
python ${python_script_dir}/plot_region_tss_log10.py ${out} ${g1} ${c1}

#Rscript /ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/plot_region_tss.R ${out} ${g1}
#Rscript /ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/plot_region_tss_log10_dist.R ${out} ${g1}
