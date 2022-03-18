export SCRIPTDIR=/ifs/scratch/cancer/Lab_RDF/ngs/scripts/python
export PYTHONPATH=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/python



# Create TSS for peaks

f1=$1

tss_bin_size=0.1

#plot color blue
c1="#2c5aa0"

bins=`echo ${f1} | sed -r "s/\.tsv/_TSS_Bins_${tss_bin_size}_kb.tsv/"`
log_bins=`echo ${f1} | sed -r "s/\.tsv/_TSS_Log10_Bins_${tss_bin_size}.tsv/"`

echo bins $f1 $log_bins

g1=`echo ${f1} | sed -r 's/^Peaks_//' | sed -r 's/_vs.+//'`

# first find the tss

out=`echo ${f1} | sed -r 's/\.tsv/_TSS.tsv/'`

# cut the closest tss column out

python ${SCRIPTDIR}/cut_col_by_name.py \
${f1} \
"TSS Closest Distance" \
> ${out}


python ${PYTHONPATH}/tss_bin.py \
${out} \
${tss_bin_size} \
> ${bins}

python ${PYTHONPATH}/xlsx_tss.py ${bins}
python ${PYTHONPATH}/xlsx_tss.py ${bins} locked

python ${PYTHONPATH}/tss_log10_bin.py \
${out} \
${tss_bin_size} \
> ${log_bins}

python ${PYTHONPATH}/xlsx_tss_log10.py ${log_bins}
python ${PYTHONPATH}/xlsx_tss_log10.py ${log_bins} locked


# plot the ss
python ${PYTHONPATH}/plot_region_tss.py ${out} ${g1} ${c1}
python ${PYTHONPATH}/plot_region_tss_log10.py ${out} ${g1} ${c1}
