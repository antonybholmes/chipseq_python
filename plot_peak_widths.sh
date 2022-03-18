export SCRIPTDIR=/ifs/scratch/cancer/Lab_RDF/ngs/scripts/python
export PYTHONPATH=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/python


# Create TSS for peaks

f1=$1

tss_bin_size=0.1

#plot color blue
c1="#2c5aa0"


g1=`echo ${f1} | sed -r 's/^Peaks_//' | sed -r 's/_vs.+//'`

# first find the tss

out=`echo ${f1} | sed -r 's/\.tsv/_W.tsv/'`

export PYTHONPATH=${PYTHONPATH}/

# cut the closest tss column out

python ${SCRIPTDIR}/cut_col_by_name.py \
${f1} \
"Peak Width" \
> ${out}

# plot the ss
#python ${PYTHONPATH}/plot_region_tss.py ${out} ${g1} ${c1}
python ${PYTHONPATH}/plot_peak_widths_log10.py ${out} ${g1} ${c1}
