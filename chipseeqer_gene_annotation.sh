PYTHONPATH=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/python

f=$1
genome=$2


new_file=`echo ${f} | sed -r 's/Peaks_/dup_/'`
echo ${f} ${new_file}
python ${PYTHONPATH}/duplicate_peaks.py ${f} > ${new_file}

	
gene_file=`echo ${f} | sed -r 's/Peaks_/Genes_/'`
	
echo ${f} ${gene_file}

if [ "$genome" == "hg19" ]
then	
  python ${PYTHONPATH}/gene_orient_peaks.py ${new_file} > ${gene_file}
else
  python ${PYTHONPATH}/gene_orient_peaks_mouse.py ${new_file} > ${gene_file}
fi


# Group peaks by the closest gene
gene_file=`echo ${f} | sed -r 's/Peaks_/All_Genes_/'`
echo ${f} ${gene_file}

if [ "$genome" == "hg19" ]
then
  python ${PYTHONPATH}/closest_gene_orient_peaks.py ${new_file} > ${gene_file}
else
  python ${PYTHONPATH}/closest_gene_orient_peaks_mouse.py ${new_file} > ${gene_file}
fi
