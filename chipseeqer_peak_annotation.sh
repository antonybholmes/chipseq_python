#annotation_path=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/annotation/entrez
#annotation_path=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/annotation/rdf

PYTHONPATH=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/python
#PYTHONPATH=${PYTHONPATH}:${annotation_path}
#PYTHONPATH=${PYTHONPATH}:/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/annotation/entrez/gene_exp/cb_vs_n_m
#PYTHONPATH=${PYTHONPATH}:/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/annotation/entrez/gene_exp/icn_vs_gfp

export PYTHONPATH

matches=$1
genome=$2

#rm annotated_*
#rm *_dup_*
#rm Peaks*
#rm Genes*
#rm Closest_Genes*
#rm All_Genes*

prom_ext_5p=2000
prom_ext_3p=1000

echo "${genome}"

if [ "$genome" == "hg19" ]
then
	echo "Human"
	python ${PYTHONPATH}/annotate_chipseeqer_peaks.py --matches=${matches} --genome=${genome} --prom5=${prom_ext_5p} --prom3=${prom_ext_3p}
else
	python ${PYTHONPATH}/annotate_chipseeqer_peaks_mouse.py ${f} "Peak" ${prom_ext_5p} ${prom_ext_3p} > ${new_file}
fi

# bed files
/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/create_peak_bed_file.sh ${genome}
	
exit(0)

for f in `find . | grep -P 'TF_targets.+(tsv|txt)'`
do
	new_file=`echo ${f} | sed -r 's/txt/tsv/' | sed -r "s/\.tsv/_5p${prom_ext_5p}_3p${prom_ext_3p}_${genome}.tsv/" | sed -r 's/TF_targets./Peaks_/'`

	echo "Found ${f} -> ${new_file}"

	if [ "$genome" == "hg19" ]
	then
		echo "Human"
		python ${PYTHONPATH}/annotate_chipseeqer_peaks.py ${f} ${new_file} ${prom_ext_5p} ${prom_ext_3p}
	else
		python ${PYTHONPATH}/annotate_chipseeqer_peaks_mouse.py ${f} "Peak" ${prom_ext_5p} ${prom_ext_3p} > ${new_file}
	fi
  
	#debug
	#break
done

exit

for f in `find . | grep 'Peaks.+${genome}.(tsv|txt)'`
do
    id=`echo ${f} | sed -r 's/TF_targets_//' | sed -r 's/_vs.*//'`
  
    /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/chipseeqer_gene_annotation.sh ${f} ${genome}
  
	# Plot TSS
	/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/plot_peak_tss.sh ${f}
  
    # Plot widths
    /ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/plot_peak_widths.sh ${f}

	#debug
	#break
done

# excel files
#/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/chipseeqer_peak_annotation_excel_peak_files.sh ${genome}
#/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/chipseeqer_peak_annotation_excel_gene_files.sh ${genome}

#debug
#exit

# pie charts
#/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/peak_type_dist.sh


# PNG files
#/ifs/scratch/cancer/Lab_RDF/ngs/tools/scripts/batch_convert_pdf_to_png.sh .
