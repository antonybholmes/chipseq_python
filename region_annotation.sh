#annotation_path=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/annotation/entrez
#annotation_path=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/annotation/rdf

DIR=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/python
export PYTHONPATH=${DIR}
export PYTHONPATH=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/python/pychipseq/:${PYTHONPATH}
#PYTHONPATH=${PYTHONPATH}:${annotation_path}
#PYTHONPATH=${PYTHONPATH}:/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/annotation/entrez/gene_exp/cb_vs_n_m
#PYTHONPATH=${PYTHONPATH}:/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts/annotation/entrez/gene_exp/icn_vs_gfp

f=$1
genome=$2


prom_ext_5p=2000
prom_ext_3p=1000


new_file=`echo ${f} | sed -r 's/txt/tsv/' | sed -r "s/\.tsv/_annotated_${genome}_5p${prom_ext_5p}_3p${prom_ext_3p}.tsv/"`

echo ${f} ${new_file}

if [ "$genome" == "hg19" ]
then
	python ${DIR}/annotate_chipseeqer_peaks.py ${f} ${new_file} ${prom_ext_5p} ${prom_ext_3p}
else
	python ${DIR}/annotate_chipseeqer_peaks_mouse.py ${f} "Peak" ${prom_ext_5p} ${prom_ext_3p} > ${new_file}
fi
