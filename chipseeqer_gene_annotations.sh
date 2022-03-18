annotation_path=/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/annotation/entrez
#annotation_path=/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/annotation/rdf

PYTHONPATH=/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/
PYTHONPATH=${PYTHONPATH}:${annotation_path}
PYTHONPATH=${PYTHONPATH}:/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/annotation/entrez/gene_exp/cb_vs_n_m
#PYTHONPATH=${PYTHONPATH}:/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/annotation/entrez/gene_exp/icn_vs_gfp

export PYTHONPATH


match=$1

pwd=`pwd`

echo $match

for f in `find . | grep -P ${match}`
do

	/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/chipseeqer_gene_annotation.sh ${f}

	# Create TSS image
	/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/peak_tss.sh ${f}

	#debug
	#break
done

# excel files
/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/chipseeqer_peak_annotation_excel_gene_files.sh
