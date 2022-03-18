# Overlap two sets of peaks and look at what is common


SCRIPT_PATH=/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/scripts
PYTHON_SCRIPT_PATH=${SCRIPT_PATH}/python
annotation_path=${PYTHON_SCRIPT_PATH}/annotation/entrez

export PYTHONPATH=${PYTHON_SCRIPT_PATH}
export PYTHONPATH=${PYTHONPATH}:${annotation_path}
export PYTHONPATH=${PYTHONPATH}:${PYTHON_SCRIPT_PATH}/annotation/entrez/gene_exp/cb_vs_n_m
#export PYTHONPATH=${PYTHONPATH}:${PYTHON_SCRIPT_PATH}/annotation/entrez/gene_exp/icn_vs_gfp

function remove_files {
  rm *.tsv
  rm *.xlsx
  rm *.pdf
  rm *.png
}

function plot {
  	id=$1
  	color=$2
  	peak_file=$3
  	peak_intersection_file=$4
  
 	 # plot the tss for the intersection
	${SCRIPT_PATH}/plot_region_tss.sh ${peak_file} ${id} ${color}
  	${SCRIPT_PATH}/plot_region_tss.sh ${peak_intersection_file} ${id} ${color}
	/ifs/scratch/cancer/Lab_RDF/abh2138/tools/scripts/batch_convert_pdf_to_png.sh .
}

function split_overlaps {
	tf1=$1
	tf2=$2
	id=$3
	color=$4
  
  	remove_files

  	prom_ext_5p=2000
  	prom_ext_3p=1000

  	peak_file=Peaks_${id}_Union.tsv
	peak_intersection_file=Peaks_${id}_Overlaps.tsv
	gene_union_file=Genes_${id}_Union.tsv
	gene_intersection_file=Genes_${id}_Overlaps.tsv
  	closest_gene_union_file=All_Genes_${id}_Union.tsv
	closest_gene_intersection_file=All_Genes_${id}_Overlaps.tsv
  
  
  # For making plots only, leave commented out otherwise
  # plot the tss for the intersection
	#${PYTHON_SCRIPT_PATH}/plot_region_tss.sh ${peak_file} ${id} ${color}
	#/ifs/scratch/cancer/Lab_RDF/abh2138/tools/scripts/batch_convert_pdf_to_png.sh .
  #plot ${id} ${color} ${peak_file} ${peak_intersection_file}
  #exit

	pwd=`pwd`
	name=`basename ${pwd}`
	annotated_peak_file=Peaks_${name}_annotated.tsv
  
  
	
	python ${PYTHON_SCRIPT_PATH}/overlapping_peaks.py \
	${tf1} \
	${tf2} > peak_overlaps.tsv
	
	echo peak_overlaps.tsv >&2
	
	types[0]=`head -1 peak_overlaps.tsv | cut -f6 | sed -r 's/\s+/_/g'`
	types[1]=`head -1 peak_overlaps.tsv | cut -f7 | sed -r 's/\s+/_/g'`
	
	
	echo ${types[0]} ${types[1]}
	cat peak_overlaps.tsv | sed 1d | cut -f1,6 | grep -P -v 'n/a' | cut -f1 | sort | uniq > ${types[0]}_peaks.tsv
	cat peak_overlaps.tsv | sed 1d | cut -f1,7 | grep -P -v 'n/a' | cut -f1 | sort | uniq > ${types[1]}_peaks.tsv
	
	python ${PYTHON_SCRIPT_PATH}/annotate_overlapping_peaks.py peak_overlaps.tsv ${prom_ext_5p} ${prom_ext_3p} ${annotated_peak_file}
  	python ${PYTHON_SCRIPT_PATH}/duplicate_regions.py ${annotated_peak_file} > duplicated_peak_overlaps.tsv
  	python ${PYTHON_SCRIPT_PATH}/gene_orient_regions_union.py duplicated_peak_overlaps.tsv > genes_peak_union.tsv
  
  	# Promoter only genes
  	python ${PYTHON_SCRIPT_PATH}/filter_promoters.py duplicated_peak_overlaps.tsv > duplicated_peak_overlaps_promoter.tsv
  	python ${PYTHON_SCRIPT_PATH}/gene_orient_regions_union.py duplicated_peak_overlaps_promoter.tsv > genes_peak_union_promoter.tsv
  
  	# Proximal only genes
  	python ${PYTHON_SCRIPT_PATH}/filter_proximal.py duplicated_peak_overlaps.tsv > duplicated_peak_overlaps_proximal.tsv
  	python ${PYTHON_SCRIPT_PATH}/gene_orient_regions_union.py duplicated_peak_overlaps_proximal.tsv > genes_peak_union_proximal.tsv
  

  	# Closest genes union
  	#python ${PYTHON_SCRIPT_PATH}/duplicate_closest_regions.py annotated_peak_overlaps.tsv > duplicated_closest_overlaps.tsv
	python ${PYTHON_SCRIPT_PATH}/closest_gene_orient_regions_union.py duplicated_peak_overlaps.tsv > closest_genes_peak_union.tsv
  
  
  
	
	# For the true intersection, remove all the non-intersecting peaks so that the closest TSS is truly the closest in the intersection etc
	#cat duplicated_peak_overlaps.tsv | awk '$5 !~ /n\/a/ && $6 !~ /n\/a/' > duplicated_peak_overlaps_intersections.tsv
	python ${PYTHON_SCRIPT_PATH}/gene_orient_regions_intersection.py duplicated_peak_overlaps.tsv > genes_peak_intersection.tsv
	
  	# Closest genes intersection
  	python ${PYTHON_SCRIPT_PATH}/closest_gene_orient_regions_intersection.py duplicated_peak_overlaps.tsv > closest_genes_peak_intersection.tsv
	
	
	
	
	#python ${PYTHON_SCRIPT_PATH}/finalize_output.py annotated_peak_overlaps.tsv > ${peak_file}
	#python ${PYTHON_SCRIPT_PATH}/finalize_output.py genes_peak_union.tsv > ${gene_union_file}
	#python ${PYTHON_SCRIPT_PATH}/finalize_output.py genes_peak_intersection.tsv > ${gene_intersection_file}
	
 	#python ${PYTHON_SCRIPT_PATH}/finalize_output.py closest_genes_peak_union.tsv > ${closest_gene_union_file}
	#python ${PYTHON_SCRIPT_PATH}/finalize_output.py closest_genes_peak_intersection.tsv > ${closest_gene_intersection_file}
	
	head -1 ${peak_file} > h
	# Only want the overlapping peaks
	cat ${peak_file} | sed 1d | awk '$5 == 2' > a
	cat h a > ${peak_intersection_file}
	
	echo ${peak_file} ${id} ${tf1} ${tf2}

	${SCRIPT_PATH}/create_region_bed_file.sh ${peak_file} \
	${id} \
	${tf1} \
	${tf2}
	
	
	#${SCRIPT_PATH}/xlsx_regions.sh ${peak_file}
	#${SCRIPT_PATH}/xlsx_regions.sh ${peak_intersection_file}
	#${SCRIPT_PATH}/xlsx_union_genes.sh ${gene_union_file}
	#${SCRIPT_PATH}/xlsx_regions_genes.sh ${gene_intersection_file}
  	#${SCRIPT_PATH}/xlsx_union_genes.sh ${closest_gene_union_file}
	#${SCRIPT_PATH}/xlsx_regions_genes.sh ${closest_gene_intersection_file}

  	#plot ${id} ${color} ${peak_file} ${peak_intersection_file}
}
