# create some outputs

echo Converting gene files to Excel files...

for f in Genes*.txt
do
	/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/xlsx_peaks_genes.sh ${f}
done

for f in All_Genes*.txt
do
	/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/xlsx_peaks_genes.sh ${f}
done
