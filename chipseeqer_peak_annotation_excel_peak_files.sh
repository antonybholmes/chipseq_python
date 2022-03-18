# create some outputs

echo Converting genes...

for f in Peaks*.txt
do
	/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/scripts/xlsx_peaks.sh ${f}
done
