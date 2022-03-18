# -*- coding: utf-8 -*-
"""
Display peaks in a gene oriented fashion

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import pychipseq.human.peaks
import pychipseq.expression


file = sys.argv[1]

if len(sys.argv) > 2:
  additional_annotations = sys.argv[2]
else:
  additional_annotations = ""


gene_orientation = pychipseq.human.peaks.GeneOrientatedPeaks("Peak")

# Add some custom annotations

# RNA-seq tpm  
#gene_orientation.add_expression("RNA-seq GCvsN TPM", pychipseq.expression.RnaSeqTPMCBvsNExpression())
#gene_orientation.add_expression("RNA-seq GCvsM TPM", pychipseq.expression.RnaSeqTPMCBvsMExpression())

gene_orientation.add_expression("RNA-seq hg19 GCvsN TPM", \
  pychipseq.expression.GeneExpression("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/ucsc_refseq_hg19_20171207/tpm_cb_vs_n/"))

gene_orientation.add_expression("RNA-seq hg19 GCvsM TPM", \
  pychipseq.expression.GeneExpression("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/ucsc_refseq_hg19_20171207/tpm_cb_vs_m/"))

gene_orientation.add_expression("RNA-seq hg19 GCvsN rmdup TPM", \
  pychipseq.expression.GeneExpression("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/ucsc_refseq_hg19_20171207/rmdup/tpm_cb_vs_n/"))

gene_orientation.add_expression("RNA-seq hg19 GCvsM rmdup TPM", \
  pychipseq.expression.GeneExpression("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/ucsc_refseq_hg19_20171207/rmdup/tpm_cb_vs_m/"))

# Deseq2 test stuff
#gene_orientation.add_expression("RNA-seq SUD10 WTvsD83V DeSeq2", pychipseq.expression.RnaSeqGeneSUD10WTvsD83vDeSeq2Expression())
#gene_orientation.add_expression("RNA-seq SUD10 WTvsSTOP DeSeq2", pychipseq.expression.RnaSeqGeneSUD10WTvsStopDeSeq2Expression())
#gene_orientation.add_expression("RNA-seq SUD10 D83VvsSTOP DeSeq2", pychipseq.expression.RnaSeqGeneSUD10D83vsStopDeSeq2Expression())


# load file
gene_orientation.load_annotations(file)

gene_orientation.print_header()

for id in gene_orientation.get_ids():
  gene_orientation.gene_orient_peak(id)

for mir in gene_orientation.get_mirs():
  gene_orientation.mir_orient_peak(mir)
  
# Print a log of the annotations we used
gene_orientation.print_log()
