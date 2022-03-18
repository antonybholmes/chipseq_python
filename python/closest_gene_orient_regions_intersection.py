# -*- coding: utf-8 -*-
"""
Display regions in a gene oriented fashion, but only those peaks
that intersect

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re

import pychipseq.text
import pychipseq.human.regions
import pychipseq.headings



def getSampleId(text):
  return re.match(r'^.+ (.+)', text).group(1)
  

def two_way(file):
  gene_orientation = pychipseq.human.regions.ClosestGeneOrientatedRegionIntersections()
  
  # add some custom annotations for katia for testing etc
  
  gene_orientation.add_expression("RNA-seq hg19 GCvsN TPM", \
    pychipseq.expression.GeneExpression("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/ucsc_refseq_hg19_20171207/tpm_cb_vs_n/"))
  
  gene_orientation.add_expression("RNA-seq hg19 GCvsM TPM", \
    pychipseq.expression.GeneExpression("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/ucsc_refseq_hg19_20171207/tpm_cb_vs_m/"))
  
  gene_orientation.add_expression("RNA-seq hg19 GCvsN rmdup TPM", \
    pychipseq.expression.GeneExpression("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/ucsc_refseq_hg19_20171207/rmdup/tpm_cb_vs_n/"))
  
  gene_orientation.add_expression("RNA-seq hg19 GCvsM rmdup TPM", \
    pychipseq.expression.GeneExpression("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/ucsc_refseq_hg19_20171207/rmdup/tpm_cb_vs_m/"))
  
  
  gene_orientation.load_annotations(file)  
  
  f = open(file, 'r')
  
  # skip header
  header = f.readline().strip().split("\t")
  
  id_column = pychipseq.text.find_index(header, pychipseq.headings.REFSEQ_ID)
  closest_id_column = pychipseq.text.find_index(header, pychipseq.headings.CLOSEST_REFSEQ_ID)
  mir_column = pychipseq.text.find_index(header, pychipseq.headings.MIR_SYMBOL)
  sample_column = pychipseq.text.find_index(header, pychipseq.headings.SAMPLE)
  
  groupsi = set()
  
  for line in f:
    ls = line.strip()
    
    if len(ls) == 0:
      continue
    
    tokens = ls.split("\t")

 
    g1 = tokens[sample_column]
    g2 = tokens[sample_column + 1]
    
    if tokens[id_column] != pychipseq.text.NA:
      id = tokens[id_column]
    elif tokens[closest_id_column] != pychipseq.text.NA:
      id = tokens[closest_id_column]
    else:
      id = pychipseq.text.NA



    if id != pychipseq.text.NA:
      if g1 != pychipseq.text.NA and g2 != pychipseq.text.NA:
        groupsi.add(id)
        
    mir = tokens[mir_column]
    
    if mir != pychipseq.text.NA:
      if g1 != pychipseq.text.NA and g2 != pychipseq.text.NA:
        groupsi.add(mir)
      
    
  f.close()
  
  # Write custom header

  gene_orientation.print_header()  
  
  for id in gene_orientation.get_ids():
    if id in groupsi:
      gene_orientation.gene_orient_peak(id)
  
  for mir in gene_orientation.get_mirs():
    if mir in groupsi:
      gene_orientation.mir_orient_peak(mir)
      

def three_way(file):
  gene_orientation = pychipseq.human.regions.ClosestGeneOrientatedRegionIntersections()
  
  gene_orientation.load_annotations(file)
  
  f = open(file, 'r')

  # skip header
  header = f.readline().strip().split("\t")
  
  id_column = pychipseq.text.find_index(header, pychipseq.headings.REFSEQ_ID)
  closest_id_column = pychipseq.text.find_index(header, pychipseq.headings.CLOSEST_REFSEQ_ID)
  mir_column = pychipseq.text.find_index(header, pychipseq.headings.MIR_SYMBOL)
  sample_column = pychipseq.text.find_index(header, pychipseq.headings.SAMPLE)
  
  groupsi = set()
  
  for line in f:
    ls = line.strip()
    
    if len(ls) == 0:
      continue
    
    tokens = ls.split("\t")

 
    g1 = tokens[sample_column]
    g2 = tokens[sample_column + 1]
    g3 = tokens[sample_column + 2]
    
    if tokens[id_column] != pychipseq.text.NA:
      id = tokens[id_column]
    elif tokens[closest_id_column] != pychipseq.text.NA:
      id = tokens[closest_id_column]
    else:
      id = pychipseq.text.NA
    
    if id != pychipseq.text.NA:
      if g1 != pychipseq.text.NA and g2 != pychipseq.text.NA and g3 != pychipseq.text.NA:
        groupsi.add(id)
        
    mir = tokens[mir_column]
    
    if mir != pychipseq.text.NA:
      if g1 != pychipseq.text.NA and g2 != pychipseq.text.NA and g3 != pychipseq.text.NA:
        groupsi.add(mir)
      
    
  f.close()
  
  # Write custom header

  gene_orientation.print_header()  
  
  for id in gene_orientation.get_ids():
    if id in groupsi:
      gene_orientation.gene_orient_peak(id)
  
  for mir in gene_orientation.get_mirs():
    if mir in groupsi:
      gene_orientation.mir_orient_peak(mir)
      
      
def gene_orient_peaks(file):
  f = open(file, 'r')

  # skip header
  header = f.readline().strip().split("\t")
  
  sample_column = pychipseq.text.find_index(header, pychipseq.headings.SAMPLE)
  
  group3 = None
  
  if pychipseq.headings.SAMPLE in header[sample_column + 2]:
    group3 = getSampleId(header[sample_column + 2])
  
  f.close()
  
  if group3 is not None:
    three_way(file)
  else:
    two_way(file)


    
file = sys.argv[1]

if len(sys.argv) > 2:
  additional_annotations = sys.argv[2]
else:
  additional_annotations = ""

gene_orient_peaks(file)
