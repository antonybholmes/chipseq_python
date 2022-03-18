# -*- coding: utf-8 -*-
"""
Display regions in a gene oriented fashion, but only those peaks
that intersect from two specific samples in a three way comparison

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re

import annotation
import genes



def getSampleId(text):
  return re.match(r'^.+ (.+)', text).group(1)
     

def gene_orient_peaks(file, sample_1, sample_2, additional_annotations):
  gene_orientation = genes.GeneOrientatedPeaks(file, "Region", additional_annotations)
  
  f = open(file, 'r')

  # skip header
  header = f.readline().strip().split("\t")
  
  entrez_column = annotation.find_heading_index(header, "Gene Entrez ID")
  mir_column = annotation.find_heading_index(header, "miR Symbol")
  sample_column = annotation.find_heading_index(header, "Sample")
  
  sample_column_1 = annotation.find_heading_index(header, sample_1)
  sample_column_2 = annotation.find_heading_index(header, sample_2)

  # determine which sample column not in use
  is_s1 = sample_column == sample_column_1 or sample_column == sample_column_2
  is_s2 = (sample_column + 1) == sample_column_1 or (sample_column + 1) == sample_column_2

  if not is_s1:
    sample_column_3 = sample_column
  elif not is_s2:
    sample_column_3 = sample_column + 1
  else:
    sample_column_3 = sample_column + 2

  sys.stderr.write(sample_1 + " " + sample_2 + " " + str(sample_column_1) + " " + str(sample_column_2) + " " + str(sample_column_3) + "\n")
  
  groupsi = set()
  
  c = 0
  for line in f:
    ls = line.strip()
    
    if len(ls) == 0:
      continue
    
    tokens = ls.split("\t")

 
    g1 = tokens[sample_column_1]
    g2 = tokens[sample_column_2]

    
    
    entrez = tokens[entrez_column]
    
    
    
    #if entrez != "n/a":
    #  if g1 != "n/a" and g2 != "n/a" and g3 == "n/a":
    #    groupsi.add(entrez)
    
    # overlap regardless of whether we
    if entrez != "n/a" and g1 != "n/a" and g2 != "n/a":
      groupsi.add(entrez)
      #sys.stderr.write(str(len(groupsi)) + " " + g1 + " " + g2 + " " + entrez + "\n")
        
    mir = tokens[mir_column]
    
    #if mir != "n/a":
    #  if g1 != "n/a" and g2 != "n/a" and g3 == "n/a":
    #    groupsi.add(mir)
    
    if mir != "n/a" and  g1 != "n/a" and g2 != "n/a":
      groupsi.add(mir)
    
  f.close()
  
  # Write custom header

  gene_orientation.print_header()  
  
  for entrez in gene_orientation.get_entrezes():
    if entrez in groupsi:
      gene_orientation.gene_orient_peak(entrez)
  
  for mir in gene_orientation.get_mirs():
    if mir in groupsi:
      gene_orientation.mir_orient_peak(mir)

    
file = sys.argv[1]
sample_1 = sys.argv[2]
sample_2 = sys.argv[3]

if len(sys.argv) > 4:
  additional_annotations = sys.argv[4]
else:
  additional_annotations = ""

gene_orient_peaks(file, sample_1, sample_2, additional_annotations)