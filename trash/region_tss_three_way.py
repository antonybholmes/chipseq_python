# -*- coding: utf-8 -*-
"""
Display regions in a gene oriented fashion, but only those peaks
that intersect from two specific samples in a three way comparison

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys

import annotation

def tss(file, sample_1, sample_2):
  f = open(file, 'r')

  # skip header
  header = f.readline().strip().split("\t")
  
  tss_column = annotation.find_heading_index(header, "Region TSS Closest Distance")
  p_column = annotation.find_heading_index(header, "Min P-value (ChIPseeqer)")
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


  sys.stdout.write("tss_dist\tp_value\n")
  
  for line in f:
    ls = line.strip()
    
    if len(ls) == 0:
      continue
    
    tokens = ls.split("\t")

 
    g1 = tokens[sample_column_1]
    g2 = tokens[sample_column_2]
    g3 = tokens[sample_column_3]
    
    #sys.stderr.write("asdjuhasdjhasdjhasd " + g1 + " " + g2 + " " + g3 + "\n")
    
    if g1 != "n/a" and g2 != "n/a" and g3 == "n/a":
      sys.stdout.write("\t".join([tokens[tss_column], tokens[p_column]]) + "\n")

  f.close()

    
file = sys.argv[1]
sample_1 = sys.argv[2]
sample_2 = sys.argv[3]

tss(file, sample_1, sample_2)