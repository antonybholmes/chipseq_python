# -*- coding: utf-8 -*-
"""
Finds the regions shared by the three intersections between three samples

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re

import lib.human.regions
import lib.headings

import lib.text


def getSampleId(text):
  return re.match(r'^.+ (.+)', text).group(1)
     

def gene_orient_peaks(file):
  gene_orientation = lib.human.regions.GeneOrientatedRegionCores()
  gene_orientation.load_annotations(file)
  
  f = open(file, 'r')

  # skip header
  header = f.readline().strip().split("\t")
  
  id_column = lib.text.find_index(header, lib.headings.REFSEQ_ID)
  mir_column = lib.text.find_index(header, lib.headings.MIR_SYMBOL)
  sample_column_1 = lib.text.find_index(header, lib.headings.SAMPLE)
  sample_column_2 = sample_column_1 + 1
  sample_column_3 = sample_column_1 + 2

  groupsi = set()
  
  for line in f:
    ls = line.strip()
    
    if len(ls) == 0:
      continue
    
    tokens = ls.split("\t")

 
    g1 = tokens[sample_column_1]
    g2 = tokens[sample_column_2]
    g3 = tokens[sample_column_3]
    
    #sys.stderr.write(g1 + " " + g2 + " " + g3 + "\n")
    
    id = tokens[id_column]
    
    if id != lib.text.NA:
      if (g1 != lib.text.NA and g2 != lib.text.NA) or (g1 != lib.text.NA and g3 != lib.text.NA) or (g2 != lib.text.NA and g3 != lib.text.NA):
        groupsi.add(id)
        
    mir = tokens[mir_column]
    
    if mir != lib.text.NA:
      if (g1 != lib.text.NA and g2 != lib.text.NA) or (g1 != lib.text.NA and g3 != lib.text.NA) or (g2 != lib.text.NA and g3 != lib.text.NA):
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

    
file = sys.argv[1]

gene_orient_peaks(file)
