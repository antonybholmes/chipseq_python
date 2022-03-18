# -*- coding: utf-8 -*-
"""
Display peaks in a gene oriented fashion using the closest gene

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import pychipseq.human.peaks

def gene_orient_peaks(file):
  gene_orientation = pychipseq.human.peaks.ClosestGeneOrientatedPeaks("Peak")
  gene_orientation.load_annotations(file)
  
  gene_orientation.print_header()

  for id in gene_orientation.get_ids():
    gene_orientation.gene_orient_peak(id)
  
  for mir in gene_orientation.get_mirs():
    gene_orientation.mir_orient_peak(mir)
  
    
file = sys.argv[1]

if len(sys.argv) > 2:
  additional_annotations = sys.argv[2]
else:
  additional_annotations = ""

gene_orient_peaks(file)
