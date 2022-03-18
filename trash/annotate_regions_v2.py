# -*- coding: utf-8 -*-
"""
Find the common overlaps between a list of files reported as a common region

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys

import genes

def annotate(file):
  peak_annotation = genes.AnnotatePeak("Region")
  
  f = open(file, "r")

  line = f.readline().strip()

  sys.stdout.write(line);
  sys.stdout.write("\t");
  
  peak_annotation.print_header()


  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue

    sys.stdout.write(line + "\t")

    tokens = line.split("\t")
  
    location = tokens[0]
    
    peak_annotation.annotate(location)
    
    #sys.stdout.write("\n")
    
  f.close()
  
  
file = sys.argv[1]

annotate(file)