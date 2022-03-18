# -*- coding: utf-8 -*-
"""
Find the common overlaps between a list of files reported as a common region

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys

import lib.human.genes
import lib.genomic

def annotate(file, prom_ext_5p=2000, prom_ext_3p=1000):
  peak_annotation = lib.human.genes.AnnotatePeak("Region", prom_ext_5p, prom_ext_3p)
  
  f = open(file, "r")

  line = f.readline().strip()

  sys.stdout.write(line);
  sys.stdout.write("\t");
  
  peak_annotation.print_header()
  sys.stdout.write("\n")

  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue

    sys.stdout.write(line + "\t")

    tokens = line.split("\t")
  
    location = lib.genomic.parse_location(tokens[0])
    
    sys.stdout.write("\t".join(peak_annotation.annotate(location)))
    
    sys.stdout.write("\n")
    
  f.close()
  
  
file = sys.argv[1]
prom_ext_5p = int(sys.argv[2])
prom_ext_3p = int(sys.argv[3])

annotate(file, prom_ext_5p, prom_ext_3p)