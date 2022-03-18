# -*- coding: utf-8 -*-
"""
Find the common overlaps between a list of files reported as a common region

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys

import lib.genes
import lib.genomic
import lib.human.genes
import lib.mouse.genes

def annotate(file, genome, prom_ext_5p, prom_ext_3p):
  
  # = lib.human.genes.AnnotatePeak("Region")
  
  if genome == "hg19":
    peak_annotation = lib.human.genes.AnnotatePeak("Region", prom_ext_5p, prom_ext_3p)
  else:
    peak_annotation = lib.mouse.genes.AnnotatePeak("Region", prom_ext_5p, prom_ext_3p)
  
  
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
  
    location = lib.genomic.parse_location(tokens[0])
    
    peak_annotation.annotate(location)
    
    #sys.stdout.write("\n")
    
  f.close()
  
  
file = sys.argv[1]
genome = sys.argv[2]
prom_ext_5p = int(sys.argv[3])
prom_ext_3p = int(sys.argv[4])

annotate(file, genome, prom_ext_5p, prom_ext_3p)
