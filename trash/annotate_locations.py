# -*- coding: utf-8 -*-
"""
Annotate ad-hoc locations

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re

import genes
import lib_genomic


def annotate(file):
  peak_annotation = genes.AnnotatePeak("Location")
   
  sys.stderr.write("Annotating " + file + "\n")  
  
  f = open(file, "r")
  
  # skip header
  f.readline()
  
  #header
  sys.stdout.write("Genomic Location (hg19)")
  sys.stdout.write("\t");
  
  peak_annotation.print_header()

  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t");

    location = lib_genomic.parse_location(tokens[0])

    #    
    # Write the initial portion
    #
    
    sys.stdout.write("\t".join([location.to_string()]) + "\t")    
    
    #
    #  Write the gene annotation
    #
    
    peak_annotation.annotate(location)
    
  f.close()

  

file = sys.argv[1]

annotate(file)
