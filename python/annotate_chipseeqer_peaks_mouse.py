# -*- coding: utf-8 -*-
"""
Annotate chipseeqer peaks

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re

import lib.genomic
import lib.mouse.genes


def annotate(file, type, prom_ext_5p, prom_ext_3p):
  sys.stderr.write("Annotating mouse " + file + "\n")  
  
  peak_annotation = lib.mouse.genes.AnnotatePeak(type, prom_ext_5p, prom_ext_3p)
   
  
  f = open(file, "r")
  
  #header
  sys.stdout.write("Genomic Location (mm10)")
  sys.stdout.write("\tP-value (ChIPseeqer)");
  #sys.stdout.write("\tMax Height (reads)");
  sys.stdout.write("\tScore (ChIPseeqer)");
  sys.stdout.write("\tPeak Width");
  sys.stdout.write("\t");
  
  peak_annotation.print_header()

  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t");

    chr = tokens[0]  
    start = int(tokens[1])
    end = int(tokens[2])
    width = end - start + 1
    
    # deal with infs
    if re.match(r'.*[Ii]nf.*', tokens[3]):
      p = -1000
    else:
      p = float(tokens[3])
    
    score = float(tokens[4])
    #max_height = int(tokens[6])
    
    location = lib.genomic.Location(chr, start, end) #chr + ":" + str(start) + "-" + str(end)

    #    
    # Write the initial portion
    #
    
    sys.stdout.write("\t".join([location.to_string(), str(p), str(score), str(width)]) + "\t")    
    
    #
    #  Write the gene annotation
    #
    
    peak_annotation.annotate(location)
    
  f.close()

  

file = sys.argv[1]
type = sys.argv[2]
prom_ext_5p = int(sys.argv[3])
prom_ext_3p = int(sys.argv[4])

annotate(file, type, prom_ext_5p, prom_ext_3p)
