# -*- coding: utf-8 -*-
"""
Generate a tss distribution for a region file

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import collections
import re
import os

import annotation

def closest(file):
  bins = collections.defaultdict(int)

  f = open(file, "r")

  sys.stdout.write("Genomic Location (hg19)")
  sys.stdout.write("\tP-value (ChIPseeqer)");
  sys.stdout.write("\tMax Height (ChIPseeqer)");
  sys.stdout.write("\tPeak Size");
  sys.stdout.write("\tTSS Distance");
  sys.stdout.write("\tTSS Entrez ID");
  sys.stdout.write("\n");
  
  total_peaks = 0
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    if tokens[3] != "2":
      continue
    
    location = tokens[0]
    
    matcher = re.match(r'(chr.+):(\d+)-(\d+)', location)
    
    chr = matcher.group(1)
    start = int(matcher.group(2))
    end = int(matcher.group(3))
    width = end - start + 1
    
    p = -1000
    
    max_height = 100
    
    sys.stdout.write("\t".join([location, str(p), str(max_height), str(width)])) 
    
    #
    #  Get the gene annotation
    #
    
    closest_tss = refseq_tss.get_tss_from_position(chr, \
      start, \
      end)
  
    tss = closest_tss.keys()[0]
    sys.stdout.write("\t" + str(tss) + "\t" + ";".join(sorted(closest_tss[tss])))
    # For testing with a single line
    #break
      
    sys.stdout.write("\n")
  
  f.close() 
  

file = sys.argv[1]

refseq_file = "/ifs/scratch/cancer/Lab_RDF/abh2138/references/ucsc/ucsc_refseq_exons_entrez_hg19.txt"
refseq_tss = annotation.RefSeqTss(refseq_file)

#
# Process the peaks
#

closest(file)
  
  # For testing one file
  #break
