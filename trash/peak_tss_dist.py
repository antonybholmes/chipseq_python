# -*- coding: utf-8 -*-
"""
Generate a tss distribution for a chipseeqer peak file

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import collections
import re
import os

import annotation

def tss(file, bin_size):
  new_file = re.sub(r'\.txt', "_refseq_hg19_bin_" + str(bin_size) + ".txt", file)
  new_file = re.sub(r'TF_targets_', "", new_file)
  new_file = "tss_" + new_file
    
  sys.stderr.write("Creating TSS from " + file + " as " + new_file + "\n")  

  bins = collections.defaultdict(int)

  f = open(file, "r")
  fout = open(new_file,'w')
  
  #header
  fout.write("Distance (" + str(bin_size) + ")")
  fout.write("\t");
  fout.write("Peak count")
  fout.write("\n");
  
  total_peaks = 0
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t");

    chr = tokens[0]  
    start = int(tokens[1])
    end = int(tokens[2])
    width = end - start + 1

    location = chr + ":" + str(start) + "-" + str(end)

    #
    #  Get the gene annotation
    #
    
    closest_tss = refseq_tss.get_tss_from_position(chr, \
      start, \
      end)
  
    for tss in closest_tss:
      tss_bin = int(tss / bin_size)
      
      if tss_bin in bins:
        bins[tss_bin] += 1
      else:
        bins[tss_bin] = 1
    
    # For testing with a single line
    #break
  
  f.close()
  
  s = min(bins)
  e = max(bins)
  
  for bin in range(s, e + 1):
    p = bin * bin_size
    
    if bin in bins:
      fout.write(str(p) + "\t" + str(bins[bin]) + "\n")
    else:
      fout.write(str(p) + "\t0\n")
    
  fout.close() 
  

dir = "."


tss_bin_size = int(sys.argv[1]) #1000

refseq_file = "/ifs/scratch/cancer/Lab_RDF/abh2138/references/ucsc/ucsc_refseq_exons_entrez_hg19.txt"
refseq_tss = annotation.RefSeqTss(refseq_file)

#
# Process the peaks
#

files = os.listdir(dir)

for file in files:
  #file  = "TF_targets_RK053_CB5_CREBBP_IEO_vs_Input_IEO_p5.txt"
  
  if not re.match(r'TF_targets.*txt', file):
    continue
  
  tss(file, tss_bin_size)
  
  # For testing one file
  #break