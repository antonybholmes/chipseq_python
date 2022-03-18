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

def closest(file):
  new_file = re.sub(r'\.txt', "_refseq_hg19.txt", file)
  new_file = re.sub(r'TF_targets_', "", new_file)
  new_file = "closest_" + new_file
    
  sys.stderr.write("Creating Closest from " + file + " as " + new_file + "\n")  

  bins = collections.defaultdict(int)

  f = open(file, "r")
  fout = open(new_file,'w')
  
  fout.write("Genomic Location (hg19)")
  fout.write("\tP-value (ChIPseeqer)");
  fout.write("\tMax Height (ChIPseeqer)");
  fout.write("\tPeak Size");
  fout.write("\tTSS Distance");
  fout.write("\tTSS Entrez ID");
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
    
    if re.match(r'.*[Ii]nf.*', tokens[3]):
      p = -1000
    else:
      p = float(tokens[3])
    
    max_height = int(tokens[6])
    
    location = chr + ":" + str(start) + "-" + str(end)

    fout.write("\t".join([location, str(p), str(max_height), str(width)])) 
    
    #
    #  Get the gene annotation
    #
    
    closest_tss = refseq_tss.get_tss_from_position(chr, \
      start, \
      end)
  
    tss =  closest_tss.keys()[0]
    fout.write("\t" + str(tss) + "\t" + ";".join(sorted(closest_tss[tss])))
    # For testing with a single line
    #break
      
    fout.write("\n")
  
  f.close()
  
  fout.close() 
  

dir = "."

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
  
  closest(file)
  
  # For testing one file
  #break
