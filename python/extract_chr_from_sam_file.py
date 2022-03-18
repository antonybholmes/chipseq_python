# -*- coding: utf-8 -*-
"""
Split reads by chromosome. Reads go to separate files.

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import re
import sys

def extract_chr(file, chrom):
  sys.stderr.write("Parsing " + file + " for " + chrom + "...\n")
  
  fout = open("reads." + chrom, 'w')
  
  f = open(file, 'r')
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
  
    if re.match(r'^@.*', line):
      continue
  
    tokens = line.split("\t")
  
    reference = tokens[2]
    
    if reference != chrom:
      continue
    
    fout.write(line + "\n")
  
  f.close()
  
  fout.close()


file = sys.argv[1]
chrom = sys.argv[2]

extract_chr(file, chrom)
