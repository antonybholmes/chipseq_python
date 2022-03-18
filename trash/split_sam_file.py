# -*- coding: utf-8 -*-
"""
Split reads by chromosome. Reads go to separate files.

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import collections
import re
import sys

def split(file):
  file_map = collections.defaultdict()

  sys.stderr.write("Parsing " + file + "...\n")

  f = open(file, 'r')
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
  
    if re.match(r'^@.*', line):
      continue
  
    tokens = line.split("\t")
  
    chr = tokens[2]
    
    if re.match(r'.*chrUn.*', chr):
      continue

    if re.match(r'.*random.*', chr):
      continue
    
    if chr not in file_map:
      fout = open("reads." +  chr, 'w')
      
      file_map[chr] = fout
    
    fout = file_map[chr]

    fout.write(line + "\n")
  
  f.close()
  
  for chr in file_map:
    file_map[chr].close()


file = sys.argv[1]

split(file)
