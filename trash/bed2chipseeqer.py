# -*- coding: utf-8 -*-
"""
Convert a peaks file into a chipseeqer file

@author: Antony Holmes
"""

import sys
import re

import lib.bed
import lib.text
import lib.headings

def create(file):
  f = open(file, 'r')
  
  hasHeader = f.readline().startswith("track")
  
  f.close()
  
  f = open(file, 'r')
  
  if hasHeader:
    header = lib.text.get_header(f)
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    chr = tokens[0]
    start = tokens[1]
    end = tokens[2]
    
    sys.stdout.write("\t".join([chr, start, end, "-1000", "100"]) + "\n")
    
  f.close()

file = sys.argv[1]

create(file)
