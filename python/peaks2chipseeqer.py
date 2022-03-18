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
  
  # skip header
  header = lib.text.get_header(f)
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    location = lib.genomic.parse_location(tokens[0])
    
    sys.stdout.write("\t".join([location.chr, str(location.start), str(location.end), "-1000", "100"]) + "\n")
    
  f.close()

file = sys.argv[1]

create(file)
