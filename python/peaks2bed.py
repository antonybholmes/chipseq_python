# -*- coding: utf-8 -*-
"""
Convert a peaks file into a bed file

@author: Antony Holmes
"""

import sys
import re

import lib.bed
import lib.text
import lib.headings

def create_bed(file, name):
  sys.stdout.write(lib.bed.create_bed_header(name) + "\n")
  
  f = open(file, 'r')
  
  # skip header
  header = lib.text.get_header(f)
  
  gene_column = lib.text.find_index(header, lib.headings.GENE_SYMBOL)
  
  c = 1
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    location = lib.genomic.parse_location(tokens[0])
    
    genes = tokens[gene_column].split(";")
    
    for gene in sorted(genes):
      id = gene + ",peak_" + str(c)
    
      lib.bed.write_bed_line(location.chr, location.start, location.end, id)
    
      c += 1
    
    
  f.close()

file = sys.argv[1]

if len(sys.argv) > 2:
  name = sys.argv[2]
else:
  name = re.sub(r'\.txt','', file)

create_bed(file, name)
