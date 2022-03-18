# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 18:16:57 2015

@author: antony
"""

import lib_dna
import re
import sys

def score(file, genomic_column):
  f = open(file, 'r')
  
  line = f.readline().strip() + "\tRepetitive Mask P"

  sys.stdout.write(line + "\n")  
  
  for line in f:
    line = line.strip()
    
    tokens = line.split("\t")
    
    location = tokens[genomic_column]
    
    matcher = re.match(r'(chr.+?):(\d+)-(\d+)', location)
    
    chr = matcher.group(1)
    start = int(matcher.group(2))
    end = int(matcher.group(3))
    
    sequence = lib_dna.get_sequence_with_mask(chr, start, end)

    score = lib_dna.score_mask_seq(sequence) 
    
    sys.stdout.write(line + "\t" + "{:.2f}".format(score) + "\n")
    
  f.close()
  
file = sys.argv[1]
genomic_column = int(sys.argv[2])

score(file, genomic_column)