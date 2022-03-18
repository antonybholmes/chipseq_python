# -*- coding: utf-8 -*-
"""
Encode read counts per base in 1 byte

@author: Antony Holmes
"""

import sys
import collections
import numpy


def encode_sam(file, chromosome, read_length, window):
  chr_sizes = collections.defaultdict(int)
  
  f = open('/ifs/scratch/cancer/Lab_RDF/abh2138/references/ucsc/hg19_chromosome_sizes.txt', 'r')

  f.readline()
  
  
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    chr = tokens[0]
    size = int(tokens[1])
    
    chr_sizes[chr] = size
  
  f.close()
  
  # size of this chromosome in windows
  size = (chr_sizes[chromosome] // window) + 1
  
  sys.stderr.write("Parsing " + file + " " + str(read_length) + " " + str(size) + " " + str(window) + "...\n")

  counts = numpy.zeros(size, int)
  max_count = 0
  max_pos = 0
  
  f = open(file, 'r')
  
  for line in f:
    line = line.strip()

    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    chr = tokens[2]
    
    if chr != chromosome:
      continue

    start = int(tokens[3]) - 1
    end = start + read_length - 1
    
    win_start = start // window
    win_end = end // window
    
    for i in xrange(win_start, win_end + 1):
      counts[i] += 1
      
      if counts[i] > max_count:
        max_count = counts[i]
        max_pos = i
      
    #break
    
  f.close()
  
  sys.stderr.write(str(max_pos) + " " + str(max_count) + "\n")
      

file = sys.argv[1]
chr = sys.argv[2]
read_length = int(sys.argv[3])
window = int(sys.argv[4])

encode_sam(file, chr, read_length, window)