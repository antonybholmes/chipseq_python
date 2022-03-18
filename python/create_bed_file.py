# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 11:35:50 2014

@author: Antony Holmes
"""

# creates a bed from from chipseeqer split sam files

import collections
import sys
import os
import re

def create_bed(window_size, color, chrom_size_file):

  pwd = os.getcwd()
  
  matcher = re.match(r'.*/(\w+_\w+_[A-Z]{2}\d+)/.*', pwd)

  name = matcher.group(1)
  
  sys.stderr.write("name: " + name + "\n")
  
  sizes = collections.defaultdict(int)
  
  f = open(chrom_size_file, 'r')
  
  # skip header
  f.readline()
  
  for line in f:
    line = line.strip()
    
    tokens = line.split("\t")
    
    chr = tokens[0]
    
    size = int(tokens[1])
    
    sizes[chr] = size
    
  f.close()

  positions = collections.defaultdict(lambda: collections.defaultdict(int))

  for chr in sizes:
    file = "reads." + chr
    
    if not os.path.exists(file):
      continue
    
    sys.stderr.write("Reading " + file + "...\n")
    
    f = open(file, 'r')
    
    for line in f:
      line = line.strip()
    
      tokens = line.split("\t")

      position = int(tokens[3])
    
      sequence = tokens[9]
    
      l = len(sequence)
    
      win = position // window_size * window_size
      end_win = (position + l - 1) // window_size * window_size
    
      while win <= end_win:
        positions[chr][win] += 1	
    
        win += window_size
    
    f.close()
  
  file = name + "_reads_bin_" + str(window_size) + ".bed"
  
  sys.stderr.write("Writing to " + file + "...\n");
  
  f = open(name + "_reads_bin_" + str(window_size) + ".bed", 'w')

  f.write("track type=bedGraph name=\"" + name + "\" description=\"" + name + " reads\" autoScale=on alwaysZero=on visibility=full color="  + color + "\n")

  for chr in sorted(positions):
    for start in sorted(positions[chr]):
      end = min(sizes[chr] - 1, start + window_size - 1)
      
      cpm = positions[chr][start]
      
      f.write("\t".join([chr, str(start), str(end), str(cpm)]) + "\n")
  
  f.close()

window_size = int(sys.argv[1])
color = sys.argv[2]
chrom_size_file = sys.argv[3]

create_bed(window_size, color, chrom_size_file)