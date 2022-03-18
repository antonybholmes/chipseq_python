# -*- coding: utf-8 -*-
"""
Create bins for a tss dist

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import collections
import math

def create_bins(file, bin_size):
  bins = collections.defaultdict(int)  

  # for binning we need the bin in bp so we can round the bin to the
  # nearest interval
 

  diff = bin_size / 100000.0
  
  f = open(file, 'r')
  
  f.readline()
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    p = float(tokens[0])
    
    p = abs(p)
    
    if p == 0:
      continue
    
    p = max(1, p - diff)
    
    #if p < 1:
    #  continue
    
    p = math.log10(p)
    
    bin = p / bin_size
    
    fbin = math.floor(bin)

    if fbin not in bins:
      bins[fbin] = 0
      
    bins[fbin] += 1
    
  f.close()
  
 
  sys.stdout.write("Bin (bp)\tCount\n")
  
  b = 0
  m = 80 #max(bins)
  
  while b <= m:
    #sys.stdout.write("{:.2f}".format(math.pow(10, b * bin_size)) + "\t")
    sys.stdout.write("{:.2E}".format(math.pow(10, b * bin_size)) + "\t")
    
    if b in bins:
      sys.stdout.write(str(bins[b]))
    else:
      sys.stdout.write("0")
    
    sys.stdout.write("\n")
    
    b += 1
    
    
file = sys.argv[1]
bin_size = float(sys.argv[2])

create_bins(file, bin_size)