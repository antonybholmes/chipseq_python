# -*- coding: utf-8 -*-
"""
Annotate chipseeqer peaks

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re
import gzip
import math


def base_qc(file, l):
  n = 0
  mean = 0.0
  m2 = 0.0
  
  base_n = [0 for i in range(l)]
  base_mean = [0.0 for i in range(l)]
  base_m2 = [0.0 for i in range(l)]
  
  indices = range(0, l)
  
  f = gzip.open(file, 'r')
  
  ml = 5
  
  lc = 0
  
  for line in f:
    lc += 1
    
    # skip everything but the 4th line
    if not lc % 4 == 0:
      continue
    
    line = line.strip()
 
    for i in indices:
      c = line[i]
      v = ord(c) - 33
      
      n += 1
      delta = v - mean
      mean += delta / n
      m2 += delta * (v - mean)
      
      base_n[i] += 1
      delta = v - base_mean[i]
      base_mean[i] += delta / base_n[i]
      base_m2[i] += delta * (v - base_mean[i])
    
    #sys.stderr.write(line + "\n")
  
  f.close()
  
  m2 /= (n - 1)
  sys.stderr.write(str(mean) + " " + str(m2) + " " + str(math.sqrt(m2)) + "\n")
  
  for i in indices:
    base_m2[i] /= (base_n[i] - 1)
    sys.stderr.write(str(i + 1) + " " + str(base_mean[i]) + " " + str(base_m2[i]) + " " + str(math.sqrt(base_m2[i])) + "\n")

file = sys.argv[1]
l = int(sys.argv[2])

base_qc(file, l)
