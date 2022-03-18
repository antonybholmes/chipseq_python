"""
Given peaks from a group, extract those that we have labeled
as not being classifiable into either the overlap or unique to
one phenotype
"""

import sys
import collections

import lib.genomic

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

MIN = -1000

file = sys.argv[1]
type = int(sys.argv[2])
min_p = float(sys.argv[3])
max_p = float(sys.argv[4])


f = open(file, 'r')

sys.stdout.write(f.readline().strip() + "\n")

for line in f:
  line = line.strip()
  
  tokens = line.split("\t")
  
  l = tokens[0]
  
  loc = lib.genomic.parse_location(l)
  w = loc.end - loc.start
  p = float(tokens[2])
  s1 = tokens[5]
  s2 = tokens[6]
  
  if p >= min_p and p <= max_p:
    if type == 1 and s1 != "n/a" and s2 == "n/a":
      sys.stdout.write(line + "\n")

    if type == 2 and s1 != "n/a" and s2 != "n/a":
      sys.stdout.write(line + "\n")
    
    if type == 3 and s1 == "n/a" and s2 != "n/a":
      sys.stdout.write(line + "\n")
  
f.close()
