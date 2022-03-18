# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 17:31:13 2014

@author: Antony Holmes
"""

# Extract moving genes in comparison tables so we can do pathway
# analysis

import sys
import collections
import re

type = sys.argv[1]
group = sys.argv[2]
file = sys.argv[3]

f = open(file, "r")

# skip header
f.readline()
  
for line in f:
  ls = line.strip()
    
  if len(ls) == 0:
    continue
    
  tokens = ls.split("\t")
  
  if tokens[2] != type:
    continue
  
  if (group == "group1" and tokens[7] == "n/a") or (group == "group2" and tokens[11] == "n/a"):
    continue
  
  gene = tokens[0]

  sys.stdout.write(gene + "\n")

f.close()