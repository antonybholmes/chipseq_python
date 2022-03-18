# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 17:48:31 2014

@author: Antony Holmes
"""

import sys
import re

file = sys.argv[1]
out = re.sub(r'\.txt', '.bed', re.sub(r'TF_targets_', '', file))

f = open(file, 'r')
fout = open(out, 'w')
  
for line in f:
  tokens = line.strip().split("\t")
    
  chr = tokens[0]
  start = int(tokens[1])
  end = int(tokens[2])

  label = chr + ":" + str(start) + "-" + str(end)

  fout.write("\t".join([chr, str(start), str(end), label]) + "\n")

f.close()
fout.close()
