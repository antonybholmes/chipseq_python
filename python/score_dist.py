"""
Extracts distributions of peak p-values to see how they
distribute
"""

import sys
import collections

import lib.genomic

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

MIN = -1000

union_file = sys.argv[1]
name_1 = sys.argv[2]
name_2 = sys.argv[3]


s_1 = []
s_only_1 = []
s_overlap = []
s_only_2 = []
s_2 = []


f = open(union_file, 'r')

f.readline()

c = 0

for line in f:
  tokens = line.strip().split("\t")
  
  l = lib.genomic.parse_location(tokens[0])
  
  overlaps = int(tokens[4])
  
  id_1 = tokens[5]
  id_2 = tokens[6]
  
  type = tokens[15]
  
  s = float(tokens[3])
  
  if overlaps == 2:
    s_overlap.append(s)
  else:
    if id_1 != "n/a":
      s_only_1.append(s)
    else:
      s_only_2.append(s)
      
  if id_1 != "n/a":
    s_1.append(s)
    
  if id_2 != "n/a":
    s_2.append(s)
      
  c += 1
  
f.close()


#
# Plot w values
# 

data = []

data.append(s_1)
data.append(s_only_1)
data.append(s_overlap)
data.append(s_2)
data.append(s_only_2)

plt.figure(1)
plt.boxplot(data)
plt.xticks([1, 2, 3, 4, 5], [name_1, name_1 + " only", "Overlap", name_2, name_2 + " only"])
plt.ylabel("Score")
#plt.show()
plt.savefig(name_1 + "_" + name_2 + "_score_dist.pdf", format="pdf")
