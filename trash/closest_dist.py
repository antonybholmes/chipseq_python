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


file_1 = sys.argv[1]
name_1 = sys.argv[2]
output = sys.argv[3]


#
# File 1
#

data = []

f = open(file_1, 'r')

f.readline()

for line in f:
  tokens = line.strip().split("\t")
  
  l = tokens[0] + ":" + tokens[1] + "-" + tokens[2]
  
  d = abs(int(tokens[3]))
  
  if d <= 1000000:
    data.append(d)
  
  
f.close()


plt.figure(1)
plt.boxplot(data)
plt.xticks([1], [name_1])
plt.ylabel("Distance")
#plt.show()
plt.savefig(output, format="pdf")
