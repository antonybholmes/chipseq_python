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


d_1 = []
d_only_1 = []
d_overlap = []
d_only_2 = []
d_2 = []


p_1 = 0
p_only_1 = 0
p_overlap = 0
p_only_2 = 0
p_2 = 0

# intronic
i_1 = 0
i_only_1 = 0
i_overlap = 0
i_only_2 = 0
i_2 = 0

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
  
  d = abs(int(tokens[16]))
  
  if overlaps == 2:
    d_overlap.append(d)
    
    if "promoter" in type:
      p_overlap += 1
    
    if "intronic" in type:
      i_overlap += 1
  else:
    if id_1 != "n/a":
      d_only_1.append(d)
      
      if "promoter" in type:
        p_only_1 += 1
      
      if "intronic" in type:
        i_only_1 += 1
    else:
      d_only_2.append(d)
      
      if "promoter" in type:
        p_only_2 += 1
        
      if "intronic" in type:
        i_only_2 += 1
      
  if id_1 != "n/a":
    d_1.append(d)
    
    if "promoter" in type:
      p_1 += 1
      
    if "intronic" in type:
      i_1 += 1
    
  if id_2 != "n/a":
    d_2.append(d)
    
    if "promoter" in type:
      p_2 += 1
    
    if "intronic" in type:
      i_2 += 1
      
  c += 1
  
f.close()


#
# Plot w values
# 

data = []

data.append(d_1)
data.append(d_only_1)
data.append(d_overlap)
data.append(d_2)
data.append(d_only_2)

plt.figure(1)
plt.boxplot(data)
plt.xticks([1, 2, 3, 4, 5], [name_1, name_1 + " only", "Overlap", name_2, name_2 + " only"])
plt.ylabel("Distance")
#plt.show()
plt.savefig(name_1 + "_" + name_2 + "_promoter_dist.pdf", format="pdf")

sys.stderr.write(name_1 + " p " + str(p_1) + "\n")
sys.stderr.write(name_1 + " only p " + str(p_only_1) + "\n")
sys.stderr.write("overlap p " + str(p_overlap) + "\n")
sys.stderr.write(name_2 + " p " + str(p_2) + "\n")
sys.stderr.write(name_2 + " only p " + str(p_only_1) + "\n")

sys.stderr.write(name_1 + " i " + str(i_1) + "\n")
sys.stderr.write(name_1 + " only i " + str(i_only_1) + "\n")
sys.stderr.write("overlap i " + str(i_overlap) + "\n")
sys.stderr.write(name_2 + " i " + str(i_2) + "\n")
sys.stderr.write(name_2 + " only i " + str(i_only_1) + "\n")


sys.stderr.write("c " + str(c) + "\n")
