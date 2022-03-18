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
file_1 = sys.argv[2]
file_2 = sys.argv[3]
name_1 = sys.argv[4]
name_2 = sys.argv[5]


#
# File 1
#

p_map_1 = collections.defaultdict(float)
w_map_1 = collections.defaultdict(int)

f = open(file_1, 'r')

f.readline()

for line in f:
  tokens = line.strip().split("\t")
  
  l = tokens[0] + ":" + tokens[1] + "-" + tokens[2]
  
  w = int(tokens[2]) - int(tokens[1])
  w_map_1[l] = w
  
  p = float(tokens[3])
  
  if p > MIN:
    p_map_1[l] = p
  
  
f.close()

#
# File 2
#

p_map_2 = collections.defaultdict(float)
w_map_2 = collections.defaultdict(int)

f = open(file_2, 'r')

f.readline()

for line in f:
  tokens = line.strip().split("\t")
  
  l = tokens[0] + ":" + tokens[1] + "-" + tokens[2]
  
  w = int(tokens[2]) - int(tokens[1])
  w_map_2[l] = w
  
  p = float(tokens[3])
  
  if p > MIN:
    p_map_2[l] = p
  
f.close()


u_p_map_1 = collections.defaultdict(float)
u_p_map_2 = collections.defaultdict(float)
u_p_map_3 = collections.defaultdict(float)

u_w_map_1 = collections.defaultdict(int)
u_w_map_2 = collections.defaultdict(int)
u_w_map_3 = collections.defaultdict(int)


f = open(union_file, 'r')

f.readline()

for line in f:
  tokens = line.strip().split("\t")
  
  l = tokens[0]
  
  loc = lib.genomic.parse_location(l)
  w = loc.end - loc.start
  p = float(tokens[2])
  s1 = tokens[5]
  s2 = tokens[6]
  
  if s1 != "n/a" and s2 != "n/a":
    u_w_map_2[l] = w
    
    if p > MIN:
      u_p_map_2[l] = p
    
  if s1 != "n/a" and s2 == "n/a" and s1 in p_map_1:
    
    u_p_map_1[s1] = p_map_1[s1]
    u_w_map_1[s1] = w_map_1[s1]
    
  if s1 == "n/a" and s2 != "n/a" and s2 in p_map_2:
    u_p_map_3[s2] = p_map_2[s2]
    u_w_map_3[s2] = w_map_2[s2]
  
f.close()

#
# Plot w values
# 

data = []

d = []

for l in w_map_1:
  d.append(w_map_1[l])

data.append(d)

d = []

for l in u_w_map_1:
  d.append(u_w_map_1[l])

data.append(d)

d = []

for l in u_w_map_2:
  d.append(u_w_map_2[l])

data.append(d)

d = []

for l in w_map_2:
  d.append(w_map_2[l])

data.append(d)

d = []

for l in u_w_map_3:
  d.append(u_w_map_3[l])

data.append(d)


plt.figure(1)
plt.boxplot(data)
plt.xticks([1, 2, 3, 4, 5], [name_1, name_1 + " only", "Overlap", name_2, name_2 + " only"])
plt.ylabel("Width")
#plt.show()
plt.savefig("width_dist.pdf", format="pdf")


#
# Plot p values
# 

data = []

d = []

for l in p_map_1:
  
  d.append(p_map_1[l])

data.append(d)

d = []

for l in u_p_map_1:
  d.append(u_p_map_1[l])

data.append(d)

d = []

for l in u_p_map_2:
  d.append(u_p_map_2[l])

data.append(d)

d = []

for l in p_map_2:
  d.append(p_map_2[l])

data.append(d)

d = []

for l in u_p_map_3:
  d.append(u_p_map_3[l])

data.append(d)

sys.stderr.write(str(max(max(data))) + "\n")

# widths p values
fig = plt.figure(2)
axes = fig.add_subplot(111)
bp = axes.boxplot(data)
axis = axes.xaxis
axis.set_ticks([1, 2, 3, 4, 5])
axis.set_ticklabels([name_1, name_1 + " only", "Overlap", name_2, name_2 + " only"])
axes.yaxis.set_label_text("log10 p-value")
#plt.show()
fig.savefig(name_1 + "_" + name_2 + "_peak_dist.pdf", format="pdf")


# distribution

fig = plt.figure(3, figsize=(20, 16), dpi=300)

ax = fig.add_subplot(231)
ax.hist(data[0], 50, range=[-100, 0], alpha = 0.5, label=name_1)
plt.xlabel("log10 p-value")
plt.ylabel("Count")
plt.title(name_1)

ax = fig.add_subplot(232)
ax.hist(data[1], 50, range=[-100, 0], alpha = 0.5, label=name_1 + " only")
plt.xlabel("log10 p-value")
plt.ylabel("Count")
plt.title(name_1 + " only")

ax = fig.add_subplot(233)
ax.hist(data[2], 50, range=[-100, 0], alpha = 0.5, label='Overlap')
plt.xlabel("log10 p-value")
plt.ylabel("Count")
plt.title("Overlap")

ax = fig.add_subplot(234)
ax.hist(data[3], 50, range=[-100, 0], alpha = 0.5, label=name_2)
plt.xlabel("log10 p-value")
plt.ylabel("Count")
plt.title(name_2)

ax = fig.add_subplot(235)
ax.hist(data[4], 50, range=[-100, 0], alpha = 0.5, label=name_2 + " only")
plt.xlabel("log10 p-value")
plt.ylabel("Count")
plt.title(name_2 + " only")

#plt.legend(loc='upper left')
#plt.show()
plt.savefig(name_1 + "_" + name_2 + "_peak_dist_hist.pdf", format="pdf")



f = open("d1.txt", "w")

for l in u_p_map_1:
  f.write(str(u_p_map_1[l]) + "\n")
  
f.close()

f = open("d2.txt", "w")

for l in u_p_map_2:
  f.write(str(u_p_map_2[l]) + "\n")
  
f.close()

f = open("d3.txt", "w")

for l in u_p_map_3:
  f.write(str(u_p_map_3[l]) + "\n")
  
f.close()


t = stats.ttest_ind(data[0], data[1], equal_var = False)
sys.stderr.write(name_1 + " vs overlap " + str(t[1]) + "\n")

t = stats.ttest_ind(data[1], data[2], equal_var = False)
sys.stderr.write(name_2 + " vs overlap " + str(t[1]) + "\n")

t = stats.ttest_ind(data[0], data[2], equal_var = False)
sys.stderr.write(name_1 + " vs " + name_2 + " " + str(t[1]) + "\n")
