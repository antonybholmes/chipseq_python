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
c1 = float(sys.argv[6])
c2 = float(sys.argv[7])
c3 = float(sys.argv[8])
c4 = float(sys.argv[9])
y1 = int(sys.argv[10])
y2 = int(sys.argv[11])
y3 = int(sys.argv[12])

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



# distribution

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.set_dpi(300)
fig.set_size_inches(8, 6)

ax1.hist(data[1], 50, range=[-100, 0], alpha = 0.5, label=name_1 + " only", color='#0055d4')
ax1.hist(data[2], 50, range=[-100, 0], alpha = 0.5, label='Overlap', color='#ff0000')
ax1.set_ylim(y2, y3)
ax1.plot((c1, c1), (0, ax1.get_ylim()[1]), 'k--')
ax1.plot((c2, c2), (0, ax1.get_ylim()[1]), 'k--')
ax1.legend(loc='upper left')

ax2.hist(data[1], 50, range=[-100, 0], alpha = 0.5, label=name_1 + " only", color='#0055d4')
ax2.hist(data[2], 50, range=[-100, 0], alpha = 0.5, label='Overlap', color='#ff0000')
ax2.set_ylim(0, y1)
ax2.plot((c1, c1), (0, ax2.get_ylim()[1]), 'k--')
ax2.plot((c2, c2), (0, ax2.get_ylim()[1]), 'k--')


ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax1.xaxis.tick_top()
ax1.tick_params(labeltop='off')  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

d = .015
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal


plt.xlabel("log10 p-value")
plt.ylabel("Peak count")
#plt.title(name_1 + " only compared to Overlap")

plt.savefig(name_1 + "_overlap_pvalue_dist_hist.pdf", format="pdf")
