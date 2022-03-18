"""
Extracts distributions of peak p-values to see how they
distribute
"""

import sys
import collections
import random
import math

import lib.genomic

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

MIN = -1000
PERMUTATIONS = 100000
MIN_D = 0.00001


def min_d(c, q):
  m = sys.float_info.max
  
  for v in c:
    d = v - q
    d *= d
    
    if d < m:
      m = d
  
  return m
  
  
def score(c1, c2, q):
  d1 = min_d(c1, q)
  d2 = min_d(c2, q)
  
  if d1 >= MIN_D:
    s1 = 1 / d1
  else:
    s1 = 0
    
  if d2 >= MIN_D:
    s2 = 1 / d2
  else:
    s2 = 0
  
  return s1 - s2
  
def score2(c1, c2, q):
  d1 = min_d(c1, q)
  d2 = min_d(c2, q)
  
  if d1 >= MIN_D:
    s1 = 1 / d1
  else:
    s1 = 0
    
  if d2 >= MIN_D:
    s2 = 1 / d2
  else:
    s2 = 0
    
  #sys.stderr.write("d1 " + str(d1) + " d2 " + str(d2) + "\n")
  
  return s1 - s2


def results(c1, c2, name_1, name_2):
  
  c1 = sorted(c1)
  c2 = sorted(c2)
  
  all_p = []
  
  all_p.extend(c1)
  all_p.extend(c2)
  
  indices = np.random.random_integers(0, len(all_p) - 1, PERMUTATIONS)
  
  score_dist = []
  
  c1_q25 = np.percentile(c1, 25)
  c1_q75 = np.percentile(c1, 75)
  c1_iqr = c1_q75 - c1_q25
  c1_min = c1_q25 - 1.5 * c1_iqr
  c1_max = c1_q75 + 1.5 * c1_iqr
  
  c2_q25 = np.percentile(c2, 25)
  c2_q75 = np.percentile(c2, 75)
  c2_iqr = c2_q75 - c2_q25
  c2_min = c2_q25 - 1.5 * c2_iqr
  c2_max = c2_q75 + 1.5 * c2_iqr
  
  sys.stderr.write("min max " + str(c1_min) + " " + str(c1_max) + "\n")
  sys.stderr.write("min max " + str(c2_min) + " " + str(c2_max) + "\n")
  
  
  for i in indices:
    q = all_p[i]
    score_dist.append(score(c1, c2, q))
  
  #for i in range(0, PERMUTATIONS):
  #  b = np.random.randint(0, len(all_p) - 1)
  #  q = all_p[b]
  #  #q = -(np.random.random() * 100)
  #  score_dist.append(score(c1, c2, q))
  
  score_dist = sorted(score_dist)
  
  t_5 = np.percentile(score_dist, 5)
  t_95 = np.percentile(score_dist, 95)
  
  c1_q = 0
  c2_q = -1000
  
  s1 = 0
  s2 = 0
  
  # Work out the cut offs for calling
  
  
  for q in -np.arange(0, 100, 0.1):
    s = score2(u_p_map_1, u_p_map_2, q)
    
    #sys.stderr.write("q " + str(q) + " s " + str(s) + "\n")
    
    # s1 (most negative score) must be the greatest value less
    # than t5
    if s > t_95:
      if q >= c1_min and q < c1_q:
        s1 = s
        c1_q = q
    
    # s2 (most postive score) must be the smalles value greater
    # than t95
    
    if s < t_5:
      if q <= c2_max and q > c2_q:
        s2 = s
        c2_q = q
        #break
        
  q1 = max(c1_q, c2_q)
  q2 = min(c1_q, c2_q)
      
  p1 = math.pow(10, q1)
  p2 = math.pow(10, q2)
  
  # See how many make the cut
  
  t1 = 0
  
  for p in c1:
    if p >= q1:
      t1 += 1
  
  t2 = 0
  
  for p in c2:
    if p <= q2:
      t2 += 1
  
  sys.stdout.write(name_1 + " vs " + name_2 + "\t" + "\t".join([str(t_5), str(t_95)]) + "\n")
  sys.stdout.write(name_1 + " score\t" + "\t".join([str(s1), str(q1), str(p1), str(t1)]) + "\n")
  sys.stdout.write(name_2 + " score\t" + "\t".join([str(s2), str(q2), str(p2), str(t2)]) + "\n")

#
# Begin
#

union_file = sys.argv[1]
file_1 = sys.argv[2]
file_2 = sys.argv[3]
name_1 = sys.argv[4]
name_2 = sys.argv[5]


#
# File 1
#

p_map_1 = collections.defaultdict(float)

f = open(file_1, 'r')

f.readline()

for line in f:
  tokens = line.strip().split("\t")
  
  l = tokens[0] + ":" + tokens[1] + "-" + tokens[2]
  
  p = float(tokens[3])
  
  if p > MIN:
    p_map_1[l] = p
  
  
f.close()

#
# File 2
#

p_map_2 = collections.defaultdict(float)

f = open(file_2, 'r')

f.readline()

for line in f:
  tokens = line.strip().split("\t")
  
  l = tokens[0] + ":" + tokens[1] + "-" + tokens[2]
  
  p = float(tokens[3])
  
  if p > MIN:
    p_map_2[l] = p
  
f.close()


u_p_map_1 = [] #collections.defaultdict(float)
u_p_map_2 = [] #collections.defaultdict(float)
u_p_map_3 = [] #collections.defaultdict(float)

f = open(union_file, 'r')

f.readline()

for line in f:
  tokens = line.strip().split("\t")
  
  l = tokens[0]
  
  loc = lib.genomic.parse_location(l)
  
  p = float(tokens[2])
  s1 = tokens[5]
  s2 = tokens[6]
  
  if s1 != "n/a" and s2 != "n/a":
    if p > MIN:
      u_p_map_2.append(p)
    
  if s1 != "n/a" and s2 == "n/a" and s1 in p_map_1:
    u_p_map_1.append(p_map_1[s1])
    
  if s1 == "n/a" and s2 != "n/a" and s2 in p_map_2:
    u_p_map_3.append(p_map_2[s2])
  
f.close()

results(u_p_map_1, u_p_map_2, name_1, "overlap")
results(u_p_map_3, u_p_map_2, name_2, "overlap")
