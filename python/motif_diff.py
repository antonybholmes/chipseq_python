"""
Look for differences between two motif lists
"""

import sys
import collections

import lib.genomic

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


diff_file = sys.argv[1]
min_q = float(sys.argv[2])
min_sens = float(sys.argv[3])
min_spec = float(sys.argv[4])
file = sys.argv[5]

min_delta = 4

special_names = set()

special_names.add("mef2")
special_names.add("mef-2")

spec_sens = set()

#
# File 1
#

ids = set()
qs = collections.defaultdict(float)

f = open(diff_file, 'r')

f.readline()

for line in f:
  tokens = line.strip().split("\t")
  
  id = tokens[1].lower()
  
  q = float(tokens[4])
  
  sens = float(tokens[9])
  spec = float(tokens[10])
  
  if "pb0035.1" in id:
    sys.stderr.write(id + "\n")
    sys.stderr.write(str(q) + " " + str(sens) + " " + str(spec) + "\n")
  
  if q >= min_q and sens >= min_sens and spec >= min_spec:
    #sys.stderr.write("add " + id + "\n")
    ids.add(id)
    
    qs[id] = q
  else:
    spec_sens.add(id)
  
f.close()


f = open(file, 'r')

sys.stdout.write(f.readline().strip() + "\tType\tDelta Q\tQ From Other Group\n")

for line in f:
  line = line.strip()
  
  tokens = line.split("\t")
  
  name = tokens[0].lower()
  
  id = tokens[1].lower()
  
  q = float(tokens[4])
  sens = float(tokens[9])
  spec = float(tokens[10])
  
  found = False
  
  for n in special_names:
    if n in name:
      found = True
      break
  
  if found and id not in ids:
    sys.stdout.write(line + "\tmef2\tn/a\tn/a\n")
  elif q >= min_q and sens >= min_sens and spec >= min_spec:
    if id not in ids:
      if id in spec_sens:
        sys.stdout.write(line + "\tspec_sens\tn/a\tn/a\n")
      else:
        sys.stdout.write(line + "\tno_id\tn/a\tn/a\n")
    else:
      q2 = qs[id]
      delta = abs(q - q2)
      
      if delta >= min_delta:
        sys.stdout.write(line + "\tdelta\t" + str(delta) + "\t" + str(q2) + "\n")
  else:
    pass
        
f.close()
