# -*- coding: utf-8 -*-
"""
Find the common overlaps between a list of files reported as a common region

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re
import collections

import lib.peaks
import lib.genomic
import lib.overlap

file = sys.argv[1]
histone_files = sys.argv[2:len(sys.argv)] #sorted(sys.argv[1:len(sys.argv)])

#
# Load all the histones
#
histones = collections.defaultdict(lib.overlap.Overlapping)
ids = []

for histone_file in histone_files:
  matcher = re.match(r'^.+\/(.+)\.', histone_file)
  
  id = matcher.group(1)
  
  overlap = lib.overlap.Overlapping(histone_file)
  
  histones[id] = overlap
  
  ids.append(id)

#
# Open the file
#

f = open(file, 'r')

header = f.readline().strip().split("\t")

# Add histone overlaps to header

for id in ids:
  header.append(id + " Overlap")
  
sys.stdout.write("\t".join(header) + "\n")

for line in f:
  line = line.strip()
  
  sys.stdout.write(line)

  tokens = line.split("\t")
    
  location = lib.genomic.parse_location(tokens[0])
  
  for id in ids:
    overlaps = histones[id].overlaps(location)
    
    sys.stdout.write("\t")
    
    if len(overlaps) > 0:
      sys.stdout.write(overlaps[0].to_string())
    else:
      sys.stdout.write("n/a")
       
  sys.stdout.write("\n")
  
f.close()
