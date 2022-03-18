# -*- coding: utf-8 -*-
"""
Find the common overlaps between a list of files reported as a common region

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import peaks
import re
import collections

import lib_genomic

files = sys.argv[1:len(sys.argv)] #sorted(sys.argv[1:len(sys.argv)])
ids = []

pvalues = collections.defaultdict(float)
scores = collections.defaultdict(float)


for file in files:
  sys.stderr.write("file: " + file + "\n")
  
  matcher = re.match(r'.*TF_targets_(.+)_vs.*', file)
    
  id = matcher.group(1)
  
  ids.append(id)
  
  sys.stderr.write("ids " + id + "\n")

  # now make a list of locations and best p-values  

  f = open(file, 'r')
    
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    chr = tokens[0]
    start = int(tokens[1])
    end = int(tokens[2])
    
    # account for inf values
    if re.match(r'.*inf.*', tokens[3]):    
      p = -1000
    else:      
      p = float(tokens[3])
      
    score = float(tokens[4])
    
    location = id + "=" + chr + ":" + str(start) + "-" + str(end)
    
    pvalues[location] = p
    scores[location] = score
    
    #sys.stderr.write(location + " " + str(p) + "\n")
    
  f.close()

  
location_core_map = peaks.overlapping_v2(files, ids)

# keep the ids sorted
#ids = sorted(ids)

sys.stdout.write("Genomic Location (hg19)\tPeak / Overlap Width\tBest P-value (ChIPseeqer)\tBest Score (ChIPseeqer)\tNumber Of Overlapping Peaks\t" + "\t".join(["Sample " + s for s in ids]) + "\n");

for core_location in sorted(location_core_map):
  location = lib_genomic.parse_location(core_location)
  
  overlap = location.end - location.start + 1
  
  c = 0
  
  for id in location_core_map[core_location]:
    for location in location_core_map[core_location][id]:
      c += 1
      
  p = 0

  score = 0  
  
  for id in ids:
    if id in location_core_map[core_location]:
      for location in location_core_map[core_location][id]:
        if pvalues[location] < p:
          p = pvalues[location]

        if scores[location] > score:
          score = scores[location]
  
  sys.stdout.write("\t".join([core_location, str(overlap), str(p), str(score), str(c)]));

  # foreach RK id
  for id in ids:
    if id in location_core_map[core_location]:
      peaks = ";".join(sorted(location_core_map[core_location][id]))
      
      # get rid of internal id
      #peaks = re.sub(r'^.+=', "", peaks)
      
      sys.stdout.write("\t" + peaks)
    else:
      sys.stdout.write("\tn/a")
  
  sys.stdout.write("\n")
