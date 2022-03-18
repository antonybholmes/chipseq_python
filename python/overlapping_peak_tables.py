# -*- coding: utf-8 -*-
"""
Find the common overlaps between a list of peaks in peak table files

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re
import collections

import lib.peaks
import lib.genomic

files = sys.argv[1:len(sys.argv)] #sorted(sys.argv[1:len(sys.argv)])
ids = []

pvalues = collections.defaultdict(float)
scores = collections.defaultdict(float)


for file in files:
  sys.stderr.write("file: " + file + "\n")
  
  matcher = re.match(r'.*Peaks_(.+?)_vs.*', file)
  
  id = matcher.group(1)
    
  ids.append(id)
  
  sys.stderr.write("id " + id + "\n")

  # now make a list of locations and best p-values  

  f = open(file, 'r')
  
  f.readline()
    
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    location = lib.genomic.parse_location(tokens[0])
    
    #sys.stderr.write(location.to_string() + "\n")
    
    p = float(tokens[1])
    score = float(tokens[2])
    
    lid = id + "=" + location.chr + ":" + str(location.start) + "-" + str(location.end)
    
    pvalues[lid] = p
    scores[lid] = score
    
    #sys.stderr.write(location + " " + str(p) + "\n")
    
  f.close()

  
location_core_map = lib.peaks.overlapping_peak_tables(files, ids)

# keep the ids sorted
#ids = sorted(ids)

sys.stdout.write("Genomic Location (hg19)\tPeak / Overlap Width\tBest P-value (ChIPseeqer)\tBest Score (ChIPseeqer)\tNumber Of Overlapping Peaks\t" + "\t".join(["Sample " + s for s in ids]) + "\n");

for core_location in sorted(location_core_map):
  location = lib.genomic.parse_location(core_location)
  
  overlap = location.end - location.start + 1
  
  c = 0
  
  for id in location_core_map[core_location]:
    #for location in location_core_map[core_location][id]:
    c += 1
      
  p = 0

  score = 0  
  
  for id in ids:
    if id in location_core_map[core_location]:
      #for location in location_core_map[core_location][id]:
      location = location_core_map[core_location][id]
      
      if pvalues[location] < p:
        p = pvalues[location]

      if scores[location] > score:
        score = scores[location]
  
  sys.stdout.write("\t".join([core_location, str(overlap), str(p), str(score), str(c)]));

  # foreach RK id
  for id in ids:
    if id in location_core_map[core_location]:
      peaks = location_core_map[core_location][id] #";".join(sorted(location_core_map[core_location][id]))
      
      # get rid of internal id
      #lib.peaks = re.sub(r'^.+=', "", lib.peaks)
      
      sys.stdout.write("\t" + peaks)
    else:
      sys.stdout.write("\tn/a")
  
  sys.stdout.write("\n")
