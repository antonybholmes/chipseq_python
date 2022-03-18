# -*- coding: utf-8 -*-
"""
Find the common overlaps between a list of files reported as a common region

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re
import collections

import lib.enhancers
import lib.genomic

files = sys.argv[1:len(sys.argv)] #sorted(sys.argv[1:len(sys.argv)])
ids = []

for file in files:
  sys.stderr.write("file: " + file + "\n")
  
  matcher = re.match(r'.*Peaks_(.+)_vs.*', file)
    
  id = matcher.group(1)
  
  ids.append(id)
  
  sys.stderr.write("ids " + id + "\n")

location_core_map = lib.enhancers.get_overlaps(files, ids)

# keep the ids sorted
#ids = sorted(ids)

sys.stdout.write("Genomic Location (hg19)\tEnhancer / Overlap Width\tNumber Of Overlapping Peaks\t" + "\t".join(["Sample " + s for s in ids]) + "\n");

for core_location in sorted(location_core_map):
  location = lib.genomic.parse_location(core_location)
  
  overlap = location.end - location.start + 1
  
  c = 0
  
  for id in location_core_map[core_location]:
    for location in location_core_map[core_location][id]:
      c += 1
  
  sys.stdout.write("\t".join([core_location, str(overlap), str(c)]));

  for id in ids:
    if id in location_core_map[core_location]:
      peaks = ";".join(sorted(location_core_map[core_location][id]))
      
      sys.stdout.write("\t" + peaks)
    else:
      sys.stdout.write("\tn/a")
  
  sys.stdout.write("\n")
