# -*- coding: utf-8 -*-
"""
Find the common overlaps between peak annotation files created by the 
pipeline. This is different to the overlapping of TF target files

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re
import collections

import annotation
import peaks

def overlap(files):
  ids = []
  
  pvalues = collections.defaultdict(float)
  

  for file in files:
    sys.stderr.write("file: " + file + "\n")
    
    matcher = re.match(r'.*Peaks_(.+)\.txt', file)
      
    id = matcher.group(1)
    
    ids.append(id)
    
    sys.stderr.write("ids " + id + "\n")
  
    # now make a list of locations and best p-values  
  
    f = open(file, 'r')
    
    header = f.readline().strip().split("\t")
      
    for line in f:
      line = line.strip()
      
      if len(line) == 0:
        continue
      
      tokens = line.split("\t")
      
      location = tokens[0]
      
      p = float(tokens[annotation.find_heading_index(header, "P-value")])
      
      pvalues[location] = p

    f.close()
  
    
  location_core_map = peaks.overlap_peak_files(files, ids)
  
  # keep the ids sorted
  #ids = sorted(ids)
  
  sys.stdout.write("Genomic Location (hg19)\tRegion / Overlap Width\tBest P-value (ChIPseeqer)\tNumber Of Overlapping Regions\t" + "\t".join(["Sample " + s for s in ids]) + "\n");
  
  for core_location in sorted(location_core_map):
    location = annotation.parse_location(core_location)
    
    overlap = location.end - location.start + 1
    
    c = 0
    
    for id in location_core_map[core_location]:
      #for location in location_core_map[core_location][id]:
      c += len(location_core_map[core_location][id])
        
    p = 0
  
    for id in ids:
      if id in location_core_map[core_location]:
        for location in location_core_map[core_location][id]:
          if pvalues[location] < p:
            p = pvalues[location]
    
    sys.stdout.write("\t".join([core_location, str(overlap), str(p), str(c)]));
  
    # foreach RK id
    for id in ids:
      if id in location_core_map[core_location]:
        p = ";".join(sorted(location_core_map[core_location][id]))
        
        # get rid of internal id
        #peaks = re.sub(r'^.+=', "", peaks)
        
        sys.stdout.write("\t" + p)
      else:
        sys.stdout.write("\tn/a")
    
    sys.stdout.write("\n")



files = sys.argv[1:len(sys.argv)] #sorted(sys.argv[1:len(sys.argv)])

overlap(files)
