# -*- coding: utf-8 -*-
"""
Generate a distribution of the peaks closest to a set of regions (peaks)

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import annotation

def closest(peak_file, region_file):
  nearest_peaks = annotation.NearestPeak(peak_file)

  f = open(region_file, "r")

  sys.stdout.write("Region")
  sys.stdout.write("\tNearest Region (bp)");
  sys.stdout.write("\n");
  
  f.readline()
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    location = tokens[0]
    
    overlap_count = int(tokens[3])
    
    nearest = nearest_peaks.get_nearest_peak_from_location(location)
    
    sys.stdout.write("\t".join([location, str(nearest)]) + "\n") 
  
  f.close() 
  

peak_file = sys.argv[2]
region_file = sys.argv[1]

closest(peak_file, region_file)