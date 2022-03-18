# -*- coding: utf-8 -*-
"""
Generate a distribution of the peaks closest to a set of regions (peaks)

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import collections
import numpy

import lib.peaks
import lib.text
import lib.headings
import lib.genomic

def closest(peak_file, region_file):
  nearest_peaks = lib.peaks.NearestPeak(peak_file)
  
  nearest_map = collections.defaultdict(list)
  location_map = collections.defaultdict(list)
  nearest_location_map = collections.defaultdict(list)

  f = open(region_file, "r")

  sys.stdout.write("Nearest Region Distance (bp)\tRegion\tNearest Region\n");
  
  header = f.readline().strip().split("\t")
  
  location_column = lib.text.find_index(header, lib.headings.LOCATION)
  overlap_column = lib.text.find_index(header, lib.headings.OVERLAPPING)
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    location = lib.genomic.parse_location(tokens[location_column])
    
    #sys.stderr.write("o " + str(overlap_column) + "\n")
    
    if overlap_column != -1:
      overlap_count = int(tokens[overlap_column])
    
      # Only keep the intersecting regions
      if overlap_count < 2:
        continue
    
    nearest_location = nearest_peaks.get_nearest_peak(location)
    
    if nearest_location is None:
      continue
    
    nearest = nearest_peaks.get_nearest_dist(location) #get_nearest_peak(location)
    
    abs_nearest = abs(nearest)
    
    location_map[abs_nearest].append(location)
    nearest_map[abs_nearest].append(nearest)
    nearest_location_map[abs_nearest].append(nearest_location)
  
  f.close()   
    
  for abs_nearest in sorted(nearest_map):
    sorted_indices = numpy.argsort(nearest_map[abs_nearest])
    
    for i in sorted_indices:
      location = location_map[abs_nearest][i]
      nearest = nearest_map[abs_nearest][i]
      nearest_location = nearest_location_map[abs_nearest][i]
      
      #sys.stderr.write(str(nearest) + " " + location.to_string() + "\n")
      
      sys.stdout.write("\t".join([str(nearest), location.to_string(), nearest_location.to_string()]) + "\n") 
  
  
  
region_file = sys.argv[1]
peak_file = sys.argv[2]


closest(peak_file, region_file)
