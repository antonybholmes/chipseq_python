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

files = sys.argv[1:len(sys.argv)] #sorted(sys.argv[1:len(sys.argv)])
ids = []

pvalues = collections.defaultdict(float)
scores = collections.defaultdict(float)


for file in files:
  sys.stderr.write("file: " + file + "\n")
  
  matcher = re.match(r'.*TF_targets.(.+)_vs.*', file)
    
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

# load the first file, then add the others
nearest_peaks = lib.peaks.NearestPeak()
nearest_peaks.load_from_cols(files[0])
nearest_peaks.load_from_cols(files[1])

nearest_peaks_s1 = lib.peaks.NearestPeak()
nearest_peaks_s1.load_from_cols(files[0])

nearest_peaks_s2 = lib.peaks.NearestPeak()
nearest_peaks_s2.load_from_cols(files[1])

  
location_core_map = lib.peaks.overlapping(files, ids)

# keep the ids sorted
#ids = sorted(ids)

sys.stdout.write("Genomic Location (hg19)\tPeak / Overlap Width\t" + "\t".join([s + " Best P-Value (ChIPseeqer)" for s in ids]) + "\t" + "\t".join([s + " Best Score (ChIPseeqer)" for s in ids]) + "\tNumber Of Overlapping Peaks\t" + "\t".join(["Sample " + s for s in ids]) + "\tClosest Region Distance\tClosest Region Absolute Distance\tClosest Region\n");

for core_location in sorted(location_core_map):
  location = lib.genomic.parse_location(core_location)
  
  overlap = location.end - location.start + 1
  
  c = 0
  
  for id in location_core_map[core_location]:
    for l in location_core_map[core_location][id]:
      c += 1
      
  sys.stdout.write("\t".join([core_location, str(overlap)]))
  
  for s in ids:
    if s in location_core_map[core_location]:
      p = 0
 
      for l in location_core_map[core_location][s]:
        if pvalues[l] < p:
          p = pvalues[l]

      sys.stdout.write("\t" + str(p))
    else:
      sys.stdout.write("\t" + lib.text.NA)
      
  for s in ids:
    if s in location_core_map[core_location]:
      score = 0 
      
      for l in location_core_map[core_location][s]:
        if scores[l] > score:
          score = scores[l]
  
      sys.stdout.write("\t" + str(score))
    else:
      sys.stdout.write("\t" + lib.text.NA)
  
  # number of peaks    
  sys.stdout.write("\t" + str(c));

  # foreach RK id
  for s in ids:
    if s in location_core_map[core_location]:
      peaks = ";".join(sorted(location_core_map[core_location][s]))
      
      # get rid of internal id
      #lib.peaks = re.sub(r'^.+=', "", lib.peaks)
      
      sys.stdout.write("\t" + peaks)
    else:
      sys.stdout.write("\tn/a")
  
  if c == len(ids):
    nearest = nearest_peaks.get_nearest_peak(location)
    d = nearest_peaks.get_nearest_dist(location)
  else:
    if ids[0] in location_core_map[core_location]:
      nearest = nearest_peaks_s2.get_nearest_peak(location)
      d = nearest_peaks_s2.get_nearest_dist(location)
    else:
      nearest = nearest_peaks_s1.get_nearest_peak(location)
      d = nearest_peaks_s1.get_nearest_dist(location)
      
  # Add the closest feature
  
  

  sys.stdout.write("\t" + str(d))
  sys.stdout.write("\t" + str(abs(d)))
  sys.stdout.write("\t" + nearest.to_string())
  
  sys.stdout.write("\n")
