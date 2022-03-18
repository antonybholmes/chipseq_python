# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 09:43:29 2014

@author: Antony Holmes
"""
import sys
import collections
import re

def parse_peaks(file):
  peaks = collections.defaultdict(lambda: collections.defaultdict(lambda: (int, int)))
  
  f = open(file, "r")
  
  # skip header
  #f.readline()
  
  for line in f:
    tokens = line.split("\t")
    
    chr = tokens[0]
    start = int(tokens[1])
    end = int(tokens[2])
    
    id = ":".join([chr, "-".join([str(start), str(end)])])
    
    #sys.stderr.write(id + "\n")
    
    peaks[chr][id] = (start, end)
    
  f.close()
  
  return peaks
  
def overlap(ref_file, query_file):
  peaks = parse_peaks(ref_file)
  
  f = open(query_file, "r")
  
  # skip header
  #f.readline()
  
  lines = []
  
  for line in f:
    tokens = line.split("\t")
    
    chr = tokens[0]
    start = int(tokens[1])
    end = int(tokens[2])
    
    features = peak_overlap(peaks, chr, start, end)
    
    lines.append("\t".join([chr + ":" + str(start) + "-" + str(end), str(len(features)), ";".join(features)]))
    
  f.close()
  
  return lines

def peak_overlap(peaks, chr, start, end):
  features = []
    
  for id in sorted(peaks[chr]):
    p_start = int(peaks[chr][id][0])
    p_end = int(peaks[chr][id][1])
      
    if start >= p_start and end <= p_end:
      features.append(";".join([id, "within"]))
    elif start < p_start and end > p_end:
      features.append(";".join([id, "over"]))
    elif start < p_start and end > p_start:
      features.append(";".join([id, "upstream"]))
    elif start < p_end and end > p_end:
      features.append(";".join([id, "downstream"]))
      
  return features
  
def overlapping(files, ids):
  """
  Takes a list of file and finds the common (if any) overlapping regions
  """
  
  location_id_map = collections.defaultdict(str)
  location_group_map = collections.defaultdict(str)
  location_chrs = collections.defaultdict(str)
  location_starts = collections.defaultdict(int)
  location_ends = collections.defaultdict(int)
  
  locations = []
  
  for i in range(0, len(files)):
    file = files[i]
    id = ids[i]

    f = open(file, 'r')
    
    c = 0
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = line.split("\t")
    
      chr = tokens[0]
      start = int(tokens[1])
      end = int(tokens[2])
    
      location = id + "=" + chr + ":" + str(start) + "-" + str(end)
    
      locations.append(location)
      
      location_id_map[location] = id
      
      location_chrs[location] = chr
      location_starts[location] = start
      location_ends[location] = end
      
      
      # default to not being in a group
      location_group_map[location] = "none"
      
      c += 1
      
    f.close()
    
    sys.stderr.write(file + " " + str(c) + "\n")

  
  # lets see what overlaps

  # debug for testing to end remove as it truncates list
  #locations = locations[1:(len(locations) / 2)]

  group_id = 0

  total = len(locations)
  
  total *= total
  
  total /= 2
  
  sys.stderr.write("Processing " + str(len(locations)) + " " + str(total) + " locations\n");

  p = 0
  
  
  
  for i in range(0, len(locations)):
    location1 = locations[i]
    
    chr1 = location_chrs[location1]
    start1 = location_starts[location1]
    end1 = location_ends[location1]
    group1 = location_group_map[location1]
    
    #if group1 != "none":
    #  continue
    
    for j in range(0, len(locations)):
      location2 = locations[j]
      
      chr2 = location_chrs[location2]
      start2 = location_starts[location2]
      end2 = location_ends[location2]
      group2 = location_group_map[location2]

      if p % 10000000 == 0:
        sys.stderr.write("p " + str(p) + "\n")

      p += 1      
      
      #if group2 != "none":
      #  continue      
      
      if chr1 != chr2:
        continue
      
      if (start1 <= start2):
        min_start = start1
        min_end = end1
        max_start = start2
        max_end = end2
      else:
        min_start = start2
        min_end = end2
        max_start = start1
        max_end = end1
        
      overlap = -1
      overlap_start = -1
      
      if max_start >= min_start and max_end <= min_end:
        overlap = max_end - max_start + 1
        overlap_start = max_start
      elif min_start < max_start and min_end > max_start:
        overlap = min_end - max_start + 1
        overlap_start = max_start
      else:
        pass
      
      if overlap == -1:
        continue
      
      
      # change the start1 and end1 coordinates to reflect the overlap
      # region so that each subsequent match must be within this region
      # this prevents long peaks that overlap two smaller peaks who
      # themselves do not overlap each other
      start1 = overlap_start
      end1 = overlap_start + overlap - 1
 
      if group1 != "none":
        group = group1
        #sys.stderr.write("reuse group 1 " + group1 + "\n")
      elif group2 != "none":
        group = group2
        #sys.stderr.write("reuse group 2 " + group2 + "\n")
      else:
        group = "group_" + str(group_id)
        group_id += 1

      # update group 1
      group1 = group

      location_group_map[location1] = group
      location_group_map[location2] = group
      
      #sys.stderr.write("overlap " + overlap_location + " " + str(overlap) + " " + location1 + " " + location2 + " " + group + " " + str(group_id) + "\n")
  
  # after iterating over everything, group locations by group

  groups = collections.defaultdict(list)  
  
  for location in location_group_map:
    group = location_group_map[location]
    
    groups[group].append(location)

  # find common overlap location

  overlaps = collections.defaultdict(lambda: collections.defaultdict(list))  

  core_id = 0  
  
  for group in groups:
    if group == "none":
      # write items as is
      for location in groups[group]:
        overlaps[location][location_id_map[location]].append(location)
      
      continue
    
    max_start = -1
    min_end = -1
    min_chr = "n/a"
    
    for location in groups[group]:
      min_chr = location_chrs[location]
      start = location_starts[location]
      end = location_ends[location]
      
      if start > max_start:
        max_start = start
      
      if end < min_end or min_end == -1:
        min_end = end
    
    core_location = "core_" + str(core_id) + "=" + min_chr + ":" + str(max_start) + "-" + str(min_end)

    core_id += 1    
    
    for location in groups[group]:
      overlaps[core_location][location_id_map[location]].append(location)
  
  return overlaps

def overlapping_v2(files, ids):
  """
  Takes a list of file and finds the common (if any) overlapping regions
  """

  bin_size = 10000  
  
  location_id_map = collections.defaultdict(str)
  location_core_map = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(bool)))
  location_chrs = collections.defaultdict(str)
  location_starts = collections.defaultdict(int)
  location_ends = collections.defaultdict(int)

  location_bins = collections.defaultdict(lambda: collections.defaultdict(bool))  
  
  locations = []
  
  for i in range(0, len(files)):
    file = files[i]
    id = ids[i]

    f = open(file, 'r')
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = line.split("\t")
    
      chr = tokens[0]
      start = int(tokens[1])
      end = int(tokens[2])
    
      location = id + "=" + chr + ":" + str(start) + "-" + str(end)
    
      locations.append(location)
      
      location_id_map[location] = id
      
      location_chrs[location] = chr
      location_starts[location] = start
      location_ends[location] = end
      
      bin_start = int(start / bin_size)
      bin_end = int(end / bin_size)
      
      for bin in range(bin_start, bin_end + 1):
        location_bins[bin][location] = True
      
    f.close()
  
  # lets see what overlaps

  # debug for testing to end remove as it truncates list
  #locations = locations[1:(len(locations) / 4)]

  total = len(locations)
  
  total *= total
  
  sys.stderr.write("Processing " + str(len(locations)) + " " + str(total) + " locations\n");

  p = 0
  
  # keep track of all locations that have been allocated at least once
  allocated = collections.defaultdict(bool)
  
  for i in range(0, len(locations)):
    location1 = locations[i]
    
    chr1 = location_chrs[location1]
    start1 = location_starts[location1]
    end1 = location_ends[location1]
    #group1 = location_group_map[location1]
    
    #if group1 != "none":
    #  continue
    
    # find possible overlapping locations

    test_locations = collections.defaultdict(bool)    
    
    bin_start = int(start1 / bin_size)
    bin_end = int(end1 / bin_size)
    
    for bin in range(bin_start, bin_end + 1):
      for location in location_bins[bin]:
        test_locations[location] = True
    
    exhausted = False

    used = collections.defaultdict(bool)    
    
    while not exhausted:
      grouped_locations = [location1]
      start1 = location_starts[location1]
      end1 = location_ends[location1]
    
      for location2 in test_locations:
        if location1 == location2:
          continue
        
        if location2 in used:
          continue

        #sys.stderr.write(location1 + " " + location2 + "\n")      
      
        chr2 = location_chrs[location2]
        start2 = location_starts[location2]
        end2 = location_ends[location2]
        #group2 = location_group_map[location2]

        if p % 10000000 == 0:
          sys.stderr.write("p " + str(p) + "\n")

        p += 1      
      
        #if group2 != "none":
        #  continue      
        
        if chr1 != chr2:
          continue
        
        if (start1 <= start2):
          min_start = start1
          min_end = end1
          max_start = start2
          max_end = end2
        else:
          min_start = start2
          min_end = end2
          max_start = start1
          max_end = end1
          
        overlap = -1
        overlap_start = -1
        
        if max_start >= min_start and max_end <= min_end:
          overlap = max_end - max_start + 1
          overlap_start = max_start
        elif min_start < max_start and min_end > max_start:
          overlap = min_end - max_start + 1
          overlap_start = max_start
        else:
          pass
        
          
        if overlap == -1:
          continue
        
      
        # change the start1 and end1 coordinates to reflect the overlap
        # region so that each subsequent match must be within this region
        # this prevents long peaks that overlap two smaller peaks who
        # themselves do not overlap each other
        start1 = overlap_start
        end1 = overlap_start + overlap - 1
      
        grouped_locations.append(location2)
    
      # now we have a list of all locations that overlap each other
      
      # if we have a group of entries, merge them. otherwise if the
      # location is by itself, only add it if it has not been allocated
      # to another group. This prevents duplicate entries of the whole
      # region by itself plus any overlapping regions
      if len(grouped_locations) > 1 or location1 not in allocated:
        overlap_location = chr1 + ":" + str(start1) + "-" + str(end1)

        for location in grouped_locations:
          id = location_id_map[location]
      
          location_core_map[overlap_location][id][location] = True
        
          used[location] = True
          allocated[location] = True
       
      if len(grouped_locations) == 1:
        # no more to add so quit looping
        exhausted = True
      
      #sys.stderr.write("overlap " + overlap_location + " " + str(overlap) + " " + location1 + " " + location2 + " " + group + " " + str(group_id) + "\n")
  
  # after iterating over everything, group locations by group
  
  return location_core_map