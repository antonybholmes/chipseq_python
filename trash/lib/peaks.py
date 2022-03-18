# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 14:13:05 2014

@author: Antony Holmes
"""

import sys
import collections
import re

import lib.text
import lib.genomic
import lib.headings

BIN_SIZE = 10000

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


def overlapping_peaks(files, ids):
  """
  Takes a list of file and finds the common (if any) overlapping regions
  """

  location_id_map = collections.defaultdict(str)
  location_chrs = collections.defaultdict(str)
  location_starts = collections.defaultdict(int)
  location_ends = collections.defaultdict(int)
  location_bins = collections.defaultdict(set)  
  locations = []
  
  for i in range(0, len(files)):
    file = files[i]
    id = ids[i]
    
    sys.stderr.write("hmm " + file + "\n")

    f = open(file, 'r')
    
    # Skip header
    if "Peaks" in file:
      f.readline()
    
    for line in f:
      tokens = line.strip().split("\t")
      
      if lib.genomic.is_location(tokens[0]):
        location = lib.genomic.parse_location(tokens[0])
      else:
        if lib.genomic.is_chr(tokens[0]):
          location = lib.genomic.Location(tokens[0], int(tokens[1]), int(tokens[2]))
        else:
          sys.stderr.write("Invalid line: " + line + "\n")
          continue

      lid = id + "=" + location.chr + ":" + str(location.start) + "-" + str(location.end)
      
      #sys.stderr.write("lid " + lid + "\n")
    
      locations.append(lid)
      
      location_id_map[lid] = id
      
      location_chrs[lid] = location.chr
      location_starts[lid] = location.start
      location_ends[lid] = location.end
      
      bin_start = int(location.start / BIN_SIZE)
      bin_end = int(location.end / BIN_SIZE)
      
      for bin in range(bin_start, bin_end + 1):
        location_bins[bin].add(lid)
      
    f.close()
    
  return overlapping(ids, location_id_map, location_chrs, location_starts, location_ends, location_bins, locations)
    

def overlapping_peak_tables(files, ids):
  """
  Takes a list of file and finds the common (if any) overlapping regions
  """

  location_id_map = collections.defaultdict(str)
  location_chrs = collections.defaultdict(str)
  location_starts = collections.defaultdict(int)
  location_ends = collections.defaultdict(int)
  location_bins = collections.defaultdict(set)  
  locations = []
  
  for i in range(0, len(files)):
    file = files[i]
    id = ids[i]

    f = open(file, 'r')
    
    f.readline()
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = line.split("\t")
    
      location = lib.genomic.parse_location(tokens[0])
    
      lid = id + "=" + location.chr + ":" + str(location.start) + "-" + str(location.end)
    
      locations.append(lid)
      
      #sys.stderr.write(id + "\n")
      location_id_map[lid] = id
      
      location_chrs[lid] = location.chr
      location_starts[lid] = location.start
      location_ends[lid] = location.end
      
      bin_start = int(location.start / BIN_SIZE)
      bin_end = int(location.end / BIN_SIZE)
      
      for bin in range(bin_start, bin_end + 1):
        location_bins[bin].add(lid)
      
    f.close()
    
  return overlapping(ids, location_id_map, location_chrs, location_starts, location_ends, location_bins, locations)
    
    
def overlapping(ids, location_id_map, location_chrs, location_starts, location_ends, location_bins, locations):
  """
  Calculates the maximum overlaps between a set of sample locations
  """
  
  # lets see what overlaps
  
  location_core_map = collections.defaultdict(lambda: collections.defaultdict(str))

  # debug for testing to end remove as it truncates list
  #locations = locations[1:(len(locations) / 4)]

  total = len(locations)
  
  total *= total
  
  sys.stderr.write("Processing " + str(len(locations)) + " locations...\n");

  p = 0
  
  # keep track of all locations that have been allocated at least once
  allocated = set()
  
  for i in range(0, len(locations)):
    # of the form id=chrN:start-end
    location1 = locations[i]
    
    chr1 = location_chrs[location1]
    start1 = location_starts[location1]
    end1 = location_ends[location1]
    #group1 = location_group_map[location1]
    
    #if group1 != "none":
    #  continue
    
    # find possible overlapping locations

    test_locations = set()
    
    bin_start = int(start1 / BIN_SIZE)
    bin_end = int(end1 / BIN_SIZE)
    
    for bin in range(bin_start, bin_end + 1):
      for location in location_bins[bin]:
        test_locations.add(location)

    used = set()
    
    # Form the largest group of overlapping peaks
    exhausted = False
    
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
        
        # Given the two locations we are testing, sort them so
        # the starts and ends are in order.
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
        
        if max_start == min_start and max_end > min_end:
          # Peak 1 and 2 start at the same point but peak 2 is wider
          # so the overlap region is peak 1
          overlap = min_end - min_start + 1
          overlap_start = min_start
        elif max_start >= min_start and max_end <= min_end:
          # Peak 1 is wider than peak 2 and contains it
          # so the overlap region is peak 2
          overlap = max_end - max_start + 1
          overlap_start = max_start
        elif min_start < max_start and min_end > max_start:
          # Peak one starts before peak 2 but ends within peak 2
          # so the overlap is the start of peak 2 to the end of
          # peak 1
          overlap = min_end - max_start + 1
          overlap_start = max_start
        else:
          pass
        
        # We have not found an overlap yet so continue
        if overlap == -1:
          continue
        
      
        # change the start1 and end1 coordinates to reflect the overlap
        # region so that each subsequent match must be within this region
        # this prevents long peaks that overlap two smaller peaks who
        # themselves do not overlap each other
        start1 = overlap_start
        end1 = start1 + overlap - 1
      
        grouped_locations.append(location2)
    
      # now we have a list of all locations that overlap each other
      
      # if we have a group of entries, merge them, otherwise if the
      # location is by itself, only add it if it has not been allocated
      # to another group. This prevents duplicate entries of the whole
      # region by itself plus any overlapping regions
      if len(grouped_locations) > 1 or location1 not in allocated:
        overlap_location = chr1 + ":" + str(start1) + "-" + str(end1)

        for location in grouped_locations:
          # id is a sample id
          id = location_id_map[location]
      
          #sys.stderr.write("overlap " + overlap_location + " " + id + " " + location + "\n")
      
          location_core_map[overlap_location][id] = location #.add(location)
        
          used.add(location)
          allocated.add(location)
       
      if len(grouped_locations) == 1:
        # no more to add so quit looping
        exhausted = True
      
  # after iterating over everything, group locations by group
  
  return location_core_map

  
def duplicate_peaks(type, file):
  f = open(file, "r")
  
  header = f.readline().strip().split("\t")
  
  entrez_column = lib.text.find_index(header, lib.headings.ENTREZ_ID)
  refseq_column = lib.text.find_index(header, lib.headings.REFSEQ_ID)
  symbol_column = lib.text.find_index(header, lib.headings.GENE_SYMBOL)
  overlap_type_column = lib.text.find_index(header, type + " Relative To Gene")
  tss_column = lib.text.find_index(header, type + " TSS Distance")
  
  mir_column = lib.text.find_index(header, "miR Symbol")
  mir_type_column = lib.text.find_index(header, type + " Relative To miR")
  mss_column = lib.text.find_index(header, type + " miR Start Distance")
     
  sys.stdout.write("\t".join(header) + "\n")
  
  for line in f:
    ls = line.strip()
      
    if len(ls) == 0:
      continue
      
    tokens = ls.split("\t");
    
    entrezes = tokens[entrez_column].split(";")
    refseqs = tokens[refseq_column].split(";")
    overlap_types = tokens[overlap_type_column].split(";")
    symbols = tokens[symbol_column].split(";")
    tsses = tokens[tss_column].split(";")

    split = False
    
    for i in range(0, len(refseqs)):
      refseq = refseqs[i]
      
      
      if refseq == lib.text.NA:
        continue
      
      entrez = entrezes[i]
      overlap_type = overlap_types[i]
      symbol = symbols[i]
      tss = tsses[i]
      
      # clone
      new_tokens = tokens[:]
      
      new_tokens[overlap_type_column] = overlap_type
      new_tokens[entrez_column] = entrez
      new_tokens[refseq_column] = refseq
      new_tokens[symbol_column] = symbol
      new_tokens[tss_column] = tss
      
      new_tokens[mir_column] = lib.text.NA
      new_tokens[mir_type_column] = lib.text.NA
      
      
      sys.stdout.write("\t".join(new_tokens) + "\n")
      
      split = True
      
    
    # Only process the mirs if they exist
    
    if mir_column != -1:
      mirs = tokens[mir_column].split(";") 
      mir_types = tokens[mir_type_column].split(";") 
      mir_mss = tokens[mss_column].split(";")
      
      for i in range(0, len(mirs)):
        mir = mirs[i]
        
        if mir == lib.text.NA:
          continue
        
        mir_type = mir_types[i]
        mss = mir_mss[i]
        
        new_tokens = tokens[:]
        
        new_tokens[overlap_type_column] = lib.text.NA
        new_tokens[entrez_column] = lib.text.NA
        new_tokens[refseq_column] = lib.text.NA
        new_tokens[symbol_column] = lib.text.NA
        new_tokens[tss_column] = lib.text.NA
        
        new_tokens[mir_column] = mir
        new_tokens[mir_type_column] = mir_type
        new_tokens[mss_column] = mss
        
        sys.stdout.write("\t".join(new_tokens) + "\n")
        
        split = True
    
    
    if split == False:
      # no lines were split so print the row as it was originally
      sys.stdout.write("\t".join(tokens) + "\n")
        
  f.close()


def filter_peaks(file, max_tss_5p_dist, max_tss_3p_dist):
  """
  Exclude peaks that are too far from a tss
  """
  
  f = open(file, "r")
  
  header = f.readline().strip().split("\t")
  
  tss_column = lib.text.find_index(header, lib.headings.TSS_DISTANCE)
      
  sys.stdout.write("\t".join(header) + "\n")
  
  for line in f:
    line = line.strip()
      
    if len(line) == 0:
      continue
      
    tokens = line.split("\t");
    
    tss = tokens[tss_column]

    if tss == lib.text.NA:
      continue
    
    d = int(tss)
    
    if d < max_tss_5p_dist or d > max_tss_3p_dist:
      continue
      
    sys.stdout.write("\t".join(tokens) + "\n")
        
  f.close()
  

class NearestPeak(object):
  """
  Finds the closest peak to another set of peaks.
  """
  
  def __init__(self, file = None):
    self.peak_starts = collections.defaultdict(lambda: collections.defaultdict(lib.genomic.Location))
    self.starts = collections.defaultdict(list)
    
    if file is not None:
      self.load(file)
      
    
  
  def load(self, file):
    sys.stderr.write("Loading peaks from " + file + "...\n")
  
    f = open(file, 'r')

    # skip header
    f.readline()
  
    for line in f:
      line = line.strip()
      
      if len(line) == 0:
        continue
      
      tokens = line.split("\t")

      location = lib.genomic.parse_location(tokens[0])      
      
      mid_point = lib.genomic.mid_point(location)
      
      self.peak_starts[location.chr][mid_point] = location
      
    f.close()
    
    for chr in self.peak_starts:
      self.starts[chr] = sorted(self.peak_starts[chr])
    
    sys.stderr.write("Finished loading peaks.\n")
    
  def load_from_cols(self, file):
    sys.stderr.write("Loading peaks from " + file + "...\n")
  
    f = open(file, 'r')

    # skip header
    f.readline()
  
    for line in f:
      line = line.strip()
      
      if len(line) == 0:
        continue
      
      tokens = line.split("\t")

      location = lib.genomic.parse_location_cols(tokens)      
      
      mid_point = lib.genomic.mid_point(location)
      
      self.peak_starts[location.chr][mid_point] = location
      
    f.close()
    
    for chr in self.peak_starts:
      self.starts[chr] = sorted(self.peak_starts[chr])
    
    sys.stderr.write("Finished loading peaks.\n")
    

  def get_nearest_peak(self, location):
    chr = location.chr
    
    if len(self.starts[chr]) == 0:
      return None

    mid_point = lib.genomic.mid_point(location)
    
    # use a binary style search to find closest peak
    
    ps = 0
    pe = len(self.starts[chr]) - 1
 
    while pe - ps > 1:
      pm = int((pe + ps) / 2)
      test_start = self.starts[chr][pm]
      
      #sys.stderr.write(str(test_start) + " " + str(mid_point) + " " + str(ps) + " " + str(pe) + " " + str(pm) + "\n")  
      
      # perfect match      
      if mid_point == test_start:
        return self.peak_starts[chr][mid_point]
      elif test_start > mid_point:      
        pe = pm
      else:
        ps = pm
    
    # point lies between two midpoints, so pick the closest

    sd = mid_point - self.starts[chr][ps]
    ed = mid_point - self.starts[chr][pe]
    
    abss = abs(sd)
    abse = abs(ed)
     
    if (abss <= abse):
      return self.peak_starts[chr][self.starts[chr][ps]]
    else:
      return self.peak_starts[chr][self.starts[chr][pe]]
  
  
  def get_nearest_abs_dist(self, location):
    return abs(self.get_nearest_dist(location))
  
    
  def get_nearest_dist(self, location):
    """
    Find the closest distance to a feature
    """
    
    chr = location.chr
    
    if len(self.starts[chr]) == 0:
      return sys.maxint
        
    mid_point = lib.genomic.mid_point(location)
    
    # use a binary style search to find closest peak
    
    ps = 0
    pe = len(self.starts[chr]) - 1
 
    while pe - ps > 1:
      pm = int((pe + ps) / 2)
      test_start = self.starts[chr][pm]
      
      #sys.stderr.write(str(test_start) + " " + str(mid_point) + " " + str(ps) + " " + str(pe) + " " + str(pm) + "\n")  
      
      # perfect match      
      if mid_point == test_start:
        return 0
      elif test_start > mid_point:      
        pe = pm
      else:
        ps = pm
    
    # point lies between two midpoints, so pick the closest
    
    # sys.stderr.write("hmm " + chr + " " + str(ps) + "\n")
    
    sd = mid_point - self.starts[chr][ps]
    ed = mid_point - self.starts[chr][pe]
    
    abss = abs(sd)
    abse = abs(ed)
     
    if (abss <= abse):
      return sd
    else:
      return ed
