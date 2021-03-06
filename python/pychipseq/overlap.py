# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 09:27:02 2015

@author: antony
"""

import collections
import sys

import pychipseq.genomic

BIN_SIZE = 1000

class Overlap:
  """
  Uses a gapped search to determine by how much a location overlaps a feature
  """
  def __init__(self, gapped_search):
    self.gapped_search = gapped_search
    
    
  def get_max_overlap(self, location):
    features = self.gapped_search.get_features(location)
      
    max_overlap_width = -1;
    max_overlap = None
    
    for feature in features:
      for l in feature.values:
        overlap = lib_genomic.overlap_locations(location, l)
        
        if overlap is not None:
          if overlap.width > max_overlap_width:
            max_overlap_width = overlap.width
            max_overlap = overlap
            
    return max_overlap


class Overlapping(object):
  """
  Finds the closest peak to another set of peaks.
  """
  
  def __init__(self, file):
    self.bins = collections.defaultdict(lambda: collections.defaultdict(set))
    
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

      location = pychipseq.genomic.parse_location(tokens[0])      
      
      sbin = location.start // BIN_SIZE
      ebin = location.end // BIN_SIZE
      
      for i in range(sbin, ebin + 1):
        self.bins[location.chr][i].add(location)
      
    f.close()
    
  def overlaps(self, location):
    sbin = location.start // BIN_SIZE
    ebin = location.end // BIN_SIZE
    
    overlaps = []
    
    for i in range(sbin, ebin + 1):
      for l in self.bins[location.chr][i]:
        if pychipseq.genomic.is_overlapping(location, l):
          overlaps.append(l)
          
    return overlaps
    
