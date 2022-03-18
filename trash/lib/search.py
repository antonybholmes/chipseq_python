# -*- coding: utf-8 -*-
"""
Binary style search to quickly indentify the closest objects to
a location.

Created on Wed Feb  4 09:58:09 2015

@author: antony
"""

import collections

class Features(object):
  def __init__(self, start):
    self.start = start
    self.values = []
      
      
class GappedSearch(object):
  def __init__(self):
    self.features = collections.defaultdict(lambda : collections.defaultdict(Features))
    self.indices = collections.defaultdict(int)
    self.locked = False
      
  def add_feature(self, location, feature):
    self.locked = False
      
    if location.chr not in self.features or location.start not in self.features[location.chr]:
      self.features[location.chr][location.start] = Features(location.start)

    if location.chr not in self.features or location.end not in self.features[location.chr]:
      self.features[location.chr][location.end] = Features(location.end)
      
    self.features[location.chr][location.start].values.append(feature)
    self.features[location.chr][location.end].values.append(feature)
  
  
  def lock(self):
    """
    Sort the indices for searching
    """
    
    if self.locked:
      return
    
    self.locked = True
    
    # We need the locations sorted in order for a binary search
    for chr in self.features:
      self.indices[chr] = sorted(self.features[chr])
  
  
  def get_closest_features(self, location):
    features = self.get_features(location);
		
    min_d = sys.maxint
    ret = None
    
    for feature in features:
      d = abs(lib_genomic.mid(location) - feature.start)
      
      if d < min_d:
        ret = feature
        min_d = d
    
    return ret
    
    
  def get_features(self, location):
    """
    Return a list of features closest to a location which then
    be tested for overlap etc
    """
    
    ret = []
    
    if location.chr not in self.features:
      return ret
      
    self.lock()
    
    indices = self.indices[location.chr]
    
    s = self.get_start_index(indices, location.start)
    e = self.get_end_index(indices, location.end)
   
    chr_features = self.features[location.chr]

    for i in range(s, e + 1):
      ret.append(chr_features[indices[i]])
	
    return ret
    
  
  def get_start_index(self, indices, start):
    if len(indices) < 2:
      return 0
 
    s = 0

    if start <= indices[s]:
      return s

    e = len(indices) - 1
    
    if start >= indices[e]:
      return e;
		
    while e - s > 1:
      im = (e + s) / 2
			
      pm = indices[im]
      
      if pm > start:
        e = im
      elif pm < start:
        s = im
      else:
        return im
		
    # If we made it this far, we have narrowed down the search to
    # two position so return the end	
    return s
	
 
  def get_end_index(self, indices, end):
    """
    Return the end
    """
    if len(indices) < 2:
      return 0
		
    s = 0
    
    if end <= indices[s]:
      return 0

    e = len(indices) - 1;
		
    if end >= indices[e]:
      return e
 
    while e - s > 1:
      im = (e + s) / 2
 
      pm = indices[im]
      
      if pm > end:
        e = im
      elif pm < end:
        s = im
      else:
        return im;
       
    return e
    
    
class BlockSearch(object):
  def __init__(self, block_size = 10000):
    self.block_size = block_size
    self.features = collections.defaultdict(lambda : collections.defaultdict(Features))
      
  def add_feature(self, location, feature):
    s = location.start // self.block_size
    e = location.end // self.block_size
     
    for b in range(e, e + 1):
      if location.chr not in self.features or b not in self.features[location.chr]:
        self.features[location.chr][b] = Features(location.start)
         
      self.features[location.chr][b].values.append(feature)
        
  
  def get_closest_features(self, location):
    features = get_features(location);
		
    min_d = sys.maxint
    ret = None
    
    for feature in features:
      d = abs(lib_genomic.mid(location) - feature.start)
      
      if d < min_d:
        ret = feature
        min_d = d
    
    return ret
    
    
  def get_features(self, location):
    """
    Return a list of features closest to a location which then
    be tested for overlap etc
    """
    
    ret = []
    
    if location.chr not in self.features:
      return ret
    
    s = location.start // self.block_size
    e = location.end // self.block_size
   
    chr_features = self.features[location.chr]

    for b in range(s, e + 1):
      if b in chr_features:
        ret.append(chr_features[b])
	
    return ret
