# -*- coding: utf-8 -*-
"""
Functions for processing long non-coding RNA

Created on Tue Aug 26 16:36:13 2014

@author: Antony Holmes
"""

import collections
import sys
import re

class LncRNA(object):
  def __init__(self):
    self.file = "/ifs/scratch/cancer/Lab_RDF/abh2138/references/gencode/gencode.v19.long_noncoding_RNAs.gtf"   
    self.transcript_starts = collections.defaultdict(lambda: collections.defaultdict(set)) 
    self.starts = collections.defaultdict(list)

    starts = collections.defaultdict(lambda: collections.defaultdict(int))
    ends = collections.defaultdict(lambda: collections.defaultdict(int))
    all_starts = collections.defaultdict(set) 

    sys.stderr.write("Loading transcripts from " + self.file + "...\n")
  
    f = open(self.file, 'r')

    for line in f:
      line = line.strip()

      if len(line) == 0:
        continue
            
      if re.match(r'^#.+', line):
        continue
      
      tokens = line.split("\t")
      
      chr = tokens[0]
      type = tokens[2]
      start = int(tokens[3])
      end = int(tokens[4])
      strand = tokens[6]
      annotation = tokens[8]
      
      if type != "transcript":
        continue
      
      matcher = re.match(r'.+transcript_id "([^"]+).+', annotation)

      transcript_id = matcher.group(1)
      
      #sys.stderr.write(transcript_id + "\n")

      starts[chr][transcript_id] = start   
      ends[chr][transcript_id] = end   
      
      all_starts[chr].add(start)
      all_starts[chr].add(end)
      
    f.close()
    
    for chr in all_starts:
      self.starts[chr] = sorted(all_starts[chr])
      
    ranks = collections.defaultdict(lambda: collections.defaultdict(int))
    
    for chr in self.starts:
      rank = 0
      
      for start in self.starts[chr]:
        ranks[chr][start] = rank
        
        rank += 1
        
    # now fill the gaps by adding transcripts to the positions between their
    # start and end
        
    for transcript_id in starts[chr]:
      start_rank = ranks[chr][starts[chr][transcript_id]]
      end_rank = ranks[chr][ends[chr][transcript_id]]
           
      for rank in range(start_rank, end_rank + 1):
        self.transcript_starts[chr][self.starts[chr][rank]].add(transcript_id)
    
    sys.stderr.write("Finished loading transcripts.\n")

  def get_lnrna_from_location(self, location):
    
    matcher = re.match(r'(chr.+):(\d+)-(\d+)', location)
    
    chr = matcher.group(1)
    start = int(matcher.group(2))
    end = int(matcher.group(3))
    
    return self.get_lnrna_from_position(chr, start, end)
  
  def get_lnrna_from_position(self, \
    chr, \
    start, \
    end):
        
    start_index = self.get_index(chr, start)
    end_index = self.get_index(chr, end)
    
    lncrna = set()

    
    # part overlap
    if start_index == -1 and end_index != -1:
      start_index = 0

    # part overlap
    if start_index != -1 and end_index == -1:
      end_index = len(self.starts[chr]) - 1
    
    if start_index != -1 and end_index != -1:
      for i in range(start_index, end_index + 1):
        s = self.starts[chr][i]
        
        for transcript_id in self.transcript_starts[chr][s]:
          lncrna.add(transcript_id)
    
    if len(lncrna) == 0:
      lncrna.add("n/a")
    
    return lncrna
    
  def get_index(self, \
    chr, \
    start):
         
    # use a binary style search to find closest peak
    
    ps = 0
    pe = len(self.starts[chr]) - 1

    # No point searching if outside range
    if start < self.starts[chr][ps] or start > self.starts[chr][pe]:
      return -1

    while pe - ps > 1:
      pm = int((pe + ps) / 2)
      test_start = self.starts[chr][pm]
      
      # perfect match      
      if test_start == start:
        return pm
      elif test_start > start:
        pe = pm
      else:
        ps = pm
          
    return ps