# -*- coding: utf-8 -*-
"""
Functions for DNA extraction.

Created on Fri Jan 30 12:37:14 2015

@author: antony
"""


import sys
import urllib
import json

#link = "http://156.145.14.248:8080/dna/api/v2/hg19/chr1/100000/101000/f/u/n"
BASE_LINK = "http://156.145.14.248:8080/dna/api/v2/hg19/"

def get_sequence_with_mask(chr, start, end):
  link = BASE_LINK + chr + "/" + str(start) + "/" + str(end) + "/f/u/n"  

  #sys.stderr.write(link + "\n")  
  
  f = urllib.urlopen(link)
  json_data = f.read()

  data = json.loads(json_data)
  
  sequence = data[0]["sequence"]
  
  return sequence


def score_mask_seq(seq):
  '''
  Score a sequence for percentage of Ns
  '''
  
  p = 0.0
  
  for c in seq:
    #sys.stderr.write(c + "\n")
    if c == 'N':
      p += 1.0
  
  p /= len(seq)
  
  #sys.stderr.write(str(p) + "\n")
  
  return p