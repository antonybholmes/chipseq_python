# -*- coding: utf-8 -*-
"""
Functions related to overlaps between peaks and genes

Created on Mon Oct  6 11:52:41 2014

@author: Antony Holmes
"""

import re

import lib.genes
import lib.text
import lib.headings

def get_sample_column_count(header):
  """
  Returns the number of sample columns in a gene file
  """
  
  ret = 0
  
  c = lib.text.find_index(header, "Sample")
  
  for i in range(c, len(header)):
    if not re.match(r'^Sample.*', header[i]):
      break
    
    ret += 1
    
  return ret
