# -*- coding: utf-8 -*-
"""
Functions related to overlaps between peaks and genes

Created on Mon Oct  6 11:52:41 2014

@author: Antony Holmes
"""

import re
from . import text

def get_sample_column_count(header:list[str]) -> int:
  """
  Returns the number of sample columns in a gene file.

  Args:
      header (list[str]): table column header

  Returns:
      int: index of first column containing "Sample"
  """
  ret = 0
  
  c = text.find_index(header, "Sample")
  
  for i in range(c, len(header)):
    if not re.match(r'^Sample.*', header[i]):
      break
    
    ret += 1
    
  return ret
