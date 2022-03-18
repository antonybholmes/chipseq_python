# -*- coding: utf-8 -*-
"""
Functions related to samples

Created on Sat Jan 31 16:46:42 2015

@author: antony
"""

import sys
import re


def get_sample_id(text):
  """
  Get the unique RK id of a sample.
  """
  
  # Id consists of two letters and a three digit number
  return re.match(r'.*?_([A-Z]{2}\d{3}).*', text).group(1)
  


