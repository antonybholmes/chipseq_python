# -*- coding: utf-8 -*-
"""
Text parsing functions

Created on Sat Jan 31 16:46:42 2015

@author: antony
"""

NA = "n/a"


def get_header(f):
  return f.readline().strip().split("\t")
  
  
def find_index(tokens, text, offset = 0):
  """
  Find the first heading in list that matches some text.
  """
  
  lt = text.lower()
  
  for i in range(offset, len(tokens)):
    if lt in tokens[i].lower():
      return i
      
  return -1


def find_indices(tokens, text, offset = 0):
  """
  Find all the headings in list that matches some text.
  """
  
  indices = []
  
  for i in range(offset, len(tokens)):
    if text.lower() in tokens[i].lower():
      indices.append(i)
      
  return indices


def empty_line(l):
  """
  Produce an empty line
  """
  
  line = NA
  
  for i in range(0, l - 1):
    line += "\t" + NA
      
  return line

