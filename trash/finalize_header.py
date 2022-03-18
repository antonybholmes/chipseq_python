# -*- coding: utf-8 -*-
"""
Rename headers to conform to Katia's specification

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import collections
import re

import annotation

def clean(text):
  ret = text
  
  ret = re.sub(r'RK\d+_', "", ret)
  ret = re.sub(r'_IEO', "", ret)
  
  return ret

def adjust_header(file):
  f = open(file, 'r')

  # skip header
  header = f.readline().strip().split("\t")
  
  new_header = []
  
  for h in header:
    if re.match(r'^RK\d.*', h):
      
      matcher = re.match(r'^((RK\d+)\w+).*', h)
      rk = matcher.group(2)
      id = matcher.group(1)
      
      id = clean(id)
      h = clean(h)
      
      h = re.sub(id, id + " " + rk, h)
    
    new_header.append(h)

  print("\t".join(new_header))  
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    print(line)
    
  f.close()


file = sys.argv[1]

adjust_header(file)