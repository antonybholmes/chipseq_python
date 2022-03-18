# -*- coding: utf-8 -*-
"""
Generate a tss distribution for a chipseeqer peak file

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import collections
import re
import os

def create_bed(file, name, color):
  
  sys.stdout.write("track type=bedGraph name=\"" + name  + "\" description=\"" + name  + " reads\" autoScale=on alwaysZero=on visibility=full color=" + color + "\n")


  bins = collections.defaultdict(int)

  f = open(file, "r")

  
  f.readline()
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    if tokens[3] != "2":
      continue
    
    location = tokens[0]
    
    matcher = re.match(r'(chr.+):(\d+)-(\d+)', location)
    
    chr = matcher.group(1)
    start = matcher.group(2)
    end = matcher.group(3)

    #default
    height = "20"    
    
    sys.stdout.write("\t".join([chr, start, end, height]) + "\n")

  f.close()


file = sys.argv[1]
name = sys.argv[2]
color = sys.argv[3]

create_bed(file, name, color)