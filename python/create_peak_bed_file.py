# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 17:48:31 2014

@author: Antony Holmes
"""

import sys
import re
import os.path

import pychipseq.bed

def create_bed(file, name, genome):
  
  # get the number of mapped reads
  
  print("xpath " + file + "\n", file=sys.stderr)
  
  found = True
  
  path = "../../alignment_quality_" + genome + ".txt"
  
  if os.path.exists(path):
    f = open(path, "r")
    sys.stderr.write("Found " + path + "\n")
  else:
    path = "../../alignment_quality.txt"
    
    if os.path.exists(path):
      f = open(path, "r")
      sys.stderr.write("Found " + path + "\n")
    else:
      found = False
  
  if found:
    f.readline()
  
    tokens = f.readline().strip().split("\t")
  
    # trim commas so number can be parsed
    mapped_reads = float(re.sub(r',', '', tokens[1]))
    
    f.close()
  else:
    mapped_reads = -1

    
  
  print(pychipseq.bed.create_bedgraph_header(name)) #"track type=bedGraph name=\"" + name  + "\" description=\"" + name + " peaks\" visibility=full autoScale=on alwaysZero=on color=255,0,0\n")
  #print "track type=bed name=\"", $name, "\" description=\"", $name, "\" itemRgb=255,0,0\n";
  
  
  print("run " + file + "\n", file=sys.stderr)
  
  f = open(file, 'r')
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    chr = tokens[0]
    
    start = int(tokens[1])
    end = int(tokens[2])

    max_peak_height = float(tokens[6])
    
    #col4 = '{}:{}-{}'.format(chr, start, end)
    col4 = max_peak_height
    
    #if mapped_reads != -1:
    #  rpm = max_peak_height / mapped_reads * 1000000
    #else:
    #  rpm = max_peak_height
    
    #sys.stderr.write(str(max_peak_height) + " " + str(mapped_reads) + " " + str(rpm) + "\n")
    #sys.stdout.write("\t".join([tokens[0], tokens[1], tokens[2], "{:.2f}".format(rpm)]) + "\n")
    #pychipseq.bed.write_bedgraph_line(chr, start, end, rpm)
    
    print('{}\t{}\t{}\t{}'.format(chr, start, end, col4))
    
  f.close()
  

file = sys.argv[1]
name = sys.argv[2]
genome = sys.argv[3]

create_bed(file, name, genome)
