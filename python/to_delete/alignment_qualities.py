# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 10:27:24 2014

@author: Antony Holmes
"""

import re
import sys


def qualities(file):

  # lets get the sample name
  #sys.stderr.write("Reading file " + file + "\n")
  
  # sample consists of id + two letter digit combo from facility
  matcher = re.match(r'.*\/(.+[A-Z]{2}\d+)\/.*', file)

  name = matcher.group(1)
  
  matcher = re.match(r'.*([A-Z]{2}\d+).*', name)

  id = matcher.group(1)

  
  
  # open the file and find the stats
  
  f = open(file, 'r')
  
  for line in f:
		line = line.strip()
		
		if "reads" in line:
			break

  matcher = re.match(r'^(\d+).*', line)
  
  total_reads = int(matcher.group(1))
  
  # skip
  f.next()
  f.next()
  
  line = f.next().strip()
  matcher = re.match(r'^(\d+).*', line)
  one_reads = int(matcher.group(1))
  
  line = f.next().strip()
  matcher = re.match(r'^(\d+).*', line)
  multiple_reads = int(matcher.group(1))
  
  total_mapped_reads = one_reads + multiple_reads

  #sys.stderr.write(str(reads) + " " + str(total_mapped_reads) + "\n"); 
  
  p = float(total_mapped_reads) / float(total_reads)
  
  #sys.stdout.write("id\treads\ttotal_mapped_reads\tp\n")
  sys.stdout.write("\t".join([id, "{:,}".format(total_reads), "{:,}".format(total_mapped_reads), "{:.2f}".format(p)]) + "\n")  
  #sys.stdout.write("\t".join([id, str(reads), str(total_mapped_reads), "{:.2f}".format(p)]) + "\n")  

file = sys.argv[1]

qualities(file)  
