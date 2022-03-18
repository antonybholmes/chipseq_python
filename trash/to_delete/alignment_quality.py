# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 10:27:24 2014

@author: Antony Holmes
"""

import re
import sys


def qualities(sample, file, unique_reads):

  # lets get the sample name
  sys.stderr.write("Reading file " + file + "\n")
  
  # sample consists of id + two letter digit combo from facility
  #matcher = re.match(r'.*\/(.+[A-Z]{2}\d+)\/.*', file)
  #id = matcher.group(1)
  
  #matcher = re.match(r'.*([A-Z]{2}\d+).*', name)

  #id = matcher.group(1)

  
  
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
  
  duplicate_reads = total_mapped_reads - unique_reads
  
  p_dup = float(duplicate_reads) / float(total_reads) * 100
  
  p = float(total_mapped_reads) / float(total_reads) * 100
  #no_pcr_dup_p = float(non_duplicate_reads) / float(total_mapped_reads)
  unique_reads_p = float(unique_reads) / float(total_reads) * 100
  
  sys.stdout.write("Name\tReads\tMapped Reads\tDuplicate Reads\t% Duplicate Reads\tUnique Reads\t% Unique Reads\n")
  #sys.stdout.write("\t".join([sample, "{:,}".format(total_reads), "{:,}".format(total_mapped_reads), "{:,}".format(non_duplicate_count), "{:.2f}".format(p)]) + "\n")  
  sys.stdout.write("\t".join([sample, "{:,}".format(total_reads), "{:,}".format(total_mapped_reads), "{:,}".format(duplicate_reads), "{:.2f}".format(p_dup), "{:,}".format(unique_reads), "{:.2f}".format(unique_reads_p)]) + "\n")  

sample = sys.argv[1]
file = sys.argv[2]
unique_reads = int(sys.argv[3])

qualities(sample, file, unique_reads)  
