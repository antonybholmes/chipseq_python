# -*- coding: utf-8 -*-
"""
Lets make some pericentromeric regions

@author: antony
"""

import sys
import re
import collections

import lib.genomic
import lib.sam

PADDING = 6000

bam_file = sys.argv[1]
locations_file = sys.argv[2]
padding = int(sys.argv[3])
window = int(sys.argv[4])
read_length = int(sys.argv[5])

sam = lib.sam.Sam(bam_file)

sys.stderr.write(bam_file + " " + str(sam.get_read_count()) + "\n")

locations, header = lib.genomic.load_locations(locations_file, True, padding, padding)

lib.sam.print_header(padding, padding, window)

bins = lib.sam.get_bins(padding, padding, window)

c = 1

for l in locations:
  #sys.stderr.write(l.to_string() + "\n")
  
  sys.stdout.write(l.to_string())
  
  #starts = sam.get_starts(l)
  
  normalized_counts = sam.get_normalized_start_bins(l, window, read_length)
  
  for i in range(len(bins)):
    sys.stdout.write("\t" + str(normalized_counts[i]))
  
  sys.stdout.write("\n")
  
  if c % 100 == 0:
    sys.stderr.write("Processed " + str(c) + " lines...\n")
  
  c += 1
  
  #break
