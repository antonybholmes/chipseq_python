# -*- coding: utf-8 -*-
"""
Find the common overlaps between a list of files

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import peaks

file1 = sys.argv[1]
file2 = sys.argv[2]

sys.stderr.write("file1: " + file1 + "\n")
sys.stderr.write("file2: " + file2 + "\n")

lines = peaks.overlap(file1, file2)

print("chr\tstart\tend\tnumber_of_overlapping_peaks\tpeaks");

for line in lines:
  print(line)