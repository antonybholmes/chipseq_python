# -*- coding: utf-8 -*-
"""
Create a peak file of the shared regions

Created on Fri Sep 12 17:48:31 2014

@author: Antony Holmes
"""

import sys
import collections

import lib.text
import lib.genomic
import lib.headings
import lib.bed
         

def create_bed(file, name, bed_file_1, bed_file_2, bed_file_3):
  sys.stderr.write("three way file " + file  + "\n")
  
  sys.stdout.write("track type=bedGraph name=\"" + name  + "\" description=\"" + name + " peaks\" visibility=full autoScale=on alwaysZero=on color=255,0,0\n")
  #print "track type=bed name=\"", $name, "\" description=\"", $name, "\" itemRgb=255,0,0\n";

  peaks_1 = lib.bed.load_bed_peaks(bed_file_1)
  peaks_2 = lib.bed.load_bed_peaks(bed_file_2)
  peaks_3 = lib.bed.load_bed_peaks(bed_file_3)

  f = open(file, 'r')

  # skip header
  header = f.readline().strip().split("\t")
  
  overlaps_column = lib.text.find_index(header, "Number Of Overlapping Peaks")
  sample_1_column = lib.text.find_index(header, lib.headings.SAMPLE)
  sample_2_column = sample_1_column + 1
  sample_3_column = sample_2_column + 1
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    overlaps = int(tokens[overlaps_column])
    
    if overlaps < 3:
      continue

    location = tokens[0]
    l = lib.genomic.parse_location(location)
    
    c = 0
    average_rpm = 0
    
    
    
    if tokens[sample_1_column] != "n/a":
      location = lib.genomic.parse_location(tokens[sample_1_column]).to_string()
      
      if location in peaks_1:
        average_rpm += peaks_1[location]
        c += 1
    
    if tokens[sample_2_column] != "n/a":
      location = lib.genomic.parse_location(tokens[sample_2_column]).to_string()
      
      if location in peaks_2:
        average_rpm += peaks_2[location]
        c += 1
    
    if tokens[sample_3_column] != "n/a":
      location = lib.genomic.parse_location(tokens[sample_3_column]).to_string()
      
      if location in peaks_3:
        average_rpm += peaks_3[location]
        c += 1
    
    #sys.stderr.write(str(average_rpm) + " " + location + " " + tokens[sample_3_column] + "\n")
    
    average_rpm /= c
    
       
    #sys.stderr.write(file + " " + line + "\n")
    
    # default to writing 100 as an arbitrary measure of the pseudo peak height
    sys.stdout.write("\t".join([l.chr, str(l.start), str(l.end), str(average_rpm)]) + "\n")
    
  f.close()

file = sys.argv[1]
name = sys.argv[2]
bed_file_1 = sys.argv[3]
bed_file_2 = sys.argv[4]
bed_file_3 = sys.argv[5]

create_bed(file, name, bed_file_1, bed_file_2, bed_file_3)
