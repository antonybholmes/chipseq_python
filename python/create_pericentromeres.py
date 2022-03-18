# -*- coding: utf-8 -*-
"""
Lets make some pericentromeric regions

@author: antony
"""

import sys
import re
import collections

import lib_dna

sizes = collections.defaultdict(int)

file = "/ifs/scratch/cancer/Lab_RDF/abh2138/references/ucsc/hg19_chromosome_sizes.txt"

f = open(file, 'r')

f.readline()

for line in f:
  line = line.strip()

  if len(line) == 0:
    continue
  
  tokens = line.split("\t")
  
  chr = tokens[0]  
  
  size = int(tokens[1])
  
  sizes[chr] = size

f.close()


file = "/ifs/scratch/cancer/Lab_RDF/abh2138/references/ucsc/ucsc_cytobands_hg19.txt"

cen_starts = collections.defaultdict(list)
cen_ends = collections.defaultdict(list)
p_cen_starts = collections.defaultdict(list)
p_cen_ends = collections.defaultdict(list)

jump = 1000
search_range = 100

f = open(file, 'r')

f.readline()

for line in f:
  line = line.strip()

  if len(line) == 0:
    continue
  
  tokens = line.split("\t")
  
  if tokens[4] != "acen":
    continue

  chr = tokens[0]
  
  sys.stdout.write(chr + "\n")
  
  cen_start_1 = int(tokens[1])
  cen_end_1 = int(tokens[2])
  
  cen_starts[chr].append(cen_start_1)
  cen_ends[chr].append(cen_end_1)
  
  tokens = f.next().strip().split("\t")
  
  cen_start_2 = int(tokens[1])
  cen_end_2 = int(tokens[2])
  
  cen_starts[chr].append(cen_start_2)
  cen_ends[chr].append(cen_end_2)


  #  
  # first the pericentromere start
  #
  
  start = cen_start_1 - jump - 1;
  end = 0
  
  while start >= 0:
    # get x bases
    sequence = lib_dna.get_sequence_with_mask(chr, start, start + jump - 1)
    
    found = False
    
    for i in range(jump - search_range - 1, -1, -1):
      if 'N' not in sequence[i:(i + search_range + 1)]:
        end = start + i + search_range
        found = True
        break
    
    if found:
      break
    
    #sys.stdout.write(sequence)
    
    start -= jump
    
  p_cen_starts[chr].append(end + 1)
  p_cen_ends[chr].append(cen_start_1 - 1)
    
  sys.stderr.write("p1 end " + str(end + 1) + " " + str(cen_start_1 - 1) + "\n")
  
  #sequence = lib_dna.get_sequence_with_mask(chr, start, end)


  #  
  # Now the end
  #

  start = cen_end_2 + 1
  end = sizes[chr] - 1
  
  while end < sizes[chr]:
    sequence = lib_dna.get_sequence_with_mask(chr, start, start + jump - 1)
    
    found = False
    
    for i in range(1, jump - search_range):
      if 'N' not in sequence[i:(i + search_range - 1)]:
        end = start + i
        found = True
        break
    
    if found:
      break
    
    start += jump
    #end += jump
  
  
  sys.stderr.write("p2 end " + str(cen_end_2 + 1) + " " + str(end - 1) + "\n")
  

  p_cen_starts[chr].append(cen_end_2 + 1)
  p_cen_ends[chr].append(end - 1)  
  
  #sequence = lib_dna.get_sequence_with_mask(chr, start, end)
  
  #sys.stdout.write(sequence + "\n")
  
  
  
  #lib_dna.get_sequence_with_mask(chr, start, end)
  
  #break

f.close()

f = open("ucsc_centromeres_hg19.bed", 'w')

for chr in sorted(cen_starts):
  for i in range(0, len(cen_starts[chr])):  
    f.write("\t".join([chr, str(cen_starts[chr][i]), str(cen_ends[chr][i]), "centromere_" + str(i + 1)]) + "\n")

f.close()

f = open("rdf_pericentromeres_hg19.bed", 'w')

for chr in sorted(cen_starts):
  for i in range(0, len(p_cen_starts[chr])):
    f.write("\t".join([chr, str(p_cen_starts[chr][i]), str(p_cen_ends[chr][i]), "pericentromere_" + str(i + 1)]) + "\n")

f.close()