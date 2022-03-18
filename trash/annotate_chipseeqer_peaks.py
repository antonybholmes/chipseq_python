# -*- coding: utf-8 -*-
"""
Annotate chipseeqer peaks

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re
import collections

import lib.genomic
import lib.human.genes

MAX_PEAKS_PER_GENE = 3

def annotate(file, prom_ext_5p, prom_ext_3p):
  peak_annotation = lib.human.genes.AnnotatePeak("Peak", prom_ext_5p, prom_ext_3p)
   
  sys.stderr.write("Annotating " + file + "\n")  
  
  f = open(file, "r")
  
  # skip header if not a chipseeqer file
  if "TF_targets" not in file:
    f.readline()
  
  #header
  sys.stdout.write("Genomic Location (hg19)")
  sys.stdout.write("\tP-value (ChIPseeqer)");
  #sys.stdout.write("\tMax Height (reads)");
  sys.stdout.write("\tScore (ChIPseeqer)");
  sys.stdout.write("\tPeak Width");
  sys.stdout.write("\t");
  
  peak_annotation.print_header() 
  
  sys.stdout.write("\tPeaks per gene")
  sys.stdout.write("\n")
  
  locations = []
  annotations = []

  for line in f:
    tokens = line.strip().split("\t")
    
    if lib.genomic.is_location(tokens[0]):
      location = lib.genomic.parse_location(tokens[0])
    else:
      location = lib.genomic.Location(tokens[0], int(tokens[1]), int(tokens[2]))

    #chr = tokens[0]
    
    # Skip chrM
    if "chrM" in location.chr:
      continue
      
    locations.append(location)
    
    #start = int(tokens[1])
    #end = int(tokens[2])
    width = location.end - location.start + 1
    
    if "TF_targets" in file:
      # deal with infs
      if re.match(r'.*[Ii]nf.*', tokens[3]):
        p = -1000
      else:
        p = float(tokens[3])
    
      score = float(tokens[4])
    else:
      p = 0
      score = 0
      
    #max_height = int(tokens[6])
    
    #location = lib.genomic.Location(chr, start, end) #chr + ":" + str(start) + "-" + str(end)
 
    # Write the initial portion
    annotation = [location.to_string(), str(p), str(score), str(width)]
    
    # Add the gene annotation
    annotation.extend(peak_annotation.annotate(location))
    
    annotations.append(annotation)
    
  f.close()
  
  # Flag when a gene has more than 3 peaks
  
  peaks_genes = collections.defaultdict(set)
  genes_peaks = collections.defaultdict(set)
  
  for i in range(0, len(locations)):
    s = locations[i].to_string()
    genes = annotations[i][6].split(";")
    
    for gene in genes:
      if gene != "n/a":
        peaks_genes[s].add(gene)
        genes_peaks[gene].add(s)
  
  for i in range(0, len(locations)):
    s = locations[i].to_string()
    
    c = 0
    
    # see which genes we are on
    if s in peaks_genes:
      genes = peaks_genes[s]
      
      # see how many peaks are on this gene
      for gene in genes:
        if gene in genes_peaks:
          c = max(c, len(genes_peaks[gene]))
    
    annotations[i].append(str(c))
    
  #
  # Output
  #
  for i in range(0, len(locations)):
    sys.stdout.write("\t".join(annotations[i]) + "\n")

  

file = sys.argv[1]
prom_ext_5p = int(sys.argv[2])
prom_ext_3p = int(sys.argv[3])

annotate(file, prom_ext_5p, prom_ext_3p)
