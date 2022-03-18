# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 20:32:59 2015

@author: antony
"""

import sys

import lib_genomic


def centromeres(file):
  centromeres = lib_genomic.SearchGenomicBedFeatures("/ifs/scratch/cancer/Lab_RDF/abh2138/references/ucsc/ucsc_centromeres_hg19.bed")
  pericentromeres = lib_genomic.SearchGenomicBedFeatures("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rdf_pericentromeres_hg19.bed")
  
  cen_overlaps = lib_genomic.GenomicFeaturesOverlap(centromeres)
  p_cen_overlaps = lib_genomic.GenomicFeaturesOverlap(pericentromeres)
  
  f = open(file, 'r')
  
  sys.stdout.write(f.readline().strip() + "\tCentromere Region\n")
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    location = lib_genomic.parse_location(line)
    
    classification = "n/a"
    
    #
    # Are we in a centromere
    #    

    in_centromere = False    
    
    overlap = cen_overlaps.get_max_overlap(location)
    
    p = 0
    
    if overlap is not None:
      p = overlap.length / location.length
      
      classification = "centromere"
      in_centromere = True
    
    
    #
    # Are we in a pericentromere
    #
    
    if not in_centromere:
      overlap = p_cen_overlaps.get_max_overlap(location)
      
      if overlap is not None:
        p = overlap.length / location.length
        classification = "pericentromere"


      
    sys.stdout.write(line + "\t" + classification + "\n")
      
  f.close()
  
  
file = sys.argv[1]

centromeres(file)