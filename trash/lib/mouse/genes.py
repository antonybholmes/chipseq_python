# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 15:12:34 2014

@author: Antony Holmes
"""

import sys
import collections
import re

import lib.annotation
import lib.expression
import tss

import lib.genomic
import lib.headings
import lib.text
import lib.sample
import lib.genes


class GeneOrientatedPeaks(object):
  """
  Core annotation for gene oriented peaks
  """
  
  def __init__(self, type):
    self.type = type

    self.collapsed_entrezes = collections.defaultdict(str)
    self.collapsed_symbols = collections.defaultdict(str)
    self.collapsed_p = collections.defaultdict(list)
    self.collapsed_scores = collections.defaultdict(list)  
    self.collapsed_locations = collections.defaultdict(list)
    self.collapsed_tss = collections.defaultdict(list)
    self.collapsed_types = collections.defaultdict(list)
 
    self.refseqs = set()
    
  
  def load_annotations(self, file):
    f = open(file, 'r')

    # skip header
    header = f.readline().strip().split("\t")
  
    location_column = lib.text.find_index(header, lib.headings.LOCATION)
    entrez_column = lib.text.find_index(header, lib.headings.ENTREZ_ID)
    refseq_column = lib.text.find_index(header, lib.headings.REFSEQ_ID)
    symbol_column = lib.text.find_index(header, lib.headings.GENE_SYMBOL)
    type_column = lib.text.find_index(header, "Relative To Gene")
    p_column = lib.text.find_index(header, lib.headings.P_VALUE)
    score_column = lib.text.find_index(header, lib.headings.SCORE)
    tss_column = lib.text.find_index(header, lib.headings.TSS_DISTANCE)
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = ls.split("\t")

      type = tokens[type_column]
      p = float(tokens[p_column])
      score = float(tokens[score_column])
      location = tokens[location_column]
      
      entrezes = tokens[entrez_column].split(";")
      symbols = tokens[symbol_column].split(";")
      refseqs = tokens[refseq_column].split(";")
      tsses = tokens[tss_column].split(";")

      for i in range(0, len(refseqs)):
        # The core id identifies the unique genes on a peak, rather than
        #        
        entrez = entrezes[i]
        symbol = symbols[i]
        refseq = refseqs[i]
        tss = tsses[i]

        if refseq != "n/a":
          self.refseqs.add(refseq)
          self.collapsed_entrezes[refseq] = entrez
          self.collapsed_symbols[refseq] = symbol
          self.collapsed_p[refseq].append(p)
          self.collapsed_scores[refseq].append(score)
          self.collapsed_locations[refseq].append(location)
          self.collapsed_tss[refseq].append(tss)
          self.collapsed_types[refseq].append(type)
    
    f.close()
    
    
  def print_header(self):
    sys.stdout.write(lib.headings.REFSEQ_ID)
    sys.stdout.write("\t" + lib.headings.ENTREZ_ID)
    sys.stdout.write("\t" + lib.headings.GENE_SYMBOL)

    sys.stdout.write("\t" + self.type + " Relative To Gene")
    sys.stdout.write("\t" + self.type + " TSS Closest Distance")
    sys.stdout.write("\t" + self.type + " " + lib.headings.TSS_DISTANCE)
    sys.stdout.write("\tBest P-value (ChIPseeqer)")
    sys.stdout.write("\tBest Score (ChIPseeqer)")
    sys.stdout.write("\t" + self.type + " Count")
    sys.stdout.write("\t" + self.type + " Genomic Locations (mm10)")
    sys.stdout.write("\n");

  
  def get_ids(self):
    return sorted(self.refseqs)

  def gene_orient_peak(self, id):
    entrez = self.collapsed_entrezes[id]
    
    sys.stdout.write(id)
    sys.stdout.write("\t" + self.collapsed_entrezes[id])
    sys.stdout.write("\t" + self.collapsed_symbols[id])

    sys.stdout.write("\t" + ";".join(self.collapsed_types[id]))
    
    # if there are some nearest tss, print the closest
    sys.stdout.write("\t" + lib.genomic.get_closest_tss(self.collapsed_tss[id]))
          
    sys.stdout.write("\t" + ";".join(self.collapsed_tss[id]))
      
    # pick the smallest p
    p = sorted(self.collapsed_p[id])
    sys.stdout.write("\t" + str(p[0]))
    
    # pick the largest score
    scores = sorted(self.collapsed_scores[id], reverse=True)
    sys.stdout.write("\t" + str(scores[0]))
      
    # peak count
    sys.stdout.write("\t" + str(len(self.collapsed_locations[id])))
      
    sys.stdout.write("\t" + ";".join(self.collapsed_locations[id]))
    
    sys.stdout.write("\n");


class ClosestGeneOrientatedPeaks(GeneOrientatedPeaks):
  """
  Produce a gene oriented version of peaks using the peaks within
  a gene/promoter first and then the intergenic ones last. This is
  a hybrid of the two annotations with preference given to 
  peaks overlapping genes
  """
  
  def load_annotations(self, file):
    f = open(file, 'r')

    # skip header
    header = f.readline().strip().split("\t")

    # independent columns
    location_column = lib.text.find_index(header, lib.headings.LOCATION)
    p_column = lib.text.find_index(header, lib.headings.P_VALUE)
    score_column = lib.text.find_index(header, lib.headings.SCORE)
    centromere_column = lib.text.find_index(header, lib.headings.CENTROMERE)

    entrez_column = lib.text.find_index(header, lib.headings.ENTREZ_ID)
    refseq_column = lib.text.find_index(header, lib.headings.REFSEQ_ID)
    symbol_column = lib.text.find_index(header, lib.headings.GENE_SYMBOL)
    type_column = lib.text.find_index(header, lib.headings.RELATIVE)
    tss_column = lib.text.find_index(header, lib.headings.TSS_DISTANCE)
    
    closest_entrez_column = lib.text.find_index(header, lib.headings.CLOSEST_ENTREZ_ID)
    closest_refseq_column = lib.text.find_index(header, lib.headings.CLOSEST_REFSEQ_ID)
    closest_symbol_column = lib.text.find_index(header, lib.headings.CLOSEST_GENE_SYMBOL)
    closest_type_column = lib.text.find_index(header, lib.headings.CLOSEST_RELATIVE)
    closest_tss_column = lib.text.find_index(header, lib.headings.CLOSEST_TSS_DISTANCE)
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = ls.split("\t")

      location = tokens[location_column]
      p = float(tokens[p_column])
      score = float(tokens[score_column])
      
      centromere = tokens[centromere_column]
      
      if tokens[entrez_column] != lib.text.NA:
        # Preference is given to peaks within a gene
        entrezes = tokens[entrez_column].split(";")
        symbols = tokens[symbol_column].split(";")
        refseqs = tokens[refseq_column].split(";")
        tsses = tokens[tss_column].split(";")
        type = tokens[type_column]
      else:
        entrezes = tokens[closest_entrez_column].split(";")
        symbols = tokens[closest_symbol_column].split(";")
        refseqs = tokens[closest_refseq_column].split(";")
        tsses = tokens[closest_tss_column].split(";")
        type = tokens[closest_type_column]
      
      
        
      for i in range(0, len(refseqs)):
        # The core id identifies the unique genes on a peak, rather than
        #        
        entrez = entrezes[i]
        symbol = symbols[i]
        refseq = refseqs[i]
        tss = tsses[i]

        if refseq != lib.text.NA:
          self.refseqs.add(refseq)
          self.collapsed_entrezes[refseq] = entrez
          self.collapsed_symbols[refseq] = symbol
          self.collapsed_p[refseq].append(p)
          self.collapsed_scores[refseq].append(score)
          self.collapsed_locations[refseq].append(location)
          self.collapsed_tss[refseq].append(tss)
          self.collapsed_types[refseq].append(type)
    
    f.close()


class RefSeqGenes(lib.genes.RefSeqGenes):
  def __init__(self):
    super(RefSeqGenes, self).__init__(lib.mouse.annotation.REFSEQ_FILE)


class AnnotatePeak(lib.genes.AnnotatePeak):
  """
  Core annotation for annotating peaks/regions
  """
  
  def __init__(self, type, prom_ext_5p=5000, prom_ext_3p=4000, bin_size=10000):
    lib.genes.AnnotatePeak.__init__(self, \
      type, \
      lib.mouse.tss.RefSeqAnnotation(prom_ext_5p, prom_ext_3p, bin_size), \
      RefSeqGenes(), \
      lib.mouse.tss.RefSeqTss(prom_ext_5p, prom_ext_3p), \
      lib.mouse.tss.RefSeqEnd(prom_ext_5p, prom_ext_3p), \
      prom_ext_5p, \
      prom_ext_3p, \
      bin_size)
