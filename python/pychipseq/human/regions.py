# -*- coding: utf-8 -*-
"""
Functions related to overlaps between peaks and genes

Created on Mon Oct  6 11:52:41 2014

@author: Antony Holmes
"""

import re

import pychipseq.human.peaks
import pychipseq.text
import pychipseq.headings
import pychipseq.regions
import pychipseq.genes
  
class GeneOrientatedRegionIntersections(pychipseq.human.peaks.GeneOrientatedPeaks):
  """
  Annotate genes that are part of intersections only.
  """
  
  def __init__(self):
    super().__init__("Region")
    
  
  def load_annotations(self, file):
    """
    Override to cope with checking for true intersections
    """
    
    f = open(file, 'r')

    # skip header
    header = f.readline().strip().split("\t")
  
    location_column = pychipseq.text.find_index(header, pychipseq.headings.LOCATION)
    entrez_column = pychipseq.text.find_index(header, pychipseq.headings.ENTREZ_ID)
    refseq_column = pychipseq.text.find_index(header, pychipseq.headings.REFSEQ_ID)
    symbol_column = pychipseq.text.find_index(header, pychipseq.headings.GENE_SYMBOL)
    type_column = pychipseq.text.find_index(header, "Relative To Gene")
    #p_column = pychipseq.text.find_index(header, pychipseq.headings.P_VALUE)
    #score_column = pychipseq.text.find_index(header, pychipseq.headings.SCORE)
    tss_column = pychipseq.text.find_index(header, pychipseq.headings.TSS_DISTANCE)
    sample_column = pychipseq.text.find_index(header, pychipseq.headings.SAMPLE)
    centromere_column = pychipseq.text.find_index(header, pychipseq.headings.CENTROMERE)
    
    sample_columns = pychipseq.regions.get_sample_column_count(header)

    mir_column = pychipseq.text.find_index(header, pychipseq.headings.MIR_SYMBOL)
    mir_type_column = pychipseq.text.find_index(header, "Relative To miR")
    mss_column = pychipseq.text.find_index(header, "mIR Start Distance")
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = ls.split("\t")

      type = tokens[type_column]
      p = pychipseq.genes.find_best_p_value(header, tokens) #float(tokens[p_column])
      score = pychipseq.genes.find_best_score(header, tokens) #float(tokens[score_column])
      location = tokens[location_column]
      centromere = tokens[centromere_column]
      
      intersection = True

      #sys.stderr.write("sample column " + str(sample_column) + " " + str(sample_columns) + "\n")      
      
      for i in range(0, sample_columns):
        if tokens[sample_column + i] == pychipseq.text.NA:
          intersection = False
          break
        
      if not intersection:
        continue
      
      entrezes = tokens[entrez_column].split(";")
      symbols = tokens[symbol_column].split(";")
      refseqs = tokens[refseq_column].split(";")
      tsses = tokens[tss_column].split(";")
      
      for i in range(0, len(refseqs)):
        refseq = refseqs[i]
        entrez = entrezes[i]
        symbol = symbols[i]
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
          self.collapsed_centromeres[refseq].append(centromere)
        
      mir = tokens[mir_column]
    
      if mir != "n/a":
        self.mirs.add(mir)
        self.collapsed_types[mir].append(tokens[mir_type_column])
        self.collapsed_p[mir].append(p)
        self.collapsed_scores[mir].append(score)
        self.collapsed_locations[mir].append(location)
        self.collapsed_tss[mir].append(tokens[mss_column])
        self.collapsed_centromeres[mir].append(centromere)
    
    f.close()
    

class ClosestGeneOrientatedRegionIntersections(pychipseq.human.peaks.GeneOrientatedPeaks):
  def __init__(self):
    super().__init__("Region")
    
  
  def load_annotations(self, file):
    """
    Override to cope with checking for true intersections
    """
    
    f = open(file, 'r')

    # skip header
    header = f.readline().strip().split("\t")
  
    location_column = pychipseq.text.find_index(header, pychipseq.headings.LOCATION)
    p_column = pychipseq.text.find_index(header, pychipseq.headings.P_VALUE)
    score_column = pychipseq.text.find_index(header, pychipseq.headings.SCORE)
    centromere_column = pychipseq.text.find_index(header, pychipseq.headings.CENTROMERE)

    entrez_column = pychipseq.text.find_index(header, pychipseq.headings.ENTREZ_ID)
    refseq_column = pychipseq.text.find_index(header, pychipseq.headings.REFSEQ_ID)
    symbol_column = pychipseq.text.find_index(header, pychipseq.headings.GENE_SYMBOL)
    type_column = pychipseq.text.find_index(header, pychipseq.headings.RELATIVE)
    tss_column = pychipseq.text.find_index(header, pychipseq.headings.TSS_DISTANCE)
    
    closest_entrez_column = pychipseq.text.find_index(header, pychipseq.headings.CLOSEST_ENTREZ_ID)
    closest_refseq_column = pychipseq.text.find_index(header, pychipseq.headings.CLOSEST_REFSEQ_ID)
    closest_symbol_column = pychipseq.text.find_index(header, pychipseq.headings.CLOSEST_GENE_SYMBOL)
    closest_type_column = pychipseq.text.find_index(header, pychipseq.headings.CLOSEST_RELATIVE)
    closest_tss_column = pychipseq.text.find_index(header, pychipseq.headings.CLOSEST_TSS_DISTANCE)
    
    
    mir_column = pychipseq.text.find_index(header, pychipseq.headings.MIR_SYMBOL)
    mir_type_column = pychipseq.text.find_index(header, "Relative To miR")
    mss_column = pychipseq.text.find_index(header, "mIR Start Distance")
    
    
    sample_column = pychipseq.text.find_index(header, pychipseq.headings.SAMPLE)
    
    sample_columns = pychipseq.regions.get_sample_column_count(header)

    mir_column = pychipseq.text.find_index(header, pychipseq.headings.MIR_SYMBOL)
    mir_type_column = pychipseq.text.find_index(header, "Relative To miR")
    mss_column = pychipseq.text.find_index(header, "mIR Start Distance")
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = ls.split("\t")

      
      intersection = True
      
      for i in range(0, sample_columns):
        if tokens[sample_column + i] == pychipseq.text.NA:
          intersection = False
          break
        
      if not intersection:
        continue
        
      type = tokens[type_column]
      p = float(tokens[p_column])
      score = float(tokens[score_column])
      location = tokens[location_column]
      centromere = tokens[centromere_column]
      
      if tokens[entrez_column] != pychipseq.text.NA:
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
        refseq = refseqs[i]
        entrez = entrezes[i]
        symbol = symbols[i]
        tss = tsses[i]
        
        if refseq != pychipseq.text.NA:
          self.refseqs.add(refseq)
          self.collapsed_entrezes[refseq] = entrez
          self.collapsed_symbols[refseq] = symbol
          self.collapsed_p[refseq].append(p)
          self.collapsed_scores[refseq].append(score)
          self.collapsed_locations[refseq].append(location)
          self.collapsed_tss[refseq].append(tss)
          self.collapsed_types[refseq].append(type)
          self.collapsed_centromeres[refseq].append(centromere)
        
      mir = tokens[mir_column]
    
      if mir != pychipseq.text.NA:
        self.mirs.add(mir)
        self.collapsed_types[mir].append(tokens[mir_type_column])
        self.collapsed_p[mir].append(p)
        self.collapsed_scores[mir].append(score)
        self.collapsed_locations[mir].append(location)
        self.collapsed_tss[mir].append(tokens[mss_column])
        self.collapsed_centromeres[mir].append(centromere)
    
    f.close()
    

class GeneOrientatedRegionCores(pychipseq.human.peaks.GeneOrientatedPeaks):
  """
  For annotating peaks that lie within any of the three way intersections
  """
  
  def __init__(self):
    super().__init__("Region")
    
  
  def load_annotations(self, file):
    """
    Override to cope with checking for true intersections
    """
    
    f = open(file, 'r')

    # skip header
    header = f.readline().strip().split("\t")
  
    location_column = pychipseq.text.find_index(header, pychipseq.headings.LOCATION)
    entrez_column = pychipseq.text.find_index(header, pychipseq.headings.ENTREZ_ID)
    refseq_column = pychipseq.text.find_index(header, pychipseq.headings.REFSEQ_ID)
    symbol_column = pychipseq.text.find_index(header, pychipseq.headings.GENE_SYMBOL)
    type_column = pychipseq.text.find_index(header, "Relative To Gene")
    p_column = pychipseq.text.find_index(header, pychipseq.headings.P_VALUE)
    score_column = pychipseq.text.find_index(header, pychipseq.headings.SCORE)
    tss_column = pychipseq.text.find_index(header, pychipseq.headings.TSS_DISTANCE)
    centromere_column = pychipseq.text.find_index(header, pychipseq.headings.CENTROMERE)
    sample_column = pychipseq.text.find_index(header, pychipseq.headings.SAMPLE)
    
    mir_column = pychipseq.text.find_index(header, pychipseq.headings.MIR_SYMBOL)
    mir_type_column = pychipseq.text.find_index(header, "Relative To miR")
    mss_column = pychipseq.text.find_index(header, "mIR Start Distance")
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = ls.split("\t")

      type = tokens[type_column]
      p = float(tokens[p_column])
      score = float(tokens[score_column])
      location = tokens[location_column]
      centromere = tokens[centromere_column]
      
      intersection = False
      
      # One of the intersections must be true
      if (tokens[sample_column] != "n/a" and tokens[sample_column + 1] != "n/a") or \
        (tokens[sample_column] != "n/a" and tokens[sample_column + 2] != "n/a") or \
        (tokens[sample_column + 1] != "n/a" and tokens[sample_column + 2] != "n/a"):
        intersection = True
        
      if not intersection:
        continue
      
      
      entrezes = tokens[entrez_column].split(";")
      symbols = tokens[symbol_column].split(";")
      refseqs = tokens[refseq_column].split(";")
      tsses = tokens[tss_column].split(";")
      
      for i in range(0, len(refseqs)):
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
          self.collapsed_centromeres[refseq].append(centromere)
        
      mir = tokens[mir_column]
    
      if mir != "n/a":
        self.mirs.add(mir)
        self.collapsed_types[mir].append(tokens[mir_type_column])
        self.collapsed_p[mir].append(p)
        self.collapsed_scores[mir].append(score)
        self.collapsed_locations[mir].append(location)
        self.collapsed_tss[mir].append(tokens[mss_column])
        self.collapsed_centromeres[mir].append(centromere)
    
    f.close()
