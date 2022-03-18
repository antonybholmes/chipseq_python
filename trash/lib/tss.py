# -*- coding: utf-8 -*-
"""
Classes to deal with finding the TSS closest to a point

Created on Mon Aug 25 17:25:56 2014

@author: Antony Holmes
"""

import collections
import sys
import re

import lib.genomic
import lib.sample
import lib.annotation


class TssGene(object):
  def __init__(self, id, types, d):
    self.id = id
    self.types = types
    self.d = d


class TssGenes(object):
  def __init__(self, d, genes):
    self.d = d
    self.genes = genes


class RefSeqTssClassification(object):
  """
  Classify genes by type
  """
  
  def __init__(self, file, prom_ext_5p, prom_ext_3p):
    self.variants = collections.defaultdict(set)    
    self.gene_strands = collections.defaultdict()    
    self.gene_starts = collections.defaultdict(int)
    self.gene_ends = collections.defaultdict(int)
    self.promoter_starts = collections.defaultdict(int)
    self.promoter_ends = collections.defaultdict(int)
    self.exon_counts = collections.defaultdict(int)
    self.exon_starts = collections.defaultdict(list)
    self.exon_ends = collections.defaultdict(list)
    
    # how to label promoters
    self.promoter_label = "promoter" #"promoter_-" + str(prom_ext_5p / 1000) + "/+" + str(prom_ext_3p / 1000)
    
    
    sys.stderr.write("Loading classifications from " + file + "...\n")
  
    f = open(file, 'r')

    # skip header
    f.readline()
  
    for line in f:
      line = line.strip()
      
      if len(line) == 0:
        continue
      
      tokens = line.split("\t")
      refseq = tokens[0]      
      entrez = tokens[1]
      symbol = tokens[2]
      chr = tokens[3]
      strand = tokens[4]
      # ucsc convention
      start = int(tokens[5]) + 1
      end = int(tokens[6])
      exon_count = int(tokens[7])

      if re.match(r'.*n/a.*', refseq):
        continue
      
      if re.match(r'.*n/a.*', entrez):
        continue
      
      if re.match(r'.*n/a.*', symbol):
        continue
      
      if re.match(r'.*MIR.*', symbol):
        continue
      
      # UCSC convention
      ex_starts = [int(p) + 1 for p in tokens[8].split(",")]
      ex_ends = [int(p) for p in tokens[9].split(",")]
      
      if strand == "+":
        promoter_start = start - prom_ext_5p
        promoter_end = start + prom_ext_3p
      else:
        promoter_start = end - prom_ext_3p
        promoter_end = end + prom_ext_5p

      variant_id = lib.sample.create_variant(refseq, chr, start, end)
      
      self.variants[refseq].add(variant_id)
      self.gene_strands[variant_id] = strand
      self.gene_starts[variant_id] = start
      self.gene_ends[variant_id] = end
      self.promoter_starts[variant_id] = promoter_start
      self.promoter_ends[variant_id] = promoter_end
      self.exon_counts[variant_id] = exon_count
        
      for p in ex_starts:
        self.exon_starts[variant_id].append(p)
          
      for p in ex_ends:
        self.exon_ends[variant_id].append(p)
    
    f.close()
    
    sys.stderr.write("Finished loading genes.\n")
    
  
  def get_classification(self, \
    variant_id, \
    mid_point):
      
    refseq = lib.sample.parse_id_from_variant(variant_id)
      
    types = set()
    
     
    # Test all variants with the same refseq and start as our variant
    # Although we picked a closest variant, we should test all variants
    # with the same start/end to determine the type of the gene

    
    min_d = -1;
    min_abs_d = -1;
    
    for var_id in self.variants[refseq]:
      
      strand = self.gene_strands[var_id]
      
      # skip variants that do not begin or end where the reference
      # variant does.
      if strand == "+":
        if self.gene_starts[var_id] != self.gene_starts[variant_id]:
          continue
      else:
        if self.gene_ends[var_id] != self.gene_ends[variant_id]:
          continue
    
      gene_start = self.gene_starts[var_id]
      gene_end = self.gene_ends[var_id]
            
      promoter_start = self.promoter_starts[var_id]
      promoter_end = self.promoter_ends[var_id]
    
      if strand == "+":
        d = mid_point - gene_start
      else:
        d = gene_end - mid_point
        
      abs_d = abs(d)
      
      if abs_d < min_abs_d or min_abs_d == -1:
        min_abs_d = abs_d
        min_d = d
      
      # Peak must fall within gene boundaries
      within_bounds = False
      
      if strand == "+":
        if mid_point >= promoter_start and mid_point <= gene_end:
          within_bounds = True
      else:
        # negative strand
        if mid_point >= gene_start and mid_point <= promoter_end:
          within_bounds = True
      
      if within_bounds:
        # Peak must be within gene boundary to be called promoter, exon or intron
        if mid_point >= promoter_start and mid_point <= promoter_end:
          types.add(self.promoter_label)
        
        in_exon = False
        
        for i in range(0, len(self.exon_starts[var_id])):  
          if mid_point >= self.exon_starts[var_id][i] and mid_point <= self.exon_ends[var_id][i]:
            types.add("exonic")
            in_exon = True
            break
        
        # if you're not exonic and you are within the gene, you must
        # be intronic
        if not in_exon and mid_point >= gene_start and mid_point <= gene_end:
          types.add("intronic")
      
    
    if "exonic" in types and "intronic" in types:
      # We favor classifying as intronic where possible
      types.remove("exonic")
    
    if len(types) == 0:
      #No types implies integenic
      types.add("intergenic")
    
    gene = TssGene(variant_id, types, min_d)
    
    # Return in sorted so that promoter comes first
    #return reversed(sorted(types))
    
    return gene


class RefSeqTss(object):
  """
  Determines the closest gene(s) to a given position.
  """
  
  def __init__(self, file, prom_ext_5p, prom_ext_3p):
    self.gene_starts = collections.defaultdict(lambda: collections.defaultdict(set)) 
    self.starts = collections.defaultdict(list)
    self.refseq = RefSeqTssClassification(file, prom_ext_5p, prom_ext_3p)
  
    sys.stderr.write("Loading RefSeq genes from " + file + "...\n")
  
    f = open(file, 'r')

    # skip header
    f.readline()
  
    # To account for multiple versions of a gene, allocate each entrez
    # id a unique index
  
    for line in f:
      line = line.strip()
      
      if len(line) == 0:
        continue
      
      tokens = line.split("\t")
      
      refseq = tokens[0]
      entrez = tokens[1]
      symbol = tokens[2]
      chr = tokens[3]
      strand = tokens[4]
      # ucsc convention
      start = int(tokens[5]) + 1
      end = int(tokens[6])
      
      if re.match(r'.*n/a.*', refseq):
        continue
        
      if re.match(r'.*n/a.*', entrez):
        continue
      
      if re.match(r'.*n/a.*', symbol):
        continue
      
      if re.match(r'.*MIR.*', symbol):
        continue
      
      # If the gene is on the negative strand, the promoter is around the
      # end coordinate on the forward strand
  
      variant_id = lib.sample.create_variant(refseq, chr, start, end)
      
      if strand == "+":
        self.gene_starts[chr][start].add(variant_id)
      else:
        self.gene_starts[chr][end].add(variant_id)
      
    f.close()
    
    # Create sorted lists of start locations for each gene
    for chr in self.gene_starts:
      self.starts[chr] = sorted(self.gene_starts[chr])
    
    sys.stderr.write("Finished loading genes.\n")


  def get_closest_gene(self, location):
    """
    Returns the genes that are closest to a given location.
    
    @param chr      The chromosome.
    @param start    The start position.
    @param end      The end position
    @return         A list of genes with the distance of their TSS from
                    the peak.
    """
    
    mid_point = lib.genomic.mid_point(location)
 
    # use a binary style search to find closest peak
    
    ps = 0
    pe = len(self.starts[location.chr]) - 1

    while pe - ps > 1:
      pm = int((pe + ps) / 2)
      
      test_start = self.starts[location.chr][pm]
      
      # perfect match      
      if mid_point == test_start:
        return self.get_annotation(location.chr, mid_point, test_start)
      elif test_start > mid_point:
        pe = pm
      else:
        ps = pm
    
    # The peak lies between two starts, pick the closest
    abss = abs(mid_point - self.starts[location.chr][ps])
    abse = abs(mid_point - self.starts[location.chr][pe])
    
    if abss <= abse:
      gene = self.get_annotation(location.chr, mid_point, self.starts[location.chr][ps])           
    else:
      gene = self.get_annotation(location.chr, mid_point, self.starts[location.chr][pe])
      
    return gene
  
    
  def get_annotation(self, chr, mid_point, start):
    # Pick the first one we find
    variant_id = sorted(self.gene_starts[chr][start])[0]
    
    gene = self.refseq.get_classification(variant_id, mid_point)

    return gene
    

class RefSeqEnd(RefSeqTss):
  """
  Find the closest gene to a given location using the genes end
  location.
  """
  
  def __init__(self, file, prom_ext_5p, prom_ext_3p):
    self.gene_starts = collections.defaultdict(lambda: collections.defaultdict(set)) 
    self.starts = collections.defaultdict(list)
    self.refseq = RefSeqTssClassification(file, prom_ext_5p, prom_ext_3p)
  
    sys.stderr.write("Loading RefSeq genes from " + file + "...\n")
  
    f = open(file, 'r')

    # skip header
    f.readline()
  
    # To account for multiple versions of a gene, allocate each entrez
    # id a unique index
  
    for line in f:
      line = line.strip()
      
      if len(line) == 0:
        continue
      
      tokens = line.split("\t")
      
      refseq = tokens[0]
      entrez = tokens[1]
      symbol = tokens[2]
      chr = tokens[3]
      strand = tokens[4]
      # ucsc convention
      start = int(tokens[5]) + 1
      end = int(tokens[6])
      
      if re.match(r'.*n/a.*', entrez):
        continue
      
      if re.match(r'.*n/a.*', symbol):
        continue
      
      if re.match(r'.*MIR.*', symbol):
        continue
      
      # If the gene is on the negative strand, the promoter is around the
      # end coordinate on the forward strand
  
      variant_id = lib.sample.create_variant(refseq, chr, start, end)
      
      #sys.stderr.write("var id " + variant_id + "\n")
      
      if strand == "+":
        self.gene_starts[chr][end].add(variant_id)
      else:
        self.gene_starts[chr][start].add(variant_id)
      
    f.close()
    
    # Create sorted lists of start locations for each gene
    for chr in self.gene_starts:
      self.starts[chr] = sorted(self.gene_starts[chr])
    
    sys.stderr.write("Finished loading genes.\n")
  

class RefSeqAnnotation(object):
  """
  Annotate peaks for all possible refseq genes they might overlap.
  """
  
  def __init__(self, file, prom_ext_5p, prom_ext_3p, bin_size):
    """
    Create a new lib.annotation object.
    
    @param prom_ext_5p   How far upstream should be considered 
                         a promoter.
    @param prom_ext_3p   How far downstream should be considered 
                         a promoter.
    @param bin_size      The size of the bins to group genes into.
    """
    
    self.bin_size = bin_size
    
    self.entrezes = collections.defaultdict(str)
    self.gene_strands = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(str)))
    self.gene_starts = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
    self.gene_ends = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
    self.promoter_starts = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
    self.promoter_ends = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
    self.exon_counts= collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
    self.exon_starts = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
    self.exon_ends = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
    self.refseq = RefSeqTssClassification(file, prom_ext_5p, prom_ext_3p)
    
    sys.stderr.write("Loading genes from " + file + "...\n")
  
    f = open(file, 'r')

    # skip header
    f.readline()
  
    for line in f:
      line = line.strip()
      
      if len(line) == 0:
        continue
      
      tokens = line.split("\t")
      refseq = tokens[0]
      entrez = tokens[1]
      symbol = tokens[2]
      chr = tokens[3]
      strand = tokens[4]
      # ucsc convention
      start = int(tokens[5]) + 1
      end = int(tokens[6])
      exon_count = int(tokens[7])

      if re.match(r'.*n/a.*', entrez):
        continue
      
      if re.match(r'.*n/a.*', symbol):
        continue
      
      if re.match(r'.*MIR.*', symbol):
        continue
      
      self.entrezes[refseq] = entrez
      
      ex_starts = [int(p) + 1 for p in tokens[8].split(",")]
      ex_ends = [int(p) for p in tokens[9].split(",")]
      
      max_offset = max(prom_ext_5p, prom_ext_3p)
      
      if strand == "+":
        promoter_start = start - prom_ext_5p
        promoter_end = start + prom_ext_3p
      else:
        promoter_start = end - prom_ext_3p
        promoter_end = end + prom_ext_5p
      
      start_bin = int((start - max_offset) / bin_size)
      end_bin = int((end + max_offset) / bin_size)
      
      variant_id = lib.sample.create_variant(refseq, chr, start, end)
      
      for bin in range(start_bin, end_bin + 1):
        # apparently refseq genes from the ucsc always report
        # coordinates on the forward strand regardless of orientation
        self.gene_strands[chr][bin][variant_id] = strand
        self.gene_starts[chr][bin][variant_id] = start
        self.gene_ends[chr][bin][variant_id] = end
        self.promoter_starts[chr][bin][variant_id] = promoter_start
        self.promoter_ends[chr][bin][variant_id] = promoter_end
        self.exon_counts[chr][bin][variant_id] = exon_count
        
        for p in ex_starts:
          self.exon_starts[chr][bin][variant_id].append(p)
          
        for p in ex_ends:
         self.exon_ends[chr][bin][variant_id].append(p)
    
    f.close()
    
    sys.stderr.write("Finished loading genes.\n")


  def annotate_location(self, location):
    if not location.chr in self.gene_starts:
      return []
        
    mid_point = int((location.start + location.end) / 2)
    
    start_bin = int(location.start / self.bin_size)
    end_bin = int(location.end / self.bin_size)  
    
    # first find all the variants we might belong to
    closest_variants = collections.defaultdict()
    closest_d = collections.defaultdict(int)
    closest_abs_d = collections.defaultdict(int)
    
    
     
    for bin in range(start_bin, end_bin + 1):
      if not bin in self.gene_starts[location.chr]:
        continue
      
     
      for variant_id in self.gene_starts[location.chr][bin]:
        refseq = lib.sample.parse_id_from_variant(variant_id)
        
        entrez = self.entrezes[refseq]
        
        strand = self.gene_strands[location.chr][bin][variant_id]
        
        gene_start = self.gene_starts[location.chr][bin][variant_id]
        gene_end = self.gene_ends[location.chr][bin][variant_id]
        
        #
        # Deal with a peak being in a promoter
        #
        
        promoter_start = self.promoter_starts[location.chr][bin][variant_id]
        promoter_end = self.promoter_ends[location.chr][bin][variant_id]
        
        absd = -1
        
        if strand == "+":
          # for short genes we need to check the max distance of
          # gene end or promoter end          
          
          if mid_point >= promoter_start and mid_point <= gene_end:
            d = mid_point - gene_start
            absd = abs(d)
        else:
          # negative strand
          if mid_point >= gene_start and mid_point <= promoter_end:
            d = gene_end - mid_point
            absd = abs(d)
           
        if absd != -1:
          # we are in a variant so record the closest in the entrez group
          if entrez in closest_abs_d:
            if absd < closest_abs_d[entrez]: # closest_d[entrez]
              closest_variants[entrez] = variant_id
              closest_d[entrez] = d
              closest_abs_d[entrez] = absd
          else:
            closest_variants[entrez] = variant_id
            closest_d[entrez] = d
            closest_abs_d[entrez] = absd

    
    # Now classify each separate closest variant (which must all have a different entrez id)

    genes = []
    
    for entrez in sorted(closest_variants):        
      variant_id = closest_variants[entrez]
      
      gene = self.refseq.get_classification(variant_id, mid_point)

      genes.append(gene)
    
    return genes
