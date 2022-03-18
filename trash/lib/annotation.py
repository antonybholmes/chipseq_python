# -*- coding: utf-8 -*-
"""
Annotation functions

Created on Thu Jun 26 09:43:29 2014

@author: Antony Holmes
"""
import sys
import collections
import re
#import radix_tree


REFSEQ_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/ucsc_refseq_exons_entrez_hg19_20150217.txt"
#REFSEQ_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/ucsc_refseq_exons_entrez_hg19.txt"
#REFSEQ_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/ucsc_refseq_exons_entrez_canonical_only_hg19.txt"

#ENSEMBL_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/ucsc_ensembl_exons_gene_symbol_hg19.txt"

#RDF_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/rdf_ucsc_refseq_ensembl_genes_hg19.txt"

# refseqs only
#RDF_FILE = "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/rdf_ucsc_refseq_genes_hg19.txt"



  


  
  
def get_mir_id(id):
  mir = id.lower()
   
  #  if not re.match('(hsa)-(.+?)-(.+?).*', s):
  #    return s
  #    
  #  matcher = re.match('(hsa)-([^\-]+)-([^\-]+)', s)
  #    
  #  mir = matcher.group(1) + "-" + matcher.group(2) + "-" + matcher.group(3)

  return mir

def get_mir_precursor_id(id):
  mir = id.lower()
   
  if not re.match('(hsa)-(.+?)-(.+?).*', mir):
    return mir
    
  matcher = re.match('(hsa)-([^\-]+)-([^\-]+)', mir)
  mir = matcher.group(1) + "-" + matcher.group(2) + "-" + matcher.group(3)

  return mir
  
  
def gene_orient_peak_locations(file, \
  entrez_gene_map, \
  collapsed_genes, \
  collapsed_locations, \
  collapsed_p):
  f = open(file, 'r')

  # skip header
  header = f.readline().split("\t")
  
  for line in f:
    ls = line.strip()
    
    if len(ls) == 0:
      continue
    
    tokens = line.split("\t")
    
    id = tokens[0]

    location = tokens[1]
    p = float(tokens[3])
    
    entrezes = tokens[find_heading_index(header, "entrez")].split(";")
    
    for entrez in entrezes:
      if entrez == "n/a":
        continue
    
      if entrez in entrez_gene_map:
        collapsed_genes[entrez] = entrez_gene_map[entrez]
      else:
        collapsed_genes[entrez] = "n/a"
        
      collapsed_locations[entrez].append(location)
        
      if entrez in collapsed_p:
        if p < collapsed_p[entrez]:
            collapsed_p[entrez] = p
      else:
        collapsed_p[entrez] = p
        
        
  f.close()
  


  



def unique_symbols(line):
  tokens = line.split(",")
  
  s = set()
  
  for item in tokens:
    s.add(item)
    
  return ",".join(sorted(s))
  
  
class EntrezGeneLookup(object):
  """
  Loads a mapping between entrez ids and gene symbols.
  """
  
  def __init__(self):
    self.entrez_gene_map = collections.defaultdict(str)
 
    sys.stderr.write("Loading genes from " + REFSEQ_FILE + "...\n")
  
    f = open(REFSEQ_FILE, 'r')

    # skip header
    f.readline()
  
    for line in f:
      line = line.strip()
    
      if len(line) == 0:
        continue
    
      tokens = line.split("\t")
    
      id = tokens[0]
      symbol = unique_symbols(tokens[1])
    
      if symbol == "n/a":
        continue
    
      self.entrez_gene_map[id] = symbol
    
    f.close()
  
    sys.stderr.write("Finished loading genes.\n")
  
  def get_symbol(self, entrez):
    if entrez not in self.entrez_gene_map:
      return "n/a"
      
    return self.entrez_gene_map[entrez]
    

    
class MirPreviousIds(object):
  def __init__(self):
    self.previous_ids = collections.defaultdict(lambda: collections.defaultdict(bool))
    
    f = open("/ifs/scratch/cancer/Lab_RDF/ngs/references/mirbase/19/mirna_mature.txt", "r")

    # skip header
    f.readline()

    for line in f:
      line = line.strip()
    
      if len(line) == 0:
        continue
  
      tokens = line.split("\t")
  
      mir = tokens[1]

      if not re.match(r'hsa.*', mir):
        continue
      
      previous = tokens[2]
      
      if len(previous) < 1:
        continue
      
      old_ids = previous.split(";")
      
      for id in old_ids:
        #sys.stderr.write(get_mir_id(mir) + " " + get_mir_id(id) + "\n")
        self.previous_ids[get_mir_id(mir)][get_mir_id(id)] = True

    f.close()
    
  def get_previous_ids(self, id):
    id = get_mir_id(id)
    
    ret = []
    
    if id in self.previous_ids:
      for old_id in self.previous_ids[id]:
        ret.append(old_id)
    
    return ret


class MirAnnotation:
  """
  Loads a mir database.
  """
  
  def __init__(self, bin_size):
    self.mir_types = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(str)))
    self.mir_starts = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
    self.mir_ends = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
    self.mir_strands = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(str)))
    self.bin_size = bin_size
    
    sys.stderr.write("Loading miRs...\n")
    
    f = open("/ifs/scratch/cancer/Lab_RDF/ngs/references/mirbase/19/mirbase_19_hg19_genomic_locations.txt", "r")

    # skip header
    f.readline()

    for line in f:
      line = line.strip()
    
      if len(line) == 0:
        continue
  
      tokens = line.split("\t")
  
      name = tokens[0]
      group = tokens[1]
      type = tokens[2]
      chrom = tokens[3]
      start = int(tokens[4])
      end = int(tokens[5])
      strand = tokens[6]
      
      #if type != "miRNA_primary_transcript":
      #  continue

      # only annotate the mature
      if type != "miRNA":
        continue
  
      start_bin = int(start / bin_size)
      end_bin = int(end / bin_size)
  
      for bin in range(start_bin, end_bin + 1):
        self.mir_types[chrom][bin][name] = type
        self.mir_starts[chrom][bin][name] = start
        self.mir_ends[chrom][bin][name] = end
        self.mir_strands[chrom][bin][name] = strand

    f.close()
  
  def annotate_location(self, location, max_gap, annotation_types, mss):
    self.annotate_position(location.chr, \
      location.start, \
      location.end, \
      max_gap, \
      annotation_types, \
      mss)
    
  def annotate_position(self, chr, start, end, max_gap, annotation_types, mss):
    if not chr in self.mir_types:
      return
    
    mid_point = int((start + end) / 2)
  
    start_bin = int((start - max_gap) / self.bin_size)
    end_bin = int((end + max_gap) / self.bin_size)

    #sys.stderr.write("mir " + chr + " " + str(mid_point) + "\n")    
    
    for bin in range(start_bin, end_bin + 1):
      if not bin in self.mir_starts[chr]:
        continue
      
      for name in self.mir_starts[chr][bin]:
        mir_start = self.mir_starts[chr][bin][name]
        mir_end = self.mir_ends[chr][bin][name]
        mir_type = self.mir_types[chr][bin][name]
        mir_strand = self.mir_strands[chr][bin][name]
        
        found = True
        
        if mir_strand == "+":
          if mir_start - mid_point <= max_gap and mir_start - mid_point > 0:
            annotation_types[name] = ",".join([mir_type, "before"])
          elif mid_point - mir_end <= max_gap and mid_point - mir_end > 0:
            annotation_types[name] = ",".join([mir_type, "after"])
          elif mid_point >= mir_start and mid_point <= mir_end:
            annotation_types[name] = ",".join([mir_type, "over"])
          else:
            found = False
        
          if found == True:
            mss[name] = str(mid_point - mir_start)
          else:
            mss[name] = "n/a"
        else:
          # neg strand
          if mir_start - mid_point <= max_gap and mir_start - mid_point > 0:
            annotation_types[name] = ",".join([mir_type, "after"])
          elif mid_point - mir_end <= max_gap and mid_point - mir_end > 0:
            annotation_types[name] = ",".join([mir_type, "before"])
          elif mid_point >= mir_start and mid_point <= mir_end:
            annotation_types[name] = ",".join([mir_type, "over"])
          else:
            found = False
        
          if found == True:
            mss[name] = str(mir_end - mid_point)
          else:
            mss[name] = "n/a"


class MirTranscriptAnnotation:
  """
  Can be used to check if a gene contains a mir or not
  """
  
  def __init__(self):
    self.precursors = collections.defaultdict(lambda: collections.defaultdict(bool))
    self.mature = collections.defaultdict(lambda: collections.defaultdict(bool))
    
    
    sys.stderr.write("Loading Transcript miRs...\n")
    
    f = open("/ifs/scratch/cancer/Lab_RDF/ngs/references/mirbase/19/mirbase_19_transcript_locations_hg19_entrez.txt", "r")

    # skip header
    f.readline()

    for line in f:
      line = line.strip()
    
      if len(line) == 0:
        continue
  
      tokens = line.split("\t")
  
      mir = tokens[0]
      entrez = tokens[6]
      
      if entrez == "n/a":
        continue
      
      self.precursors[entrez][mir] = True

    f.close()

    #    
    # open a mapping between precursor and mature
    #
    
    f = open("/ifs/scratch/cancer/Lab_RDF/ngs/references/mirbase/19/hsa.gff3", "r")


    id_precursor_map = collections.defaultdict(str)    
    
    for line in f:
      if re.match(r'^#.*', line):
        continue
      
      line = line.strip()
    
      if len(line) == 0:
        continue
  
      tokens = line.split("\t")
      
      type = tokens[2]
      
      matcher = re.match(r'^ID=([^\;]+).*Name=([^\;]+).*', tokens[8])
      
      id = matcher.group(1)
      mir = matcher.group(2)
      
      if type == "miRNA_primary_transcript":
        id_precursor_map[id] = mir
      else:
        matcher = re.match(r'.*derives_from=([^\;]+).*', tokens[8])
        
        precursor_id = matcher.group(1)
        
        precursor = id_precursor_map[precursor_id]
        
        self.mature[precursor][mir] = True

    f.close()
    
  
  def find_mirs(self, entrez):
    ret = []
    
    if entrez not in self.precursors:
      return ret
      
    for precursor in sorted(self.precursors[entrez]):
      for mature in self.mature[precursor]:
        ret.append(mature)
      
    return ret
    
    
    
    

