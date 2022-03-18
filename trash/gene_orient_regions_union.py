# -*- coding: utf-8 -*-
"""
Display regions in a gene oriented fashion

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re

import lib.human.genes
import lib.text
import lib.headings
import lib.sample
import lib.text


def get_sample_id(text):
  return re.match(r'^.+ (.+)', text).group(1)
  
def two_way(file):
  gene_orientation = lib.human.genes.GeneOrientatedPeaks("Region")
  
  # add some custom annotations for katia for testing etc
  
  gene_orientation.add_expression("RNA-seq hg19 GCvsN TPM", \
    lib.expression.GeneExpression("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/ucsc_refseq_hg19_20171207/tpm_cb_vs_n/"))
  
  gene_orientation.add_expression("RNA-seq hg19 GCvsM TPM", \
    lib.expression.GeneExpression("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/ucsc_refseq_hg19_20171207/tpm_cb_vs_m/"))
  
  gene_orientation.add_expression("RNA-seq hg19 GCvsN rmdup TPM", \
    lib.expression.GeneExpression("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/ucsc_refseq_hg19_20171207/rmdup/tpm_cb_vs_n/"))
  
  gene_orientation.add_expression("RNA-seq hg19 GCvsM rmdup TPM", \
    lib.expression.GeneExpression("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/ucsc_refseq_hg19_20171207/rmdup/tpm_cb_vs_m/"))
   
  gene_orientation.load_annotations(file)
  
  f = open(file, 'r')

  # skip header
  header = f.readline().strip().split("\t")
  
  sample_column = lib.text.find_index(header, lib.headings.SAMPLE)
  
  group1 = get_sample_id(header[sample_column])
  group2 = get_sample_id(header[sample_column + 1])
  
  id_column = lib.text.find_index(header, lib.headings.REFSEQ_ID)
  mir_column = lib.text.find_index(header, lib.headings.MIR_SYMBOL)
  
  groups1 = set()
  groups2 = set()
  groupsi = set()
  
  for line in f:
    ls = line.strip()
    
    if len(ls) == 0:
      continue
    
    tokens = ls.split("\t")
    
    id = tokens[id_column]
 
    g1 = tokens[sample_column]
    g2 = tokens[sample_column + 1]
    
    if id != lib.text.NA:
      if g1 != lib.text.NA and g2 == lib.text.NA:
        groups1.add(id)
      elif g1 == lib.text.NA and g2 != lib.text.NA:
        groups2.add(id)
      elif g1 != lib.text.NA and g2 != lib.text.NA:
        groupsi.add(id)
      else:
        pass
        
    mir = tokens[mir_column]
    
    if mir != lib.text.NA:
      if g1 != lib.text.NA and g2 == lib.text.NA:
        groups1.add(mir)
      elif g1 == lib.text.NA and g2 != lib.text.NA:
        groups2.add(mir)
      elif g1 != lib.text.NA and g2 != lib.text.NA:
        groupsi.add(mir)
      else:
        pass
      
    
  f.close()
  
  # Write custom header

  sys.stdout.write(group1 + " Unique");
  sys.stdout.write("\tOverlap");
  sys.stdout.write("\t" + group2 + " Unique\t");

  gene_orientation.print_header()  
  
  for id in gene_orientation.get_ids():
    u1 = 0
    u2 = 0
    inters = 0
    
    if id in groupsi:
      inters = 1
    
    if id in groups1:
      u1 = 1
      
    if id in groups2:
      u2 = 1
    
    sys.stdout.write("\t".join([str(u1), str(inters), str(u2)]) + "\t")
    
    gene_orientation.gene_orient_peak(id)
  
  for mir in gene_orientation.get_mirs():
    u1 = 0
    u2 = 0
    inters = 0
    
    if mir in groupsi:
      inters = 2
    
    if mir in groups1:
      u1 = 2
      
    if mir in groups2:
      u2 = 2
    
    sys.stdout.write("\t".join([str(u1), str(inters), str(u2)]) + "\t")
    
    gene_orientation.mir_orient_peak(mir)
    
def three_way(file):
  gene_orientation = lib.human.genes.GeneOrientatedPeaks("Region")
  gene_orientation.load_annotations(file)
  
  f = open(file, 'r')

  # skip header
  header = f.readline().strip().split("\t")
  
  sample_column = lib.text.find_index(header, lib.headings.SAMPLE)
  
  group1 = get_sample_id(header[sample_column])
  group2 = get_sample_id(header[sample_column + 1])
  group3 = get_sample_id(header[sample_column + 2])
  
  sys.stderr.write("groups " + group1 + " " + group2 + " " + group3 + "\n")
  
  id_column = lib.text.find_index(header, lib.headings.REFSEQ_ID)
  mir_column = lib.text.find_index(header, lib.headings.MIR_SYMBOL)
  
  g1 = set()
  g2 = set()
  g3 = set()
  g12 = set()
  g13 = set()
  g23 = set()
  gi = set()
  
  for line in f:
    ls = line.strip()
    
    if len(ls) == 0:
      continue
    
    tokens = ls.split("\t")
    
    id = tokens[id_column]
 
    p1 = tokens[sample_column]
    p2 = tokens[sample_column + 1]
    p3 = tokens[sample_column + 2]    
    
    if id != lib.text.NA:
      if p1 != lib.text.NA and p2 == lib.text.NA and p3 == lib.text.NA:
        g1.add(id)
      
      if p1 == lib.text.NA and p2 != lib.text.NA and p3 == lib.text.NA:
        g2.add(id)
      
      if p1 == lib.text.NA and p2 == lib.text.NA and p3 != lib.text.NA:
        g3.add(id)
      
      if p1 != lib.text.NA and p2 != lib.text.NA:
        g12.add(id)

      if p1 != lib.text.NA and p3 != lib.text.NA:
        g13.add(id)
      
      if p2 != lib.text.NA and p3 != lib.text.NA:
        g23.add(id)
        
      if p1 != lib.text.NA and p2 != lib.text.NA and p3 != lib.text.NA:
        gi.add(id)
        
    mir = tokens[mir_column]
    
    if mir != lib.text.NA:
      if p1 != lib.text.NA and p2 == lib.text.NA and p3 == lib.text.NA:
        g1.add(mir)
      
      if p1 == lib.text.NA and p2 != lib.text.NA and p3 == lib.text.NA:
        g2.add(mir)
      
      if p1 == lib.text.NA and p2 == lib.text.NA and p3 != lib.text.NA:
        g3.add(mir)
      
      if p1 != lib.text.NA and p3 != lib.text.NA:
        g13.add(mir)
      
      if p2 != lib.text.NA and p3 != lib.text.NA:
        g23.add(mir)
      
      if p1 != lib.text.NA and p2 != lib.text.NA:
        g12.add(mir)

      if p1 != lib.text.NA and p3 != lib.text.NA:
        g13.add(mir)
        
      if p1 != lib.text.NA and p2 != lib.text.NA and p3 != lib.text.NA:
        gi.add(mir)
        
    
  f.close()
  
  # Write custom header
  
  rk1 = lib.sample.get_sample_id(group1)
  rk2 = lib.sample.get_sample_id(group2)
  rk3 = lib.sample.get_sample_id(group3)

  sys.stdout.write(group1 + " Unique\t");
  sys.stdout.write(group2 + " Unique\t");
  sys.stdout.write(group3 + " Unique\t");
  sys.stdout.write(rk1 + " " + rk2 + " Overlap\t");
  sys.stdout.write(rk1 + " " + rk3 + " Overlap\t");
  sys.stdout.write(rk2 + " " + rk3 + " Overlap\t");
  sys.stdout.write("Overlap\t");

  gene_orientation.print_header()  
  
  for id in gene_orientation.get_ids():
    u1 = 0
    u2 = 0
    u3 = 0
    u12 = 0
    u13 = 0
    u23 = 0
    inters = 0
    
    if id in gi:
      inters = 1
    
    if id in g1:
      u1 = 1
      
    if id in g2:
      u2 = 1

    if id in g3:
      u3 = 1
    
    if id in g12:
      u12 = 1
      
    if id in g13:
      u13 = 1
    
    if id in g23:
      u23 = 1
      
    sys.stdout.write("\t".join([str(u1), str(u2), str(u3), str(u12), str(u13), str(u23), str(inters)]) + "\t")
    
    gene_orientation.gene_orient_peak(id)
  
  for mir in gene_orientation.get_mirs():
    u1 = 0
    u2 = 0
    u3 = 0
    u12 = 0
    u13 = 0
    u23 = 0
    inters = 0
    
    if mir in gi:
      inters = 2
    
    if mir in g1:
      u1 = 2
      
    if mir in g2:
      u2 = 2

    if mir in g3:
      u3 = 2
    
    if mir in g12:
      u12 = 2
      
    if mir in g13:
      u13 = 2
    
    if mir in g23:
      u23 = 2
    
    sys.stdout.write("\t".join([str(u1), str(u2), str(u3), str(u12), str(u13), str(u23), str(inters)]) + "\t")
    
    gene_orientation.mir_orient_peak(mir)
    
    
def gene_orient_peaks(file):
  f = open(file, 'r')

  # skip header
  header = f.readline().strip().split("\t")
  f.close()
  
  sample_column = lib.text.find_index(header, lib.headings.SAMPLE)
  
  group3 = None
  
  if lib.headings.SAMPLE in header[sample_column + 2]:
    group3 = get_sample_id(header[sample_column + 2])
    
  if group3 is not None:
    three_way(file)
  else:
    two_way(file)
  
    
file = sys.argv[1]

if len(sys.argv) > 2:
  additional_annotations = sys.argv[2]
else:
  additional_annotations = ""

gene_orient_peaks(file)
