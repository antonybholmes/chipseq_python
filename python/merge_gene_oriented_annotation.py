# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 09:43:29 2014

@author: Antony Holmes
"""
import sys
import collections
import re

import peaks

# Take two gene oriented annotation files and merge them together
# also check if the genes are up or down-regulated in expression
# data

def load_gene_exp(file, genes_map, type):
  f = open(file, "r")
  
  for line in f:
    ls = line.strip()
    
    if len(ls) == 0:
      continue
    
    tokens = ls.split("\t");
    
    if len(tokens[0]) == 0 or re.match(r'n\/a', tokens[0]):
      continue
    
    genes_map[tokens[0]] = type
  
  f.close()

def parse_merged(file,
                 entrez_map,
                 peaks_map,
                 peaks_p_map,
                 peak_counts_map,
                 annotation_map):
                   
  sys.stderr.write(file + "\n")
  
  f = open(file, "r")

  # skip header
  f.readline()

  for line in f:
    # remove newlines etc before splitting
    tokens = line.strip().split("\t")
    gene = tokens[0]
  
    entrezes = re.split(r'[,;]', tokens[1])
    p = float(tokens[2])
    peaks = tokens[3].split(";")
  
    entrez_map[gene] = ";".join(sorted(entrezes))
    peaks_map[gene] = peaks
    peaks_p_map[gene] = p
    peak_counts_map[gene] = len(peaks)
  
    for entrez in entrezes:
      if entrez in affy_regulation_entrez:
        annotation_map[gene]["affy_reg"][affy_regulation_entrez[entrez]] = 1
    
      if entrez in affy_dz_regulation_entrez:
        annotation_map[gene]["affy_dz_reg"][affy_dz_regulation_entrez[entrez]] = 1
    
      if entrez in rna_regulation_entrez:
        annotation_map[gene]["rna_reg"][rna_regulation_entrez[entrez]] = 1

  f.close()

def test_ig_peaks(ig_peaks, peaks_map, annotation_map):
  for gene in peaks_map:
    for peak in peaks_map[gene]:
       match = re.search(r'(chr.+?):(\d+)-(\d+)', peak)
       
       chr = match.group(1)
       start = int(match.group(2))
       end = int(match.group(3))
       
       # see if we overlap a peak
       features = peaks.peak_overlap(ig_peaks, chr, start, end)
       
       if (len(features) == 0):
         continue
       
       annotation_map[gene] = ";".join(features)

# start of script 
 
file1 = sys.argv[1]
file2 = sys.argv[2]

n1 = file1
n1 = re.sub(r'\.txt', '', n1)
n1 = re.sub(r'_p', '_p_10E-', n1)
n1 = re.sub(r'gene_mir_', "", n1)

n2 = file2
n2 = re.sub(r'\.txt', '', n2)
n2 = re.sub(r'_p', '_p_10E-', n2)
n2 = re.sub(r'gene_mir_', "", n2)

#
# Load the Ig vs input peaks so we can test if any of these peaks
# overlap them
#

ig_peaks = peaks.parse_peaks("/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/rdf/20140616/hg19/align_2_mismatches/analysis/chipseeqer/p10E-5_fold2/TF_targets_CB4_Ig_RK026_vs_Input_RK029_p5.txt")


#
# First load affymetrix up or down regulated
#

affy_regulation_entrez = collections.defaultdict(str)

#load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/hg-u133_plus_2/hg-u133_plus_2_entrez_genes.txt",
#  affy_regulation_entrez,
#  "not_moving")
#
## now mark as down-regulated
#
#load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/hg-u133_plus_2/cb_vs_n_m_hg-u133_plus_2_entrez_genes_down_regulated.txt",
#  affy_regulation_entrez,
#  "down")
#
## now mark as up-regulated
#
#load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/hg-u133_plus_2/cb_vs_n_m_hg-u133_plus_2_entrez_genes_up_regulated.txt",
#  affy_regulation_entrez,
#  "up");

# log 2 version

load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/gene_expression/cb_vs_n_m_hg-u133_plus_2_entrez.txt",
  affy_regulation_entrez,
  "not_moving")

# now mark as down-regulated

load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/gene_expression/cb_vs_n_m_hg_u133_plus_2_down_entrez_p5_z2.txt",
  affy_regulation_entrez,
  "down")

# now mark as up-regulated

load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/gene_expression/cb_vs_n_m_hg_u133_plus_2_up_entrez_p5_z2.txt",
  affy_regulation_entrez,
  "up");

  

#
# dz up/down
#

affy_dz_regulation_entrez = collections.defaultdict(str)

load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/dark_zone_vs_light_zone/dark_zone_vs_light_zone_entrez_genes.txt",
  affy_dz_regulation_entrez,
  "not_moving")

load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/dark_zone_vs_light_zone/dark_zone_vs_light_zone_entrez_down_regulated.txt",
  affy_dz_regulation_entrez,
  "down")


load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/dark_zone_vs_light_zone/dark_zone_vs_light_zone_entrez_up_regulated.txt",
  affy_dz_regulation_entrez,
  "up")



#
# Load rna-seq annotation
#

rna_regulation_entrez = collections.defaultdict(str)

#load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/rna_seq_entrez_genes.txt",
#  rna_regulation_entrez,
#  "not_moving");
#  
#load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/cb_vs_n_m_rna_seq_entrez_genes_down_regulated.txt",
#  rna_regulation_entrez,
#  "down");
#
#load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq/cb_vs_n_m_rna_seq_entrez_genes_up_regulated.txt",
#  rna_regulation_entrez,
#  "up");

# newer version using log transformed data

load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/gene_expression/rna_seq_fpkm_expression_entrez.txt",
  rna_regulation_entrez,
  "not_moving");
  
load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/gene_expression/rna_seq_fpkm_expression_cb_m_n_down_entrez_p5_z2.txt",
  rna_regulation_entrez,
  "down");

load_gene_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/gene_expression/rna_seq_fpkm_expression_cb_m_n_up_entrez_p5_z2.txt",
  rna_regulation_entrez,
  "up");

#
# Now annotate
#


entrez_map = collections.defaultdict(str)
annotation_map = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(int)))

peaks_map_1 = collections.defaultdict(str)
peaks_p_map_1 = collections.defaultdict(float)
peak_counts_map_1 = collections.defaultdict(int)


parse_merged(file1,
             entrez_map,
             peaks_map_1,
             peaks_p_map_1,
             peak_counts_map_1,
             annotation_map)

ig_map_1 = collections.defaultdict(str)             
test_ig_peaks(ig_peaks, peaks_map_1, ig_map_1)
             
peaks_map_2 = collections.defaultdict(str)
peaks_p_map_2 = collections.defaultdict(float)
peak_counts_map_2 = collections.defaultdict(int)

parse_merged(file2,
             entrez_map,
             peaks_map_2,
             peaks_p_map_2,
             peak_counts_map_2,
             annotation_map)

ig_map_2 = collections.defaultdict(str)             
test_ig_peaks(ig_peaks, peaks_map_2, ig_map_2)



print("gene\tentrez_id\taffy_gc_gene_exp\trna-seq_gc_gene_exp\taffy_dz_gene_exp\tgroup_1\tgroup_2\tlog10_p-value_group_1\tpeak_counts_group_1\tpeaks_group_1\tig_peak_1\tlog10_p-value_group_2\tpeak_counts_group_2\tpeaks_group_2\tig_peak_2");

for gene in sorted(entrez_map):
  sys.stdout.write("\t".join([gene, entrez_map[gene]]) + "\t")
  
  if "affy_reg" in annotation_map[gene]:
    sys.stdout.write(";".join(sorted(annotation_map[gene]["affy_reg"])))
  else:
    sys.stdout.write("n/a")
  
  sys.stdout.write("\t")
  
  if "rna_reg" in annotation_map[gene]:
    sys.stdout.write(";".join(sorted(annotation_map[gene]["rna_reg"])))
  else:
    sys.stdout.write("n/a")
  
  sys.stdout.write("\t")
  
  if "affy_dz_reg" in annotation_map[gene]:
    sys.stdout.write(";".join(sorted(annotation_map[gene]["affy_dz_reg"])))
  else:
    sys.stdout.write("n/a")
  
  sys.stdout.write("".join(["\t", n1, "\t", n2, "\t"]))
  
  if gene in peaks_map_1:
    sys.stdout.write("\t".join([str(peaks_p_map_1[gene]), str(peak_counts_map_1[gene]), ";".join(peaks_map_1[gene])]))
  else:
    sys.stdout.write("n/a\t0\tn/a")
  
  sys.stdout.write("\t")
  
  if gene in ig_map_1:
    sys.stdout.write(ig_map_1[gene])
  else:
    sys.stdout.write("n/a")
  
  sys.stdout.write("\t")
  
  if gene in peaks_map_2:
    sys.stdout.write("\t".join([str(peaks_p_map_2[gene]), str(peak_counts_map_2[gene]), ";".join(peaks_map_2[gene])]))
  else:
    sys.stdout.write("n/a\t0\tn/a")
  
  sys.stdout.write("\t")
  
  if gene in ig_map_2:
    sys.stdout.write(ig_map_2[gene])
  else:
    sys.stdout.write("n/a")
  
  sys.stdout.write("\n")