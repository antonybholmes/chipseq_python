# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 15:55:15 2014

@author: Antony Holmes
"""

# -*- coding: utf-8 -*-
"""
Find the common overlaps between a list of files reported as a common region

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import collections
import re

import annotation

file = sys.argv[1]

promoter_offset_up = 5000
promoter_offset_down = 4000

bin_size = 10000

refseq_file = "/ifs/scratch/cancer/Lab_RDF/abh2138/references/ucsc/ucsc_refseq_exons_entrez_hg19.txt"
entrez_gene_map = annotation.load_entrez_gene_map(refseq_file)

affy_gene_expression = annotation.AffyGeneExpression()

rna_gene_expression = annotation.RnaSeqGeneExpression()


collapsed_genes = collections.defaultdict(str)
collapsed_locations = collections.defaultdict(list)
collapsed_p = collections.defaultdict(float)

annotation.gene_orient_peak_locations(file, \
  entrez_gene_map, \
  collapsed_genes, \
  collapsed_locations, \
  collapsed_p)
  
entrezes = sorted(collapsed_genes)

sys.stdout.write("entrez\tgene_symbol\tp\tlocation_count\tlocations\taffy_gc_vs_n_m_gene_exp\trna_seq_gc_vs_n_m_gene_exp\n")

for entrez in entrezes:
  sys.stdout.write(entrez);
  
  sys.stdout.write("\t" + collapsed_genes[entrez])
  sys.stdout.write("\t" + str(collapsed_p[entrez]))
  sys.stdout.write("\t" + str(len(collapsed_locations[entrez])))
  sys.stdout.write("\t" + ";".join(sorted(collapsed_locations[entrez])))

  sys.stdout.write("\t" + affy_gene_expression.get_expression(entrez))
  sys.stdout.write("\t" + rna_gene_expression.get_expression(entrez))
  
  sys.stdout.write("\n");