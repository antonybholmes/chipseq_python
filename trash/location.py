# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 17:37:34 2014

@author: Antony Holmes
"""

import sys

import tss


location = sys.argv[1]

annotation = tss.RefSeqTss(5000, 4000)

gene = annotation.get_tss_from_location(location)

sys.stdout.write(gene.entrez + " " + str(gene.d) + "\n")

ref = tss.RefSeqAnnotation(5000, 4000, 10000)

genes = ref.annotate_location(location)

for gene in genes:
  sys.stdout.write(gene.variant_id + " " + ":".join(gene.types) + "\n")