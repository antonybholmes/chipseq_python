# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 15:12:34 2014

@author: Antony Holmes
"""

import sys
import collections

import lib.expression
import lib.annotation
import lib.headings


class RnaSeqGeneICNvsGFPExpression(lib.expression.RnaSeqExpression):
  def __init__(self):
    super(RnaSeqGeneICNvsGFPExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq_icn_vs_gfp/")
    

class CustomICNvsGFPExpression(lib.expression.CustomExpression):
  """
  Annotate up or down in MO cells.
  """
  def __init__(self, type):
    super(RnaSeqGeneICNvsGFPExpression, self).__init__("RNAseq ICNvsGFP", RnaSeqGeneICNvsGFPExpression())
