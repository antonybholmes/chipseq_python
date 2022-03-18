# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 15:12:34 2014

@author: Antony Holmes
"""

import sys
import collections

import pychipseq.expression
import pychipseq.annotation
import pychipseq.headings


class RnaSeqGeneICNvsGFPExpression(pychipseq.expression.RnaSeqExpression):
  def __init__(self):
    super().__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq_icn_vs_gfp/")
    

class CustomICNvsGFPExpression(pychipseq.expression.CustomExpression):
  """
  Annotate up or down in MO cells.
  """
  def __init__(self, type):
    super().__init__("RNAseq ICNvsGFP", RnaSeqGeneICNvsGFPExpression())
