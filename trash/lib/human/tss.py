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
import lib.tss


class RefSeqTssClassification(lib.tss.RefSeqTssClassification):
  def __init__(self, prom_ext_5p, prom_ext_3p):
    super(RefSeqTssClassification, self).__init__(lib.annotation.REFSEQ_FILE, prom_ext_5p, prom_ext_3p)

class RefSeqTss(lib.tss.RefSeqTss):
  def __init__(self, prom_ext_5p, prom_ext_3p):
    super(RefSeqTss, self).__init__(lib.annotation.REFSEQ_FILE, prom_ext_5p, prom_ext_3p)
   
    
class RefSeqEnd(lib.tss.RefSeqEnd):
  def __init__(self, prom_ext_5p, prom_ext_3p):
    super(RefSeqEnd, self).__init__(lib.annotation.REFSEQ_FILE, prom_ext_5p, prom_ext_3p)
   
   
class RefSeqAnnotation(lib.tss.RefSeqAnnotation):
  def __init__(self, prom_ext_5p, prom_ext_3p, bin_size):
    super(RefSeqAnnotation, self).__init__(lib.annotation.REFSEQ_FILE, prom_ext_5p, prom_ext_3p, bin_size)
