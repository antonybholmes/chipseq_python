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
import lib.mouse.annotation
import lib.tss


class RefSeqTssClassification(lib.tss.RefSeqTssClassification):
  def __init__(self, promoter_offset_up, promoter_offset_down):
    super(RefSeqTssClassification, self).__init__(lib.mouse.annotation.REFSEQ_FILE, promoter_offset_up, promoter_offset_down)


class RefSeqTss(lib.tss.RefSeqTss):
  def __init__(self, promoter_offset_up, promoter_offset_down):
    super(RefSeqTss, self).__init__(lib.mouse.annotation.REFSEQ_FILE, promoter_offset_up, promoter_offset_down)
    

class RefSeqEnd(lib.tss.RefSeqEnd):
  def __init__(self, promoter_offset_up, promoter_offset_down):
    super(RefSeqEnd, self).__init__(lib.mouse.annotation.REFSEQ_FILE, promoter_offset_up, promoter_offset_down)
   
   
class RefSeqAnnotation(lib.tss.RefSeqAnnotation):
  def __init__(self, promoter_offset_up, promoter_offset_down, bin_size):
    super(RefSeqAnnotation, self).__init__(lib.mouse.annotation.REFSEQ_FILE, promoter_offset_up, promoter_offset_down, bin_size)
