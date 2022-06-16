# -*- coding: utf-8 -*-
"""
Classes to deal with finding the TSS closest to a point

Created on Mon Aug 25 17:25:56 2014

@author: Antony Holmes
"""

from . import mouse
from .. import tss


class RefSeqTssClassification(tss.RefSeqTssClassification):
  def __init__(self, promoter_offset_up, promoter_offset_down):
    super().__init__(mouse.REFSEQ_FILE, promoter_offset_up, promoter_offset_down)


class RefSeqTssStart(tss.RefSeqTss):
  def __init__(self, promoter_offset_up, promoter_offset_down):
    super().__init__(mouse.REFSEQ_FILE, promoter_offset_up, promoter_offset_down)
    

class RefSeqEnd(tss.RefSeqEnd):
  def __init__(self, promoter_offset_up, promoter_offset_down):
    super().__init__(mouse.REFSEQ_FILE, promoter_offset_up, promoter_offset_down)
   
   
class RefSeqAnnotation(tss.RefSeqAnnotation):
  def __init__(self, promoter_offset_up, promoter_offset_down, bin_size):
    super().__init__(mouse.REFSEQ_FILE, promoter_offset_up, promoter_offset_down, bin_size)
