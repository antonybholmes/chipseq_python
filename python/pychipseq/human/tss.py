# -*- coding: utf-8 -*-
"""
Classes to deal with finding the TSS closest to a point

@author: Antony Holmes
"""

import pychipseq.genomic
import pychipseq.sample
import pychipseq.human.annotation
import pychipseq.tss
import pychipseq.genes

class RefSeqTssClassification(pychipseq.tss.RefSeqTssClassification):
  def __init__(self, prom_ext_5p, prom_ext_3p):
    super().__init__(pychipseq.human.annotation.REFSEQ_FILE, prom_ext_5p, prom_ext_3p)


class RefSeqStart(pychipseq.tss.RefSeqStart):
  def __init__(self, refseq_genes: pychipseq.genes.RefSeqGenes, prom_ext_5p, prom_ext_3p):
    super().__init__(pychipseq.human.annotation.REFSEQ_FILE, refseq_genes, prom_ext_5p, prom_ext_3p)
   
    
class RefSeqEnd(pychipseq.tss.RefSeqEnd):
  def __init__(self, refseq_genes: pychipseq.genes.RefSeqGenes, prom_ext_5p, prom_ext_3p):
    super().__init__(pychipseq.human.annotation.REFSEQ_FILE, refseq_genes, prom_ext_5p, prom_ext_3p)
   
   
class RefSeqAnnotation(pychipseq.tss.RefSeqAnnotation):
  def __init__(self, prom_ext_5p, prom_ext_3p, bin_size):
    super().__init__(pychipseq.human.annotation.REFSEQ_FILE, prom_ext_5p, prom_ext_3p, bin_size)


class OverlapTss(pychipseq.tss.OverlapTss):
  def __init__(self, block_size=100):
    super().__init__(pychipseq.human.annotation.REFSEQ_FILE, block_size)