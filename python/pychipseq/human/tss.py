# -*- coding: utf-8 -*-
"""
Classes to deal with finding the TSS closest to a point.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.
This program is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with 
this program. If not, see <https://www.gnu.org/licenses/>. 

Copyright (C) 2022 Antony Holmes.
"""

import collections
from . import human
from .. import tss
from .. import genes


class RefSeqTssClassification:
    @staticmethod
    def getInstance(prom_ext_5p: int, prom_ext_3p: int):
        return tss.RefSeqTssClassification.getInstance(human.REFSEQ_FILE, prom_ext_5p, prom_ext_3p)

    def __init__(self):  # , prom_ext_5p, prom_ext_3p):
        # super().__init__(annotation.REFSEQ_FILE, prom_ext_5p, prom_ext_3p)
        raise Exception("Call RefSeqTssClassification.getInstance()")


class RefSeqStart(tss.RefSeqStart):
    def __init__(self, refseq_genes: genes.RefSeqGenes, prom_ext_5p, prom_ext_3p):
        super().__init__(human.REFSEQ_FILE, refseq_genes, prom_ext_5p, prom_ext_3p)


class RefSeqEnd(tss.RefSeqEnd):
    def __init__(self, refseq_genes: genes.RefSeqGenes, prom_ext_5p, prom_ext_3p):
        super().__init__(human.REFSEQ_FILE, refseq_genes, prom_ext_5p, prom_ext_3p)


# class RefSeqAnnotation(tss.RefSeqAnnotation):
#   def __init__(self, prom_ext_5p, prom_ext_3p, bin_size):
#     super().__init__(annotation.REFSEQ_FILE, prom_ext_5p, prom_ext_3p, bin_size)


class RefSeqAnnotationFactory:
    _classifiers = collections.defaultdict(
        lambda: collections.defaultdict(lambda: {}))

    @staticmethod
    def getInstance(prom_ext_5p: int, prom_ext_3p: int, bin_size: int = 1000):
        if bin_size not in RefSeqAnnotationFactory._classifiers[prom_ext_5p][prom_ext_3p]:
            RefSeqAnnotationFactory._classifiers[prom_ext_5p][prom_ext_3p][bin_size] = tss.RefSeqAnnotationFactory.getInstance(
                human.REFSEQ_FILE, prom_ext_5p, prom_ext_3p, bin_size)

        return RefSeqAnnotationFactory._classifiers[prom_ext_5p][prom_ext_3p][bin_size]

    def __init__(self):
        raise Exception("Call RefSeqAnnotation.getInstance()")


class OverlapTss(tss.OverlapTss):
    def __init__(self, block_size=100):
        super().__init__(human.REFSEQ_FILE, block_size)
