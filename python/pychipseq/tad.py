# -*- coding: utf-8 -*-
"""
Classes to deal with finding the TSS closest to a point

Created on Mon Aug 25 17:25:56 2014

@author: Antony Holmes
"""

import collections
import sys
import re
from typing import Any, List, Mapping, Union
import numpy as np

import pychipseq.genomic
import pychipseq.text
import pychipseq.headings

TADS_HEADING = 'TAD domains'
TAD_HEADING = 'Genes in same TAD domain'
IS_TAD_HEADING = 'Gene in TAD'
IS_CLOSEST_TAD_HEADING = 'Closest gene in TAD'

class TADAnnotation(pychipseq.genomic.Annotation):
    """
    Annotate peaks for all possible refseq genes they might overlap.
    """

    def __init__(self, file, bin_size: int = 100):
        """
        Create a new pychipseq.annotation object.

        @param prom_ext_5p   How far upstream should be considered 
                             a promoter.
        @param prom_ext_3p   How far downstream should be considered 
                             a promoter.
        @param bin_size      The size of the bins to group genes into.
        """

        self._bin_size = bin_size

        self._tads = collections.defaultdict(list)
        self._tad_bins = collections.defaultdict(
            lambda: collections.defaultdict(int))

        print(f'Loading TADs from {file}...', file=sys.stderr)

        f = open(file, 'r')

        # skip header
        f.readline()

        for line in f:
            line = line.strip()

            tokens = line.split("\t")

            chr = tokens[0]
            start = int(tokens[1]) + 1
            end = int(tokens[2])
            location = pychipseq.genomic.Location(chr, start, end)

            ids = tokens[3]
            gene_symbols = tokens[4]

            start_bin = int(start / bin_size)
            end_bin = int(end / bin_size)

            for bin in range(start_bin, end_bin + 1):
                # apparently refseq genes from the ucsc always report
                # coordinates on the forward strand regardless of orientation
                self._tad_bins[chr][bin] = len(self._tads[chr])

                # if chr == 'chr1':
                #print(chr, location, bin, len(self._tads), file=sys.stderr)

            self._tads[chr].append([location, ids, gene_symbols])

            # if chr == 'chr1':
            #    print(self._tads[chr], file=sys.stderr)

        f.close()

        #print(self._tad_bins, file=sys.stderr)
        #print(self._tad_bins['chr10'], file=sys.stderr)

        # exit(0)

        print('Finished loading TADs.', file=sys.stderr)

    def get_names(self):
        return [TADS_HEADING, TAD_HEADING]

    def update_row(self, location:pychipseq.genomic.Location, row_map:Mapping[str, Any]):
        chr = location.chr

        if chr not in self._tad_bins:
            return pychipseq.text.NA

        start_bin = int(location.start / self._bin_size)
        end_bin = int(location.end / self._bin_size)

        # first find all the variants we might belong to
        genes = set()
        domains = []

        for bin in range(start_bin, end_bin + 1):
            if bin not in self._tad_bins[chr]:
                continue

            tad = self._tads[chr][self._tad_bins[chr][bin]]

            

            # if (chr == 'chr10'):
            #    print('lll', bin, location, tad[0], self._tad_bins[chr][bin], pychipseq.genomic.is_overlapping(location, tad[0]), file=sys.stderr)

            if pychipseq.genomic.is_overlapping(location, tad[0]):
                loc = str(tad[0])

                if loc not in domains:
                    domains.append(loc)

                #print('hug', location, tad[0], file=sys.stderr)
                genes.update(tad[2].split(';'))

        domains = ';'.join(domains)

        if domains == '':
            domains = pychipseq.text.NA

        genes = ';'.join(sorted(genes))

        if genes == '':
            genes = pychipseq.text.NA

        row_map[TADS_HEADING] = domains
        row_map[TAD_HEADING] = genes

        return [domains, genes]


class IsTADAnnotation(pychipseq.genomic.Annotation):
    """
    Annotate peaks for all possible refseq genes they might overlap.
    """

    def get_names(self):
        return [IS_TAD_HEADING]

    def update_row(self, location:pychipseq.genomic.Location, row_map:Mapping[str, Any]):
        
        #print(row_map, file=sys.stderr)
        #print(row_map[pychipseq.headings.GENE_SYMBOL], file=sys.stderr)

        if TAD_HEADING not in row_map:
            return '0'
        
        if pychipseq.headings.GENE_SYMBOL not in row_map:
            return '0'

        gene = row_map[pychipseq.headings.GENE_SYMBOL]
        ret = 1 if gene in row_map[TAD_HEADING] else 0

        return [ret]


class IsClosestTADAnnotation(pychipseq.genomic.Annotation):
    """
    Annotate peaks for all possible refseq genes they might overlap.
    """

    def get_names(self):
        return [IS_CLOSEST_TAD_HEADING]

    def update_row(self, location:pychipseq.genomic.Location, row_map:Mapping[str, Any]):
        if TAD_HEADING not in row_map:
            return [0]
        
        if pychipseq.headings.CLOSEST_GENE_SYMBOL not in row_map:
            return [0]

        gene = row_map[pychipseq.headings.CLOSEST_GENE_SYMBOL]
        ret = 1 if gene in row_map[TAD_HEADING] else 0

        return [ret]
