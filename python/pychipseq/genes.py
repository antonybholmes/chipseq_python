# -*- coding: utf-8 -*-
"""
Created on Mon Sep    8 15:12:34 2014

@author: Antony Holmes
"""

import sys
import collections
import re

import pychipseq.annotation
import pychipseq.expression

import pychipseq.genomic
import pychipseq.human.genomic
import pychipseq.headings
import pychipseq.text
import pychipseq.sample

GENE_ID = 'gene_id'
GENE_SYMBOL = 'gene_symbol'
GENE_REFSEQ = "refseq"
GENE_ENTREZ = 'entrez'
GENE_ENSEMBL = 'ensembl'


def parse_rdf_gene_id(text):
    """
    An RDF id contains a core gene id plus a decimal to indicate the variant
    """

    return re.match(r'(RDF\d+).*', text).group(1)


def parse_rdf_gene_variant_id(text):
    return int(re.match(r'.*(\d+)$', text).group(1))


def create_variant(id, chr, start, end):
    return f'{id}#{pychipseq.genomic.location_string(chr, start, end)}'


def parse_location_from_variant(location):
    matcher = re.match(r'.+(chr.+):(\d+)-(\d+)', location)

    chr = matcher.group(1)
    start = int(matcher.group(2))
    end = int(matcher.group(3))

    location = pychipseq.genomic.Location(chr, start, end)

    return location


def parse_id_from_variant(variant_id):
    matcher = re.match(r'^([^#]+).*', variant_id)

    return matcher.group(1)


def find_best_p_value(header, tokens):
    # Required for overlap files where both p-values from each
    # replicate are maintained. We need to select one.

    min_p = 1

    indices = pychipseq.text.find_indices(header, pychipseq.headings.P_VALUE)

    for i in indices:
        if tokens[i] != pychipseq.text.NA:
            p = float(tokens[i])

            if p < min_p:
                min_p = p

    return min_p


def find_best_score(header, tokens):
    # Required for overlap files where both scores from each
    # replicate are maintained. We need to select one.

    max_score = 0

    indices = pychipseq.text.find_indices(header, pychipseq.headings.SCORE)

    for i in indices:
        if tokens[i] != pychipseq.text.NA:
            s = float(tokens[i])

            if s > max_score:
                max_score = s

    return max_score


class Gene(pychipseq.genomic.Location):
    def __init__(self, id, symbol, strand, chr, start, end):
        super().__init__(chr, start, end)
        self._strand = strand
        self._id_map = {}
        self._id_map[GENE_ID] = id
        self._id_map[GENE_SYMBOL] = symbol

    def get_id(self, name: str) -> str:
        """
        Returns a named id associated with the gene or n/a if not present

        Args:
            name:   name of id, e.g. "refseq"
        Returns:
            Named argument.
        """
        return self._id_map.get(name, pychipseq.text.NA)

    @property
    def id(self):
        return self.get_id(GENE_ID)

    @property
    def symbol(self):
        return self.get_id(GENE_SYMBOL)

    @property
    def name(self):
        return self.symbol

    @property
    def strand(self):
        return self._strand


class RefSeqGene(Gene):
    def __init__(self, refseq, entrez, symbol, strand, chr, start, end):
        super().__init__(refseq, symbol, strand, chr, start, end)
        self._id_map[GENE_REFSEQ] = refseq
        self._id_map[GENE_ENTREZ] = entrez

    @property
    def refseq(self):
        return self.get_id(GENE_REFSEQ)

    @property
    def entrez(self):
        return self.get_id(GENE_ENTREZ)


class RefSeqGenes:
    """
    Gene lookup by variant id, other dbs store variant ids
    # instead of objects. This is to prevent redun
    """

    def __init__(self, file):
        self.genes = collections.defaultdict(Gene)

        print(f'Loading genes from {file}...', file=sys.stderr)

        # pychipseq.annotation.REFSEQ_FILE
        f = open(file, 'r')

        # skip header
        f.readline()

        # To account for multiple versions of a gene, allocate each entrez
        # id a unique index

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            tokens = line.split('\t')

            refseq = tokens[0]
            entrez = tokens[1]
            symbol = tokens[2]
            chr = tokens[3]
            strand = tokens[4]
            # ucsc convention
            start = int(tokens[5]) + 1
            end = int(tokens[6])

            if re.match(r'.*n/a.*', entrez):
                continue

            if re.match(r'.*n/a.*', symbol):
                continue

            if re.match(r'.*MIR.*', symbol):
                continue

            # index gene on an id and coordinates to keep it unique

            variant_id = create_variant(refseq, chr, start, end)

            self.genes[variant_id] = RefSeqGene(
                refseq, entrez, symbol, strand, chr, start, end)

        f.close()

    def get_gene(self, variant_id: str):
        return self.genes[variant_id]