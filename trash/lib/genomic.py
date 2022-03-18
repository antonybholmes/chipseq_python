# -*- coding: utf-8 -*-
"""
Functions for genomics

Created on Sat Jan 31 16:49:34 2015

@author: antony
"""

import collections
import re

import lib.search
import lib.text

TELOMERE_SIZE = 100000


class Location(object):
    def __init__(self, chr, start, end):
        self.chr = chr
        self.start = start
        self.end = end
        self.length = end - start + 1

    def to_string(self):
        return self.chr + ":" + str(self.start) + "-" + str(self.end)


def is_location(location):
    '''
    Returns true if location is a genomic location
    '''
    return re.match(r'(chr.+):(\d+)-(\d+)', location)


def is_chr(location):
    '''
    Returns true if location is a chr
    '''
    return re.match(r'(chr.+)', location)


def parse_location(location, padding5p=0, padding3p=0):
    matcher = re.match(r'.*(chr.+):(\d+)-(\d+).*', location)

    chr = matcher.group(1)
    start = int(matcher.group(2))
    end = int(matcher.group(3))

    return Location(chr, start - padding5p, end + padding3p)


def parse_location_cols(tokens, offset=0):
    return Location(tokens[offset], int(tokens[offset + 1]), int(tokens[offset + 2]))


def mid_point(location):
    return int((location.start + location.end) / 2)


def center_location(location):
    c = mid_point(location)

    return Location(location.chr, c, c)


def pad_location(location, padding5p=0, padding3p=0):
    return Location(location.chr, location.start - padding5p, location.end + padding3p)


def overlap_locations(location1, location2):
    if location1.chr != location2.chr:
        return None

    if location1.end < location2.start or \
            location2.end < location1.start or \
            location1.start > location2.end or \
            location2.start > location1.end:
        return None

    start = max(location1.start, location2.start)
    end = min(location1.end, location2.end)

    return Location(location1.chr, start, end)


def is_overlapping(location1, location2):
    if location1.chr != location2.chr:
        return False

    max_start = max(location1.start, location2.start)
    max_end = max(location1.end, location2.end)

    min_start = min(location1.start, location2.start)
    min_end = min(location1.end, location2.end)

    return max_start >= min_start and max_start < min_end


def mid(location):
    return (location.start + location.end) / 2


def get_closest_tss(tss):
    closest = "n/a"

    min_d = 1000000
    min_abs_d = 1000000

    for p in tss:
        if p == "n/a":
            continue

        d = int(p)
        abs_d = abs(d)

        if abs_d < min_abs_d:
            min_d = d
            min_abs_d = abs_d

    if min_d != 1000000:
        closest = str(min_d)

    return closest


def load_locations(file, center=False, padding5p=0, padding3p=0):
    """
    Load a set of locations
    """

    f = open(file, 'r')

    header = f.readline().strip()

    ret = []

    for line in f:
        l = parse_location(line.strip())

        if center:
            l = center_location(l)

        l = pad_location(l, padding5p, padding3p)

        ret.append(l)

    f.close()

    return ret, header


class SearchGenomicFeatures(lib.search.GappedSearch):
    """
    Instance of gapped search for quickly identifying if a region
    is close to centromere
    """

    def __init__(self, file):
        lib.search.GappedSearch.__init__(self)

        f = open(file, 'r')

        f.readline()

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            location = lib.genomic.parse_location(line)

            # we store locations for determining if they overlap or not
            self.add_feature(location, location)

        f.close()

        self.lock()


class SearchBed(lib.search.BlockSearch):  # (lib.search.GappedSearch):
    """
    Instance of gapped search for bed files allowing for named regions
    """

    def __init__(self, file):
        lib.search.BlockSearch.__init__(self)  # GappedSearch

        f = open(file, 'r')

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            if "track" in line:
                continue

            tokens = line.split("\t")

            location = Location(tokens[0], int(tokens[1]), int(tokens[2]))

            if len(tokens) > 3:
                name = tokens[3]
            else:
                name = ""

            bed = (location, name)

            # we store locations for determining if they overlap or not
            self.add_feature(location, bed)

        f.close()


# (lib.search.GappedSearch):
class SearchGenomicBedFeatures(lib.search.BlockSearch):
    """
    Instance of gapped search for bed files
    """

    def __init__(self, file):
        lib.search.BlockSearch.__init__(self)  # GappedSearch

        f = open(file, 'r')

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            if "track" in line:
                continue

            tokens = line.split("\t")

            location = Location(tokens[0], int(tokens[1]), int(tokens[2]))

            # we store locations for determining if they overlap or not
            self.add_feature(location, location)

        f.close()


class GenomicFeaturesOverlap:
    """
    Uses a gapped search to determine by how much a location overlaps a feature
    """

    def __init__(self, gapped_search):
        self.gapped_search = gapped_search

    def get_features(self, location):
        return self.gapped_search.get_features(location)

    def get_max_overlap(self, location):
        features = self.get_features(location)

        max_overlap_width = -1
        max_overlap = None

        for feature in features:
            for l in feature.values:
                overlap = overlap_locations(location, l)

                if overlap is not None:
                    if overlap.length > max_overlap_width:
                        max_overlap_width = overlap.length
                        max_overlap = overlap

        return max_overlap


class GenomicBedOverlap(GenomicFeaturesOverlap):
    """
    Uses a gapped search to determine by how much a location overlaps a feature
    """

    def __init__(self, gapped_search):
        super().__init__(gapped_search)

    def get_regions(self, location):
        features = self.get_features(location)

        ret = []

        for feature in features:
            for l in feature.values:
                ret.append(l)

        return ret

    def get_overlapping_regions(self, location):
        regions = self.get_regions(location)

        ret = []

        for r in regions:
            overlap = overlap_locations(location, r[0])
            if overlap is not None:
                ret.append(r)

        return ret


class Chromosomes(object):
    """
    Chromosome sizes
    """

    def __init__(self):
        self.sizes = collections.defaultdict(int)

        f = open(
            "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/hg19_chromosome_sizes.txt", "r")

        f.readline()

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            tokens = line.split("\t")

            chr = tokens[0]
            size = int(tokens[1])

            self.sizes[chr] = size

        f.close()

    def get_size(self, chr):
        return self.sizes[chr]


class Telomeres(object):
    """
    Determine whether a location overlaps a centromere
    """

    def __init__(self):
        self.chromosomes = Chromosomes()

    def get_classification(self, location):
        end = self.chromosomes.get_size(location.chr) - TELOMERE_SIZE

        if location.start <= TELOMERE_SIZE:
            classification = "peri_telomeric"
        elif location.end >= end:
            classification = "peri_telomeric"
        else:
            classification = lib.text.NA

        return classification


class Centromeres(object):
    """
    Determine whether a location overlaps a centromere
    """

    def __init__(self):
        centromeres = SearchGenomicBedFeatures(
            "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/ucsc_centromeres_hg19.bed")
        pericentromeres = SearchGenomicBedFeatures(
            "/ifs/scratch/cancer/Lab_RDF/ngs/references/rdf/rdf_pericentromeres_hg19.bed")

        self.cen_overlaps = GenomicFeaturesOverlap(centromeres)
        self.p_cen_overlaps = GenomicFeaturesOverlap(pericentromeres)

    def get_classification(self, location):
        classification = lib.text.NA

        #
        # Are we in a centromere
        #

        in_centromere = False

        overlap = self.cen_overlaps.get_max_overlap(location)

        p = 0

        if overlap is not None:
            p = overlap.length / location.length

            classification = "centromeric"
            in_centromere = True

        #
        # Are we in a pericentromere
        #

        if not in_centromere:
            overlap = self.p_cen_overlaps.get_max_overlap(location)

            if overlap is not None:
                p = overlap.length / location.length
                classification = "peri_centromeric"

        return classification


class Annotation(object):
    def get_name(self):
        raise NotImplementedError()

    def annotate(self, location):
        """
        Should return a list of items
        """
        raise NotImplementedError()


class Repetitive(Annotation):
    """
    Determine whether a location overlaps a repetitive region
    """

    def __init__(self):
        self.centromeres = Centromeres()
        self.telomeres = Telomeres()

    def get_classification(self, location):
        classifications = set()

        classifications.add(self.centromeres.get_classification(location))
        classifications.add(self.telomeres.get_classification(location))

        # If there are multiple classifications, get rid of the n/a
        if len(classifications) > 1:
            classifications.remove(lib.text.NA)

        return sorted(classifications)

    def get_name(self):
        return "Centromere/Telomere"

    def annotate(self, location):
        return ",".join(sorted(self.get_classification(location)))
