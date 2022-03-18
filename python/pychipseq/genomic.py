# -*- coding: utf-8 -*-
"""
Functions for genomics

Created on Sat Jan 31 16:49:34 2015

@author: antony
"""

import collections
import re
import sys
from typing import Any, List, Mapping, Optional, Union

import pychipseq.text

TELOMERE_SIZE = 100000


class Location:
    def __init__(self, chr, start, end):
        self._chr = chr
        self._start = start
        self._end = end
        self._length = end - start + 1

    @property
    def chr(self) -> str:
        return self._chr

    @property
    def start(self) -> int:
        return self._start

    @property
    def end(self) -> int:
        return self._end

    @property
    def length(self) -> int:
        return self._length

    def to_string(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        return location_string(self.chr, self.start, self.end)

    def __repr__(self) -> str:
        return self.__str__()


class FeatureSet:
    def __init__(self, start):
        self._start: int = start
        self._values: List[Any] = []

    @property
    def start(self):
        return self._start

    @property
    def values(self):
        return self._values

    def add(self, value: Any):
        self._values.append(value)

    def __iter__(self):
        return iter(self._values)


class GappedSearch:
    def __init__(self):
        self._features = collections.defaultdict(
            lambda: collections.defaultdict(FeatureSet))
        self._starts = collections.defaultdict(list)
        self._locked: bool = False

    def add_feature(self, location: Location, feature: Any):
        self._locked = False

        if location.chr not in self.features or location.start not in self.features[location.chr]:
            self.features[location.chr][location.start] = FeatureSet(
                location.start)

        if location.chr not in self.features or location.end not in self.features[location.chr]:
            self.features[location.chr][location.end] = FeatureSet(
                location.end)

        self._features[location.chr][location.start].add(feature)
        self._features[location.chr][location.end].add(feature)

    def _sort_starts(self):
        """
        Sort the starts for binary searching
        """

        if self._locked:
            return

        # We need the locations sorted in order for a binary search
        for chr in self._features:
            self._starts[chr] = sorted(self.features[chr])

        self._locked = True

    def get_closest_features(self, location: Location) -> FeatureSet:
        featureSets: List[FeatureSet] = self.get_features(location)

        min_d = sys.maxint
        ret = None

        for featureSet in featureSets:
            d = abs(mid(location) - featureSet.start)

            if d < min_d:
                ret = featureSet
                min_d = d

        return ret

    def get_features(self, location: Location) -> List[FeatureSet]:
        """
        Return a list of features closest to a location which then
        be tested for overlap etc
        """

        ret = []

        if location.chr not in self._features:
            return ret

        self._sort_starts()

        starts = self._starts[location.chr]

        s = self._get_start_index(starts, location.start, True)
        e = self._get_start_index(starts, location.end, False)

        chr_features = self._features[location.chr]

        for i in range(s, e + 1):
            ret.append(chr_features[starts[i]])

        return ret

    def _get_start_index(self, indices: List[int], start: int, startMode: bool):
        if len(indices) < 2:
            return 0

        s = 0

        if start <= indices[s]:
            return s

        e = len(indices) - 1

        if start >= indices[e]:
            return e

        while e - s > 1:
            im = int((s + e) / 2) #s + (e - s) // 2

            pm = indices[im]

            if pm > start:
                e = im # - 1
            elif pm < start:
                s = im # + 1
            else:
                return im

        # If we made it this far, we have narrowed down the search to
        # two position so return based on trying to cover the region
        # of interest. We can narrow down the actual closest bin later.

        if startMode:
            # return the closest before the start
            return s
        else:
            # return the closest after the end
            return e


class BlockSearch:
    def __init__(self, block_size: int = 10000):
        self._block_size = block_size
        self._features = collections.defaultdict(
            lambda: collections.defaultdict(FeatureSet))

    def add_feature(self, location: Location, feature: Any):
        s = int(location.start / self._block_size)
        e = int(location.end / self._block_size)

        for b in range(s, e + 1):
            if location.chr not in self._features or b not in self._features[location.chr]:
                self._features[location.chr][b] = FeatureSet(location.start)

            self._features[location.chr][b].add(feature)

    def get_closest_features(self, location: Location) -> FeatureSet:
        featureSets = self.get_features(location)

        min_d = sys.maxint
        ret = None

        for featureSet in featureSets:
            d = abs(mid(location) - featureSet.start)

            if d < min_d:
                ret = featureSet
                min_d = d

        return ret

    def get_features(self, location: Location) -> List[FeatureSet]:
        """
        Return a list of features closest to a location which then
        be tested for overlap etc
        """

        ret = []

        if location.chr not in self._features:
            return ret

        s = int(location.start / self._block_size)
        e = int(location.end / self._block_size)

        chr_features = self._features[location.chr]

        for b in range(s, e + 1):
            if b in chr_features:
                ret.append(chr_features[b])

        return ret


def location_string(chr: str, start: Union[int, str], end: Union[int, str]) -> str:
    """
    Returns a standardized string representation of a genomic location.

    Args:
        chr:    chromosome
        start:  1-based start location
        end:    1-based end location
    Returns:
        A chr:start-end formatted string.
    """
    return f'{chr}:{str(start)}-{str(end)}'


def is_location(location: str):
    """
    Returns true if location is a genomic location
    """
    return re.match(r'(chr.+):(\d+)-(\d+)', location)


def is_chr(location: str):
    """
    Returns true if location is a chr
    """
    return re.match(r'(chr.+)', location)


def parse_location(location: str, padding5p: int = 0, padding3p: int = 0) -> Location:
    """
    Parses a location string, e.g. 'chr1:1-100' into a location object.

    Args:
        location:       A location string.
        padding5p:      Add padding to the 5p end.

    Returns:
        A location object representing the location string.
    """

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


def overlap_locations(location1: Location, location2: Location) -> Union[None, Location]:
    if location1.chr != location2.chr:
        return None

    max_start = max(location1.start, location2.start)
    #max_end = max(location1.end, location2.end)
    #min_start = min(location1.start, location2.start)
    min_end = min(location1.end, location2.end)

    if min_end < max_start:
        return None

    # if location1.end < location2.start or \
    #         location2.end < location1.start or \
    #         location1.start > location2.end or \
    #         location2.start > location1.end:
    #     return None

    #start = max(location1.start, location2.start)
    #end = min(location1.end, location2.end)

    return Location(location1.chr, max_start, min_end)


def is_overlapping(location1: Location, location2: Location):
    if location1.chr != location2.chr:
        return False

    max_start = max(location1.start, location2.start)
    #max_end = max(location1.end, location2.end)
    #min_start = min(location1.start, location2.start)
    min_end = min(location1.end, location2.end)

    return max_start <= min_end  # return max_start >= min_start and max_start < min_end


def mid(location: Location):
    return (location.start + location.end) / 2


def get_closest_tss(tss):
    closest = pychipseq.text.NA

    min_d = 1000000
    min_abs_d = 1000000

    for p in tss:
        if p == "n/a" or p == "NA":
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


class SearchGenomicFeatures(GappedSearch):
    """
    Instance of gapped search for quickly identifying if a region
    is close to centromere
    """

    def __init__(self, file: str):
        super().__init__()

        f = open(file, 'r')

        f.readline()

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            location: Location = parse_location(line)

            # we store locations for determining if they overlap or not
            self.add_feature(location, location)

        f.close()

        self.lock()


class SearchBed(BlockSearch):  # (lib.search.GappedSearch):
    """
    Instance of gapped search for bed files allowing for named regions
    """

    def __init__(self, file: str):
        super().__init__()  # GappedSearch

        f = open(file, 'r')

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            if 'track' in line:
                continue

            tokens = line.split('\t')

            location = Location(tokens[0], int(tokens[1]), int(tokens[2]))

            if len(tokens) > 3:
                name = tokens[3]
            else:
                name = ''

            bed = (location, name)

            # we store locations for determining if they overlap or not
            self.add_feature(location, bed)

        f.close()


# (lib.search.GappedSearch):
class SearchGenomicBedFeatures(BlockSearch):
    """
    Instance of gapped search for bed files
    """

    def __init__(self, file: str):
        super().__init__()  # GappedSearch

        f = open(file, 'r')

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            if 'track' in line:
                continue

            tokens = line.split('\t')

            location = Location(tokens[0], int(tokens[1]), int(tokens[2]))

            # we store locations for determining if they overlap or not
            self.add_feature(location, location)

        f.close()


class GenomicFeaturesOverlap:
    """
    Uses a gapped search to determine by how much a location overlaps a feature
    """

    def __init__(self, gapped_search: GappedSearch):
        self._gapped_search = gapped_search

    def get_features(self, location: Location):
        return self._gapped_search.get_features(location)

    def get_max_overlap(self, location: Location):
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

    def __init__(self, gapped_search: GappedSearch):
        super().__init__(gapped_search)

    def get_regions(self, location: Location):
        features = self.get_features(location)

        ret = []

        for feature in features:
            for l in feature.values:
                ret.append(l)

        return ret

    def get_overlapping_regions(self, location: Location):
        regions = self.get_regions(location)

        ret = []

        for r in regions:
            overlap = overlap_locations(location, r[0])
            if overlap is not None:
                ret.append(r)

        return ret


class Annotation:
    def get_names(self) -> List[str]:
        """
        Should return the name of the column or an array of columns
        """
        raise NotImplementedError()

    def annotate(self, location: Location) -> List[Union[str, int, float]]:
        """
        Should return a list of items
        """
        raise NotImplementedError()

    def update_row(self, location: Location, row_map: Mapping[str, Any]) -> List[Union[str, int, float]]:
        return self.annotate(location)


class ClassifyRegion:
    """
    Returns a class label for a region, such as identifying
    a telomere.
    """

    def get_classification(self, location: Location) -> str:
        raise NotImplementedError()
