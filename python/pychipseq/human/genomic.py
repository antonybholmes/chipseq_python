import collections
from typing import Any, Mapping, Union
import pychipseq.text
import pychipseq.genomic


class Chromosomes:
    """
    Chromosome sizes
    """

    def __init__(self):
        self._sizes = collections.defaultdict(int)

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

            self._sizes[chr] = size

        f.close()

    def get_size(self, chr):
        return self._sizes[chr]


class Telomeres(pychipseq.genomic.ClassifyRegion):
    """
    Determine whether a location overlaps a centromere
    """

    def __init__(self):
        self._chromosomes = Chromosomes()

    def get_classification(self, location: pychipseq.genomic.Location):
        end = self._chromosomes.get_size(
            location.chr) - pychipseq.genomic.TELOMERE_SIZE

        if location.start <= pychipseq.genomic.TELOMERE_SIZE:
            classification = "peri_telomeric"
        elif location.end >= end:
            classification = "peri_telomeric"
        else:
            classification = pychipseq.text.NA

        return classification


class Centromeres(pychipseq.genomic.ClassifyRegion):
    """
    Determine whether a location overlaps a centromere
    """

    def __init__(self):
        centromeres = pychipseq.genomic.SearchGenomicBedFeatures(
            "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/ucsc_centromeres_hg19.bed")
        pericentromeres = pychipseq.genomic.SearchGenomicBedFeatures(
            "/ifs/scratch/cancer/Lab_RDF/ngs/references/rdf/rdf_pericentromeres_hg19.bed")

        self.cen_overlaps = pychipseq.genomic.GenomicFeaturesOverlap(
            centromeres)
        self.p_cen_overlaps = pychipseq.genomic.GenomicFeaturesOverlap(
            pericentromeres)

    def get_classification(self, location: pychipseq.genomic.Location):
        classification = pychipseq.text.NA

        # Are we in a centromere
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


class Repetitive(pychipseq.genomic.Annotation):
    """
    Determine whether a location overlaps a repetitive region
    """

    def __init__(self):
        self._centromeres = Centromeres()
        self._telomeres = Telomeres()

    def get_classification(self, location: pychipseq.genomic.Location):
        classifications = set()

        classifications.add(self._centromeres.get_classification(location))
        classifications.add(self._telomeres.get_classification(location))

        # If there are multiple classifications, get rid of the n/a
        if len(classifications) > 1:
            classifications.remove(pychipseq.text.NA)

        return sorted(classifications)

    def get_names(self):
        return ["Centromere/Telomere"]

    def update_row(self, location: pychipseq.genomic.Location, row_map: Mapping[str, Union[str, int, float]]):
        ret = ','.join(sorted(self.get_classification(location)))

        return [ret]


class SimpleTandemRepeats(pychipseq.genomic.Annotation):
    """
    Determine whether a location overlaps a centromere
    """

    def __init__(self):
        trf = pychipseq.genomic.SearchGenomicBedFeatures(
            "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/assembly/hg19/simple_tandem_repeats_hg19.bed")
        self._trf_overlaps = pychipseq.genomic.GenomicFeaturesOverlap(trf)

    def get_names(self):
        return ["Simple Tandem Repeats"]

    def update_row(self, location: pychipseq.genomic.Location, row_map: Mapping[str, Union[str, int, float]]):
        overlap = self._trf_overlaps.get_max_overlap(location)

        if overlap is not None:
            classification = "tandem_repeat"
        else:
            classification = pychipseq.text.NA

        return [classification]


class EncodeBlacklist(pychipseq.genomic.Annotation):
    """
    Determine whether a location overlaps a centromere
    """

    def __init__(self):
        f = pychipseq.genomic.SearchGenomicBedFeatures(
            '/ifs/scratch/cancer/Lab_RDF/ngs/references/encode/chipseq/blacklist.bed')

        self._trf_overlaps = pychipseq.genomic.GenomicFeaturesOverlap(f)

    def get_names(self):
        return ["ENCODE blacklist"]

    def update_row(self, location: pychipseq.genomic.Location, row_map: Mapping[str, Union[str, int, float]]):
        overlap = self._trf_overlaps.get_max_overlap(location)

        if overlap is not None:
            classification = "encode_blacklist"
        else:
            classification = pychipseq.text.NA

        return [classification]


class GiuliaBlacklist(pychipseq.genomic.Annotation):
    """
    Determine whether a location overlaps a centromere
    """

    def __init__(self):
        f = pychipseq.genomic.SearchGenomicBedFeatures(
            '/ifs/scratch/cancer/Lab_RDF/ngs/references/rdf/rdf_giulia_blacklist_hg19.bed')

        self._trf_overlaps = pychipseq.genomic.GenomicFeaturesOverlap(f)

    def get_names(self):
        return ["Giulia blacklist"]

    def update_row(self, location: pychipseq.genomic.Location, row_map: Mapping[str, Any]):
        overlap = self._trf_overlaps.get_max_overlap(location)

        if overlap is not None:
            classification = "giulia_blacklist"
        else:
            classification = pychipseq.text.NA

        return [classification]


class Nnnn(pychipseq.genomic.Annotation):
    """
    Determine whether a location overlaps a centromere
    """

    def __init__(self):
        nnnn = pychipseq.genomic.SearchGenomicBedFeatures(
            "/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/assembly/hg19/nnnn_hg19.bed")

        self._nnnn_overlaps = pychipseq.genomic.GenomicFeaturesOverlap(nnnn)

    def get_names(self):
        return ["NNNNs"]

    def update_row(self, location: pychipseq.genomic.Location, row_map: Mapping[str, Any]):
        overlap = self._nnnn_overlaps.get_max_overlap(location)

        if overlap is not None:
            classification = "nnnn"
        else:
            classification = pychipseq.text.NA

        return [classification]
