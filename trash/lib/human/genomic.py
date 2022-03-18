import lib.text
import lib.genomic

class SimpleTandemRepeats(lib.genomic.Annotation):
    """
    Determine whether a location overlaps a centromere
    """
    
    def __init__(self):
        trf = lib.genomic.SearchGenomicBedFeatures("/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/assembly/hg19/simple_tandem_repeats_hg19.bed")

        self.trf_overlaps = lib.genomic.GenomicFeaturesOverlap(trf)

    def get_name(self):
        return "Simple Tandem Repeats"
    
    def annotate(self, location):
        overlap = self.trf_overlaps.get_max_overlap(location)
        
        if overlap is not None:
            classification = "tandem_repeat"
        else:
            classification = lib.text.NA
        
        return classification


class EncodeBlacklist(lib.genomic.Annotation):
    """
    Determine whether a location overlaps a centromere
    """
    
    def __init__(self):
        f = lib.genomic.SearchGenomicBedFeatures('/ifs/scratch/cancer/Lab_RDF/ngs/references/encode/chipseq/blacklist.bed')

        self.trf_overlaps = lib.genomic.GenomicFeaturesOverlap(f)

    def get_name(self):
        return "ENCODE blacklist"
    
    def annotate(self, location):
        overlap = self.trf_overlaps.get_max_overlap(location)
        
        if overlap is not None:
            classification = "encode_blacklist"
        else:
            classification = lib.text.NA
        
        return classification
    
class GiuliaBlacklist(lib.genomic.Annotation):
    """
    Determine whether a location overlaps a centromere
    """
    
    def __init__(self):
        f = lib.genomic.SearchGenomicBedFeatures('/ifs/scratch/cancer/Lab_RDF/ngs/references/rdf/rdf_giulia_blacklist_hg19.bed')

        self.trf_overlaps = lib.genomic.GenomicFeaturesOverlap(f)

    def get_name(self):
        return "Giulia blacklist"
    
    def annotate(self, location):
        overlap = self.trf_overlaps.get_max_overlap(location)
        
        if overlap is not None:
            classification = "giulia_blacklist"
        else:
            classification = lib.text.NA
        
        return classification


class Nnnn(lib.genomic.Annotation):
    """
    Determine whether a location overlaps a centromere
    """
    
    def __init__(self):
        nnnn = lib.genomic.SearchGenomicBedFeatures("/ifs/scratch/cancer/Lab_RDF/ngs/references/ucsc/assembly/hg19/nnnn_hg19.bed")

        self.nnnn_overlaps = lib.genomic.GenomicFeaturesOverlap(nnnn)

    def get_name(self):
        return "NNNNs"
    
    def annotate(self, location):
        overlap = self.nnnn_overlaps.get_max_overlap(location)
        
        if overlap is not None:
            classification = "nnnn"
        else:
            classification = lib.text.NA
        
        return classification
