from .. import tad

TAD_FILE = '/ifs/scratch/cancer/Lab_RDF/ngs/references/rdf/tad/hg19/tads.gencode.v38lift37.genes.approved.tsv'


class GencodeTADAnnotation(tad.TADAnnotation):
    def __init__(self, bin_size=100):
        super().__init__(TAD_FILE, bin_size=bin_size)
