# -*- coding: utf-8 -*-
"""
Created on Mon Sep    8 15:12:34 2014

@author: Antony Holmes
"""

import sys
import collections
import re

import lib.annotation
import lib.expression
import lib.tss

import lib.genomic
import lib.human.genomic
import lib.headings
import lib.text
import lib.sample


class Gene(object):
    def __init__(self, refseq, entrez, symbol, strand, chr, start, end):
        self.refseq = refseq        
        self.entrez = entrez
        self.symbol = symbol
        self.strand = strand
        self.chr = chr
        self.start = start
        self.end = end
        
        
class RefSeqGenes(object):
    """
    Gene lookup
    """
    
    def __init__(self, file):
        self.genes = collections.defaultdict(Gene)
     
    
        sys.stderr.write("Loading genes from " + file + "...\n")
    
        # lib.annotation.REFSEQ_FILE
        f = open(file, 'r')

        # skip header
        f.readline()
    
        # To account for multiple versions of a gene, allocate each entrez
        # id a unique index
    
        for line in f:
            line = line.strip()
            
            if len(line) == 0:
                continue
            
            tokens = line.split("\t")
            
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
            
            # If the gene is on the negative strand, the promoter is around the
            # end coordinate on the forward strand
            
            variant_id = lib.sample.create_variant(refseq, chr, start, end)
 
            self.genes[variant_id] = Gene(refseq, entrez, symbol, strand, chr, start, end)
            
        f.close()
        
        
    def get_gene(self, variant_id):
        return self.genes[variant_id]
        
        
def find_best_p_value(header, tokens):
    # Required for overlap files where both p-values from each
    # replicate are maintained. We need to select one.
    
    min_p = 1
    
    indices = lib.text.find_indices(header, lib.headings.P_VALUE)
    
    for i in indices:
        if tokens[i] != lib.text.NA:
            p = float(tokens[i])
            
            if p < min_p:
                min_p = p
                
    return min_p


def find_best_score(header, tokens):
    # Required for overlap files where both scores from each
    # replicate are maintained. We need to select one.
    
    max_score = 0
    
    indices = lib.text.find_indices(header, lib.headings.SCORE)
    
    for i in indices:
        if tokens[i] != lib.text.NA:
            s = float(tokens[i])
            
            if s > max_score:
                max_score = s
                
    return max_score
    

class AnnotatePeak(object):
    """
    Core annotation for annotating peaks/regions
    """
    
    def __init__(self, \
        type, \
        refseq_annotation, \
        refseq_genes, \
        refseq_tss, \
        refseq_end, \
        prom_ext_5p=5000, \
        prom_ext_3p=4000, \
        bin_size=10000):
        
        self.type = type
        #self.prom_ext_5p = prom_ext_5p
        #self.prom_ext_3p = prom_ext_3p
        #self.bin_size = 10000

        self.mir_max_gap = prom_ext_5p
        self.annotation = refseq_annotation #lib.human.tss.RefSeqAnnotation(prom_ext_5p, prom_ext_3p, self.bin_size)
        self.genes = refseq_genes #RefSeqGenes()

        self.tss_refseq = refseq_tss #lib.human.tss.RefSeqTss(prom_ext_5p, prom_ext_3p)
        self.refseq_end = refseq_end #lib.human.tss.RefSeqEnd(prom_ext_5p, prom_ext_3p)

        # For annotating if in centromere or not
        #self.repetitive = lib.genomic.Repetitive()

        #
        # Create a mir database
        #
        self.mir_annotation = lib.annotation.MirAnnotation(bin_size)
        self.mir_transcript_annotation = lib.annotation.MirTranscriptAnnotation()
        
        self.promoter_type = "(prom=-" + str(prom_ext_5p / 1000) + "/+" + str(prom_ext_3p / 1000) + "kb)"
        
        #
        # Modules are self contained annotation blocks to make it easier
        # to update them
        #
        
        self.annotation_modules = []
        
        self.annotation_modules.append(lib.genomic.Repetitive())
        #self.annotation_modules.append(lib.human.genomic.SimpleTandemRepeats())
        #self.annotation_modules.append(lib.human.genomic.Nnnn())


    def print_header(self):
        sys.stdout.write(lib.headings.REFSEQ_ID)
        sys.stdout.write("\t" + lib.headings.ENTREZ_ID)
        sys.stdout.write("\t" + lib.headings.GENE_SYMBOL)
        sys.stdout.write("\t" + self.type + " Relative To Gene " + self.promoter_type)
        sys.stdout.write("\t" + self.type + " " + lib.headings.TSS_DISTANCE)
        
        #
        # Closest
        #
        sys.stdout.write("\tClosest " + lib.headings.REFSEQ_ID)
        sys.stdout.write("\tClosest " + lib.headings.ENTREZ_ID) 
        sys.stdout.write("\tClosest " + lib.headings.GENE_SYMBOL)
        sys.stdout.write("\t" + self.type + " Relative To Closest Gene " + self.promoter_type)
        sys.stdout.write("\t" + self.type + " TSS Closest Distance")
        
        #
        # Closest End
        #
        
        sys.stdout.write("\tClosest End " + lib.headings.REFSEQ_ID)
        sys.stdout.write("\tClosest End " + lib.headings.ENTREZ_ID) 
        sys.stdout.write("\tClosest End " + lib.headings.GENE_SYMBOL)
        sys.stdout.write("\t" + self.type + " Relative To Closest Gene End " + self.promoter_type)
        sys.stdout.write("\t" + self.type + " End Closest Distance")
        
        #
        # Other
        #
        
        #sys.stdout.write("\t" + self.type + " Centromere Telomere")
        sys.stdout.write("\tmiR Symbol")
        sys.stdout.write("\t" + self.type + " Relative To miR")
        sys.stdout.write("\t" + self.type + " miR Start Closest Distance")
        sys.stdout.write("\t" + self.type + " miR Start Distance")
        
        #
        # Modules
        #
        
        for module in self.annotation_modules:
            sys.stdout.write("\t" + module.get_name())
        
        # End the header
        #sys.stdout.write("\n")


    def annotate(self, location):
        #sys.stderr.write(location + "\n")
        
        ret = []
        
        genes = self.annotation.annotate_location(location)

        refseqs = []
        symbols = []
        entrezes = []
        chrs = []
        strands = []
        starts = []
        ends = []
        types = []
        tss_list = []
    
        for gene in genes:
            gene_annotation = self.genes.get_gene(gene.id)
            
            chrs.append(gene_annotation.chr)
            starts.append(gene_annotation.start)
            strands.append(gene_annotation.strand)
            ends.append(gene_annotation.end)
            entrezes.append(gene_annotation.entrez)
            refseqs.append(gene_annotation.refseq)
            symbols.append(gene_annotation.symbol)

            # include tss distance for all peaks regardless
            tss_list.append(str(gene.d))

            types.append(",".join(reversed(sorted(gene.types))))
        
        
        if len(entrezes) > 0:
            ret.append(";".join(refseqs))
            ret.append(";".join(entrezes))
            ret.append(";".join(symbols))
            ret.append(";".join(types))
            ret.append(";".join(tss_list))
        else:
            ret.append("n/a")
            ret.append("n/a")
            ret.append("n/a")
            ret.append("n/a")
            ret.append("n/a")
        
        #
        # Look at which gene TSS is closest
        #
        
        #sys.stderr.write(location.to_string() + "\n")
        
        tss_gene = self.tss_refseq.get_closest_gene(location)
        
        gene_annotation = self.genes.get_gene(tss_gene.id)

        # Closest TSS gene symbol
        ret.append(gene_annotation.refseq)
        
        # Closest TSS gene entrez
        ret.append(gene_annotation.entrez)
    
        # Closest TSS gene symbol
        ret.append(gene_annotation.symbol)
        
        # Closest TSS gene types
        ret.append(",".join(reversed(sorted(tss_gene.types))))
    
        # Closest TSS gene distances, since there can only be one closest distance
        # all entries will have the same d so take the first. This is for
        # odd entries where there might be two or more genes with the same starts
    
        ret.append(str(tss_gene.d))
        
        
        #
        # Look at which gene end is closest
        #
        
        tss_gene = self.refseq_end.get_closest_gene(location)
        
        gene_annotation = self.genes.get_gene(tss_gene.id)

        # Closest TSS gene symbol
        ret.append(gene_annotation.refseq)
        
        # Closest TSS gene entrez
        ret.append(gene_annotation.entrez)
    
        # Closest TSS gene symbol
        ret.append(gene_annotation.symbol)
        
        # Closest TSS gene types
        ret.append(",".join(reversed(sorted(tss_gene.types))))
    
        # Closest TSS gene distances, since there can only be one closest distance
        # all entries will have the same d so take the first. This is for
        # odd entries where there might be two or more genes with the same starts
    
        ret.append(str(tss_gene.d))
        
        #
        # miRs
        #
 
        mir_annotations = collections.defaultdict(str)
        mss = collections.defaultdict(str)
        self.mir_annotation.annotate_location(location, self.mir_max_gap, mir_annotations, mss)
    
        mirs = []
        mir_types = []
        mss_list = []
    
        #
        # If the peak is a promoter, check to see if there is a mir in the gene
        #
    
        for gene in genes:
            gene_annotation = self.genes.get_gene(gene.id)
            
            promoter = False
        
            for type in gene.types:
                if re.match(r'.*promoter.*', type):
                    promoter = True
                    break
        
            if promoter == False:
                continue
        
            transcript_mirs = self.mir_transcript_annotation.find_mirs(gene_annotation.entrez)
        
            if len(transcript_mirs) == 0:
                continue

            symbol = gene_annotation.symbol
        
            if symbol == "n/a":
                symbol = "gene"
        
            for mir in sorted(transcript_mirs):
                mirs.append(mir);
                mir_types.append(symbol + "_promoter")
                mss_list.append("n/a")
    
        # Standard mir annotation
    
        for mir in sorted(mir_annotations):
            mirs.append(mir)
            mir_types.append(mir_annotations[mir])
            mss_list.append(mss[mir])
     
        if len(mirs) > 0:
            closest_mss = lib.genomic.get_closest_tss(mss_list)
        
            ret.append(";".join(mirs)) 
            ret.append(";".join(mir_types))
            ret.append(closest_mss)
            ret.append(";".join(mss_list))
        else:
            ret.append("n/a")
            ret.append("n/a")
            ret.append("n/a")
            ret.append("n/a")
        
        #
        # Modules
        #
        
        #repetitive_classification = self.repetitive.annotate(location)
        #ret.append(",".join(sorted(repetitive_classification)))
        
        for module in self.annotation_modules:
            ret.append(module.annotate(location))
        
        # Finish the annotation
        return ret
