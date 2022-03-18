import sys
import collections
import re
from typing import Mapping, Union

import pychipseq.annotation
import pychipseq.expression
import pychipseq.tss
import pychipseq.human.tss
import pychipseq.human.annotation

import pychipseq.genomic
import pychipseq.headings
import pychipseq.text
import pychipseq.sample
import pychipseq.peaks
import pychipseq.genes
import pychipseq.human.genomic
import pychipseq.tad
import pychipseq.human.tad
import pychipseq.human.mir


class GeneOrientatedPeaks(object):
    """
    Core annotation for gene oriented peaks
    """

    def __init__(self, type):
        self.type = type

        self.expression_list = []
        self.expression_list_headers = []

        #self.affy_gene_cb_vs_m_expression = pychipseq.expression.AffyGeneCBvsMExpression()
        #self.affy_gene_cb_vs_n_expression = pychipseq.expression.AffyGeneCBvsNExpression()
        #self.rna_gene_cb_vs_m_expression = pychipseq.expression.RnaSeqGeneCBvsMExpression()
        #self.rna_gene_cb_vs_n_expression = pychipseq.expression.RnaSeqGeneCBvsNExpression()

        self.expression_list_headers.append('GEP Affy GCvsN')
        self.expression_list_headers.append('GEP Affy GCvsM')
        self.expression_list.append(
            pychipseq.expression.AffyGeneCBvsNExpression())
        self.expression_list.append(
            pychipseq.expression.AffyGeneCBvsMExpression())

        self.expression_list_headers.append('RNA-seq GCvsN')
        self.expression_list_headers.append('RNA-seq GCvsM')
        self.expression_list.append(
            pychipseq.expression.RnaSeqGeneCBvsNExpression())
        self.expression_list.append(
            pychipseq.expression.RnaSeqGeneCBvsMExpression())

        # RNA-seq from deseq2
        #self.expression_list_headers.append('RNA-seq GCvsN DeSeq2')
        #self.expression_list_headers.append('RNA-seq GCvsM DeSeq2')
        # self.expression_list.append(pychipseq.expression.RnaSeqGeneCBvsNDeSeq2Expression())
        # self.expression_list.append(pychipseq.expression.RnaSeqGeneCBvsMDeSeq2Expression())

        #
        # SUD10 stuff
        #
        self.expression_list_headers.append('RNA-seq SUD10 WTvsD83V')
        self.expression_list_headers.append('RNA-seq SUD10 WTvsSTOP')
        self.expression_list_headers.append('RNA-seq SUD10 D83VvsSTOP')
        self.expression_list.append(
            pychipseq.expression.RnaSeqGeneSUD10WTvsD83vExpression())
        self.expression_list.append(
            pychipseq.expression.RnaSeqGeneSUD10WTvsStopExpression())
        self.expression_list.append(
            pychipseq.expression.RnaSeqGeneSUD10D83vsStopExpression())

        #
        # Mouse MEF2B vs WT
        #

        self.expression_list_headers.append('GEP Affy MEF2B Mouse WTvsKO')
        self.expression_list.append(
            pychipseq.expression.GEPMEF2BMouseWTvsKOExpression())

        #
        # miR annotations
        #

        self.solid_mir_expression = pychipseq.expression.SolidMirExpression()
        self.solid_mir_cb_vs_n_expression = pychipseq.expression.SolidMirCBvsNExpression()
        self.agilent_mir_cb_vs_m_expression = pychipseq.expression.AgilentMirCBvsMExpression()
        self.agilent_mir_cb_vs_n_expression = pychipseq.expression.AgilentMirCBvsNExpression()

        #
        # Sets for annotations
        #

        self.collapsed_entrezes = collections.defaultdict(str)
        self.collapsed_symbols = collections.defaultdict(str)
        self.collapsed_p = collections.defaultdict(list)
        self.collapsed_scores = collections.defaultdict(list)
        self.collapsed_locations = collections.defaultdict(list)
        self.collapsed_tss = collections.defaultdict(list)
        self.collapsed_types = collections.defaultdict(list)
        self.collapsed_centromeres = collections.defaultdict(list)

        self.mirs = set()
        self.refseqs = set()

        #self.collapsed_mir_types = collections.defaultdict(list)
        #self.collapsed_mirs_p = collections.defaultdict(list)
        #self.collapsed_mirs_scores = collections.defaultdict(list)
        #self.collapsed_mirs_locations = collections.defaultdict(list)
        #self.collapsed_mss = collections.defaultdict(list)

    def add_expression(self, name, expression):
        self.expression_list_headers.append(name)
        self.expression_list.append(expression)

    def load_annotations(self, file):
        f = open(file, 'r')

        # skip header
        header = f.readline().strip().split('\t')

        location_column = pychipseq.text.find_index(
            header, pychipseq.headings.LOCATION)
        entrez_column = pychipseq.text.find_index(
            header, pychipseq.headings.ENTREZ_ID)
        refseq_column = pychipseq.text.find_index(
            header, pychipseq.headings.REFSEQ_ID)
        symbol_column = pychipseq.text.find_index(
            header, pychipseq.headings.GENE_SYMBOL)
        type_column = pychipseq.text.find_index(header, 'Relative To Gene')
        #p_column = pychipseq.text.find_index(header, pychipseq.headings.P_VALUE)
        #score_column = pychipseq.text.find_index(header, pychipseq.headings.SCORE)
        tss_column = pychipseq.text.find_index(
            header, pychipseq.headings.TSS_DISTANCE)
        centromere_column = pychipseq.text.find_index(
            header, pychipseq.headings.CENTROMERE)

        mir_column = pychipseq.text.find_index(
            header, pychipseq.headings.MIR_SYMBOL)
        mir_type_column = pychipseq.text.find_index(header, 'Relative To miR')
        mss_column = pychipseq.text.find_index(header, 'mIR Start Distance')

        for line in f:
            ls = line.strip()

            if len(ls) == 0:
                continue

            tokens = ls.split('\t')

            type = tokens[type_column]
            p = pychipseq.genes.find_best_p_value(
                header, tokens)  # float(tokens[p_column])
            score = pychipseq.genes.find_best_score(
                header, tokens)  # = float(tokens[score_column])
            location = tokens[location_column]

            entrezes = tokens[entrez_column].split(';')
            symbols = tokens[symbol_column].split(';')
            refseqs = tokens[refseq_column].split(';')
            tsses = tokens[tss_column].split(';')
            centromere = tokens[centromere_column]

            for i in range(0, len(refseqs)):
                # The core id identifies the unique genes on a peak, rather than
                #
                entrez = entrezes[i]
                symbol = symbols[i]
                refseq = refseqs[i]
                tss = tsses[i]

                if refseq != pychipseq.text.NA:
                    self.refseqs.add(refseq)
                    self.collapsed_entrezes[refseq] = entrez
                    self.collapsed_symbols[refseq] = symbol
                    self.collapsed_p[refseq].append(p)
                    self.collapsed_scores[refseq].append(score)
                    self.collapsed_locations[refseq].append(location)
                    self.collapsed_tss[refseq].append(tss)
                    self.collapsed_types[refseq].append(type)
                    self.collapsed_centromeres[refseq].append(centromere)

            mir = tokens[mir_column]

            if mir != pychipseq.text.NA:
                self.mirs.add(mir)
                self.collapsed_types[mir].append(tokens[mir_type_column])
                self.collapsed_p[mir].append(p)
                self.collapsed_scores[mir].append(score)
                self.collapsed_locations[mir].append(location)
                self.collapsed_tss[mir].append(tokens[mss_column])
                self.collapsed_centromeres[mir].append(centromere)

        f.close()

    def print_header(self):
        print('\t'.join(self.get_header()))

    def get_header(self):
        ret = []
        ret.append(pychipseq.headings.REFSEQ_ID)
        ret.append(pychipseq.headings.ENTREZ_ID)
        ret.append(pychipseq.headings.GENE_SYMBOL)
        #ret.append(f'\tGEP Affy GCvsN")
        #ret.append(f'\tGEP Affy GCvsM")
        #ret.append(f'\tRNAseq GCvsN")
        #ret.append(f'\tRNAseq GCvsM")

        for header in self.expression_list_headers:
            ret.append(header)

        ret.append(f'{self.type} Relative To Gene')
        ret.append(f'{self.type} TSS Closest Distance')
        ret.append(f'{self.type} {pychipseq.headings.TSS_DISTANCE}')
        ret.append(f'{self.type} {pychipseq.headings.CENTROMERE}')
        ret.append('miR Symbol')
        ret.append('miREP Agilent GCvsN')
        ret.append('miREP Agilent GCvsM')
        ret.append('sRE GCvsNM')
        ret.append('sRE GCvsN')
        ret.append(f'{self.type} Relative To miR')
        ret.append(f'{self.type} miR Start Closest Distance')
        ret.append(f'{self.type} miR Start Distance')
        ret.append('Best P-value (ChIPseeqer)')
        ret.append('Best Score (ChIPseeqer)')
        ret.append('{self.type} Count')
        ret.append('{self.type} Genomic Locations (hg19)')
        
        return ret

    def get_ids(self):
        return sorted(self.refseqs)

    def get_mirs(self):
        return sorted(self.mirs)

    def gene_orient_peak(self, id):
        entrez = self.collapsed_entrezes[id]

        ret = []

        ret.append(id)
        ret.append(self.collapsed_entrezes[id])
        ret.append(self.collapsed_symbols[id])

        #ret.append(f'\t{self.affy_gene_cb_vs_n_expression.get_expression(entrez))
        #ret.append(f'\t{self.affy_gene_cb_vs_m_expression.get_expression(entrez))
        #ret.append(f'\t{self.rna_gene_cb_vs_n_expression.get_expression(entrez))
        #ret.append(f'\t{self.rna_gene_cb_vs_m_expression.get_expression(entrez))

        for expression in self.expression_list:
            ret.append(expression.get_expression(entrez))

        ret.append(';'.join(self.collapsed_types[id]))

        # if there are some nearest tss, print the closest
        ret.append(
            pychipseq.genomic.get_closest_tss(self.collapsed_tss[id]))

        ret.append(';'.join(self.collapsed_tss[id]))

        # Centromeres
        ret.append(';'.join(self.collapsed_centromeres[id]))

        # no mir symbol
        ret.append(pychipseq.text.NA)

        # no agilent expression
        ret.append(pychipseq.text.NA)
        ret.append(pychipseq.text.NA)

        # no small rna
        ret.append(pychipseq.text.NA)
        ret.append(pychipseq.text.NA)

        # no peak relative to mir
        ret.append(pychipseq.text.NA)
        ret.append(pychipseq.text.NA)
        ret.append(pychipseq.text.NA)

        # pick the smallest p
        p = sorted(self.collapsed_p[id])
        ret.append(p[0])

        # pick the largest score
        scores = sorted(self.collapsed_scores[id], reverse=True)
        ret.append(scores[0])

        # peak count
        ret.append(len(self.collapsed_locations[id]))

        ret.append(';'.join(self.collapsed_locations[id]))

        print('\t'.join([str(x) for x in ret]))

        return ret

    def mir_orient_peak(self, mir):
        ret = []

        ret.append(pychipseq.text.NA)

        ret.extend([pychipseq.text.NA] * len(self.expression_list_headers))

        # Fill in the gap
        ret.extend([pychipseq.text.NA] * 5)

        ret.append(';'.join(self.collapsed_centromeres[mir]))
        ret.append(mir)
        ret.append(
            self.agilent_mir_cb_vs_n_expression.get_expression(mir))
        ret.append(
            self.agilent_mir_cb_vs_m_expression.get_expression(mir))
        ret.append(self.solid_mir_expression.get_expression(mir))
        ret.append(self.solid_mir_cb_vs_n_expression.get_expression(mir))
        ret.append(';'.join(self.collapsed_types[mir]))
        ret.append(
            pychipseq.genomic.get_closest_tss(self.collapsed_tss[mir]))
        ret.append(';'.join(self.collapsed_tss[mir]))

        # pick the smallest p
        p = sorted(self.collapsed_p[mir])
        ret.append(p[0])

        # pick the largest score
        scores = sorted(self.collapsed_scores[mir], reverse=True)
        ret.append(scores[0])

        ret.append(len(self.collapsed_locations[mir]))
        ret.append(';'.join(self.collapsed_locations[mir]))

        print('\t'.join([str(x) for x in ret]))

        return ret

    def print_log(self):
        """
        Print a log to indicate what was used for annotations.
        """

        f = open('genes.log', 'w')

        f.write('Expression Type\tSource\n')

        for i in range(0, len(self.expression_list_headers)):
            f.write(
                f'{self.expression_list_headers[i]}\t{self.expression_list[i].get_file()}\n')

        f.close()


class ClosestGeneOrientatedPeaks(GeneOrientatedPeaks):
    """
    Produce a gene oriented version of peaks using the peaks within
    a gene/promoter first and then the intergenic ones last. This is
    a hybrid of the two annotations with preference given to 
    peaks overlapping genes
    """

    def load_annotations(self, file):
        f = open(file, 'r')

        # skip header
        header = f.readline().strip().split('\t')

        # independent columns
        location_column = pychipseq.text.find_index(
            header, pychipseq.headings.LOCATION)
        p_column = pychipseq.text.find_index(
            header, pychipseq.headings.P_VALUE)
        score_column = pychipseq.text.find_index(
            header, pychipseq.headings.SCORE)
        centromere_column = pychipseq.text.find_index(
            header, pychipseq.headings.CENTROMERE)

        entrez_column = pychipseq.text.find_index(
            header, pychipseq.headings.ENTREZ_ID)
        refseq_column = pychipseq.text.find_index(
            header, pychipseq.headings.REFSEQ_ID)
        symbol_column = pychipseq.text.find_index(
            header, pychipseq.headings.GENE_SYMBOL)
        type_column = pychipseq.text.find_index(
            header, pychipseq.headings.RELATIVE)
        tss_column = pychipseq.text.find_index(
            header, pychipseq.headings.TSS_DISTANCE)

        closest_entrez_column = pychipseq.text.find_index(
            header, pychipseq.headings.CLOSEST_ENTREZ_ID)
        closest_refseq_column = pychipseq.text.find_index(
            header, pychipseq.headings.CLOSEST_REFSEQ_ID)
        closest_symbol_column = pychipseq.text.find_index(
            header, pychipseq.headings.CLOSEST_GENE_SYMBOL)
        closest_type_column = pychipseq.text.find_index(
            header, pychipseq.headings.CLOSEST_RELATIVE)
        closest_tss_column = pychipseq.text.find_index(
            header, pychipseq.headings.CLOSEST_TSS_DISTANCE)

        mir_column = pychipseq.text.find_index(
            header, pychipseq.headings.MIR_SYMBOL)
        mir_type_column = pychipseq.text.find_index(header, 'Relative To miR')
        mss_column = pychipseq.text.find_index(header, 'mIR Start Distance')

        for line in f:
            ls = line.strip()

            if len(ls) == 0:
                continue

            tokens = ls.split('\t')

            location = tokens[location_column]
            p = float(tokens[p_column])
            score = float(tokens[score_column])

            centromere = tokens[centromere_column]

            if tokens[entrez_column] != pychipseq.text.NA:
                # Preference is given to peaks within a gene
                entrezes = tokens[entrez_column].split(';')
                symbols = tokens[symbol_column].split(';')
                refseqs = tokens[refseq_column].split(';')
                tsses = tokens[tss_column].split(';')
                type = tokens[type_column]
            else:
                entrezes = tokens[closest_entrez_column].split(';')
                symbols = tokens[closest_symbol_column].split(';')
                refseqs = tokens[closest_refseq_column].split(';')
                tsses = tokens[closest_tss_column].split(';')
                type = tokens[closest_type_column]

            for i in range(0, len(refseqs)):
                # The core id identifies the unique genes on a peak, rather than
                #
                entrez = entrezes[i]
                symbol = symbols[i]
                refseq = refseqs[i]
                tss = tsses[i]

                if refseq != pychipseq.text.NA:
                    self.refseqs.add(refseq)
                    self.collapsed_entrezes[refseq] = entrez
                    self.collapsed_symbols[refseq] = symbol
                    self.collapsed_p[refseq].append(p)
                    self.collapsed_scores[refseq].append(score)
                    self.collapsed_locations[refseq].append(location)
                    self.collapsed_tss[refseq].append(tss)
                    self.collapsed_types[refseq].append(type)
                    self.collapsed_centromeres[refseq].append(centromere)

            mir = tokens[mir_column]

            if mir != pychipseq.text.NA:
                self.mirs.add(mir)
                self.collapsed_types[mir].append(tokens[mir_type_column])
                self.collapsed_p[mir].append(p)
                self.collapsed_scores[mir].append(score)
                self.collapsed_locations[mir].append(location)
                self.collapsed_tss[mir].append(tokens[mss_column])
                self.collapsed_centromeres[mir].append(centromere)

        f.close()


class RefSeqGenes(pychipseq.genes.RefSeqGenes):
    def __init__(self):
        super().__init__(pychipseq.human.annotation.REFSEQ_FILE)


REFSEQ_GENES = RefSeqGenes()


class AnnotatePeak(pychipseq.peaks.AnnotatePeak):
    """
    Core annotation for annotating peaks/regions
    """

    def __init__(self,
                 type,
                 prom_ext_5p: int = 5000,
                 prom_ext_3p: int = 4000,
                 bin_size: int = 10000,
                 n_closest: int = 5):
        super().__init__(type,
                         pychipseq.human.tss.RefSeqAnnotation(
                             prom_ext_5p, prom_ext_3p, bin_size),
                         REFSEQ_GENES,
                         prom_ext_5p,
                         prom_ext_3p)
        # annotations specific to human

        # default closest gene
        refseq_start = pychipseq.human.tss.RefSeqStart(
            REFSEQ_GENES, prom_ext_5p, prom_ext_3p)

        self.add_module(refseq_start)
        self.add_module(pychipseq.human.tss.RefSeqEnd(
            REFSEQ_GENES, prom_ext_5p, prom_ext_3p))

        self.add_module(
            pychipseq.peaks.NClosestGenes(REFSEQ_GENES, refseq_start, n=n_closest))

        self.add_module(pychipseq.human.mir.MirAnnotation(
            REFSEQ_GENES, prom_ext_5p, bin_size))
        self.add_module(pychipseq.human.genomic.Repetitive())
        self.add_module(pychipseq.human.genomic.SimpleTandemRepeats())
        self.add_module(pychipseq.human.genomic.EncodeBlacklist())
        self.add_module(pychipseq.human.genomic.GiuliaBlacklist())
        self.add_module(pychipseq.human.tad.GencodeTADAnnotation())
        self.add_module(pychipseq.tad.IsTADAnnotation())
        self.add_module(pychipseq.tad.IsClosestTADAnnotation())
        self.add_module(pychipseq.human.tss.OverlapTss())
