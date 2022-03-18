import sys
import collections
import re

import lib.annotation
import lib.expression
import lib.tss
import lib.human.tss
import lib.human.annotation

import lib.genomic
import lib.headings
import lib.text
import lib.sample
import lib.genes
import lib.human.genomic

class GeneOrientatedPeaks(object):
  """
  Core annotation for gene oriented peaks
  """
  
  def __init__(self, type):
    self.type = type
    
    self.expression_list = []
    self.expression_list_headers = []

    #self.affy_gene_cb_vs_m_expression = lib.expression.AffyGeneCBvsMExpression() 
    #self.affy_gene_cb_vs_n_expression = lib.expression.AffyGeneCBvsNExpression() 
    #self.rna_gene_cb_vs_m_expression = lib.expression.RnaSeqGeneCBvsMExpression() 
    #self.rna_gene_cb_vs_n_expression = lib.expression.RnaSeqGeneCBvsNExpression()
  
    self.expression_list_headers.append("GEP Affy GCvsN")
    self.expression_list_headers.append("GEP Affy GCvsM")
    self.expression_list.append(lib.expression.AffyGeneCBvsNExpression())
    self.expression_list.append(lib.expression.AffyGeneCBvsMExpression())
  
    self.expression_list_headers.append("RNA-seq GCvsN")
    self.expression_list_headers.append("RNA-seq GCvsM")
    self.expression_list.append(lib.expression.RnaSeqGeneCBvsNExpression())
    self.expression_list.append(lib.expression.RnaSeqGeneCBvsMExpression())
    
    
    
    
    # RNA-seq from deseq2
    #self.expression_list_headers.append("RNA-seq GCvsN DeSeq2")
    #self.expression_list_headers.append("RNA-seq GCvsM DeSeq2")
    #self.expression_list.append(lib.expression.RnaSeqGeneCBvsNDeSeq2Expression())
    #self.expression_list.append(lib.expression.RnaSeqGeneCBvsMDeSeq2Expression())
    
    #
    # SUD10 stuff
    #
    self.expression_list_headers.append("RNA-seq SUD10 WTvsD83V")
    self.expression_list_headers.append("RNA-seq SUD10 WTvsSTOP")
    self.expression_list_headers.append("RNA-seq SUD10 D83VvsSTOP")
    self.expression_list.append(lib.expression.RnaSeqGeneSUD10WTvsD83vExpression())
    self.expression_list.append(lib.expression.RnaSeqGeneSUD10WTvsStopExpression())
    self.expression_list.append(lib.expression.RnaSeqGeneSUD10D83vsStopExpression())
    

    
    
    #
    # Mouse MEF2B vs WT
    #
    
    self.expression_list_headers.append("GEP Affy MEF2B Mouse WTvsKO")
    self.expression_list.append(lib.expression.GEPMEF2BMouseWTvsKOExpression())
    
    #
    # miR annotations
    #
    
    self.solid_mir_expression = lib.expression.SolidMirExpression()
    self.solid_mir_cb_vs_n_expression = lib.expression.SolidMirCBvsNExpression()
    self.agilent_mir_cb_vs_m_expression = lib.expression.AgilentMirCBvsMExpression()
    self.agilent_mir_cb_vs_n_expression = lib.expression.AgilentMirCBvsNExpression()

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
    header = f.readline().strip().split("\t")
  
    location_column = lib.text.find_index(header, lib.headings.LOCATION)
    entrez_column = lib.text.find_index(header, lib.headings.ENTREZ_ID)
    refseq_column = lib.text.find_index(header, lib.headings.REFSEQ_ID)
    symbol_column = lib.text.find_index(header, lib.headings.GENE_SYMBOL)
    type_column = lib.text.find_index(header, "Relative To Gene")
    #p_column = lib.text.find_index(header, lib.headings.P_VALUE)
    #score_column = lib.text.find_index(header, lib.headings.SCORE)
    tss_column = lib.text.find_index(header, lib.headings.TSS_DISTANCE)
    centromere_column = lib.text.find_index(header, lib.headings.CENTROMERE)

    mir_column = lib.text.find_index(header, lib.headings.MIR_SYMBOL)
    mir_type_column = lib.text.find_index(header, "Relative To miR")
    mss_column = lib.text.find_index(header, "mIR Start Distance")
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = ls.split("\t")

      type = tokens[type_column]
      p = lib.genes.find_best_p_value(header, tokens) #float(tokens[p_column])
      score = lib.genes.find_best_score(header, tokens) #= float(tokens[score_column])
      location = tokens[location_column]
      
      entrezes = tokens[entrez_column].split(";")
      symbols = tokens[symbol_column].split(";")
      refseqs = tokens[refseq_column].split(";")
      tsses = tokens[tss_column].split(";")
      centromere = tokens[centromere_column]
      
      for i in range(0, len(refseqs)):
        # The core id identifies the unique genes on a peak, rather than
        #        
        entrez = entrezes[i]
        symbol = symbols[i]
        refseq = refseqs[i]
        tss = tsses[i]

        if refseq != "n/a":
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
    
      if mir != "n/a":
        self.mirs.add(mir)
        self.collapsed_types[mir].append(tokens[mir_type_column])
        self.collapsed_p[mir].append(p)
        self.collapsed_scores[mir].append(score)
        self.collapsed_locations[mir].append(location)
        self.collapsed_tss[mir].append(tokens[mss_column])
        self.collapsed_centromeres[mir].append(centromere)
    
    f.close()
    
    
  def print_header(self):
    sys.stdout.write(lib.headings.REFSEQ_ID)
    sys.stdout.write("\t" + lib.headings.ENTREZ_ID)
    sys.stdout.write("\t" + lib.headings.GENE_SYMBOL)
    #sys.stdout.write("\tGEP Affy GCvsN")
    #sys.stdout.write("\tGEP Affy GCvsM")
    #sys.stdout.write("\tRNAseq GCvsN")
    #sys.stdout.write("\tRNAseq GCvsM")
    
    for header in self.expression_list_headers:
      sys.stdout.write("\t" + header)

    sys.stdout.write("\t" + self.type + " Relative To Gene")
    sys.stdout.write("\t" + self.type + " TSS Closest Distance")
    sys.stdout.write("\t" + self.type + " " + lib.headings.TSS_DISTANCE)
    sys.stdout.write("\t" + self.type + " " + lib.headings.CENTROMERE)
    sys.stdout.write("\tmiR Symbol")
    sys.stdout.write("\tmiREP Agilent GCvsN")
    sys.stdout.write("\tmiREP Agilent GCvsM")
    sys.stdout.write("\tsRE GCvsNM")
    sys.stdout.write("\tsRE GCvsN")
    sys.stdout.write("\t" + self.type + " Relative To miR")
    sys.stdout.write("\t" + self.type + " miR Start Closest Distance")
    sys.stdout.write("\t" + self.type + " miR Start Distance")
    sys.stdout.write("\tBest P-value (ChIPseeqer)")
    sys.stdout.write("\tBest Score (ChIPseeqer)")
    sys.stdout.write("\t" + self.type + " Count")
    sys.stdout.write("\t" + self.type + " Genomic Locations (hg19)")
    sys.stdout.write("\n");

  
  def get_ids(self):
    return sorted(self.refseqs)
    
    
  def get_mirs(self):
    return sorted(self.mirs)
    

  def gene_orient_peak(self, id):
    entrez = self.collapsed_entrezes[id]
    
    sys.stdout.write(id)
    sys.stdout.write("\t" + self.collapsed_entrezes[id])
    sys.stdout.write("\t" + self.collapsed_symbols[id])
      
    #sys.stdout.write("\t" + self.affy_gene_cb_vs_n_expression.get_expression(entrez))
    #sys.stdout.write("\t" + self.affy_gene_cb_vs_m_expression.get_expression(entrez))
    #sys.stdout.write("\t" + self.rna_gene_cb_vs_n_expression.get_expression(entrez))
    #sys.stdout.write("\t" + self.rna_gene_cb_vs_m_expression.get_expression(entrez))

    for expression in self.expression_list:
      sys.stdout.write("\t" + expression.get_expression(entrez))

    sys.stdout.write("\t" + ";".join(self.collapsed_types[id]))
    
    # if there are some nearest tss, print the closest
    sys.stdout.write("\t" + lib.genomic.get_closest_tss(self.collapsed_tss[id]))
          
    sys.stdout.write("\t" + ";".join(self.collapsed_tss[id]))
    
    # Centromeres
    sys.stdout.write("\t" + ";".join(self.collapsed_centromeres[id]))
    
    # no mir symbol
    sys.stdout.write("\tn/a")
    
    # no agilent expression
    sys.stdout.write("\tn/a")
    sys.stdout.write("\tn/a")    
    
    #no small rna    
    sys.stdout.write("\tn/a")
    sys.stdout.write("\tn/a")
    
    #no peak relative to mir
    sys.stdout.write("\tn/a")
    sys.stdout.write("\tn/a")
    sys.stdout.write("\tn/a")
      
    # pick the smallest p
    p = sorted(self.collapsed_p[id])
    sys.stdout.write("\t" + str(p[0]))
    
    # pick the largest score
    scores = sorted(self.collapsed_scores[id], reverse=True)
    sys.stdout.write("\t" + str(scores[0]))
      
    # peak count
    sys.stdout.write("\t" + str(len(self.collapsed_locations[id])))
      
    sys.stdout.write("\t" + ";".join(self.collapsed_locations[id]))
    
    sys.stdout.write("\n");
    
    
  def mir_orient_peak(self, mir):
    sys.stdout.write("n/a")

    for i in range(0, len(self.expression_list_headers)):
      sys.stdout.write("\tn/a")

    # Fill in the gap
    for i in range(0, 5):
      sys.stdout.write("\tn/a")

    sys.stdout.write("\t" + ";".join(self.collapsed_centromeres[mir]))
    sys.stdout.write("\t" + mir)
    sys.stdout.write("\t" + self.agilent_mir_cb_vs_n_expression.get_expression(mir))
    sys.stdout.write("\t" + self.agilent_mir_cb_vs_m_expression.get_expression(mir))
    sys.stdout.write("\t" + self.solid_mir_expression.get_expression(mir))
    sys.stdout.write("\t" + self.solid_mir_cb_vs_n_expression.get_expression(mir))
    sys.stdout.write("\t" + ";".join(self.collapsed_types[mir]))
    sys.stdout.write("\t" + lib.genomic.get_closest_tss(self.collapsed_tss[mir]))
    sys.stdout.write("\t" + ";".join(self.collapsed_tss[mir]))
    
    # pick the smallest p
    p = sorted(self.collapsed_p[mir])
    sys.stdout.write("\t" + str(p[0]))
    
    # pick the largest score
    scores = sorted(self.collapsed_scores[mir], reverse=True)
    sys.stdout.write("\t" + str(scores[0]))
    
    sys.stdout.write("\t" + str(len(self.collapsed_locations[mir])))
    sys.stdout.write("\t" + ";".join(self.collapsed_locations[mir]))
    
    sys.stdout.write("\n");
  
  def print_log(self):
    """
    Print a log to indicate what was used for annotations.
    """
    
    f = open("genes.log", "w")
    
    f.write("Expression Type\tSource\n")
    
    for i in range(0, len(self.expression_list_headers)):
      f.write(self.expression_list_headers[i] + "\t" + self.expression_list[i].get_file() + "\n")
      
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
    header = f.readline().strip().split("\t")

    # independent columns
    location_column = lib.text.find_index(header, lib.headings.LOCATION)
    p_column = lib.text.find_index(header, lib.headings.P_VALUE)
    score_column = lib.text.find_index(header, lib.headings.SCORE)
    centromere_column = lib.text.find_index(header, lib.headings.CENTROMERE)

    entrez_column = lib.text.find_index(header, lib.headings.ENTREZ_ID)
    refseq_column = lib.text.find_index(header, lib.headings.REFSEQ_ID)
    symbol_column = lib.text.find_index(header, lib.headings.GENE_SYMBOL)
    type_column = lib.text.find_index(header, lib.headings.RELATIVE)
    tss_column = lib.text.find_index(header, lib.headings.TSS_DISTANCE)
    
    closest_entrez_column = lib.text.find_index(header, lib.headings.CLOSEST_ENTREZ_ID)
    closest_refseq_column = lib.text.find_index(header, lib.headings.CLOSEST_REFSEQ_ID)
    closest_symbol_column = lib.text.find_index(header, lib.headings.CLOSEST_GENE_SYMBOL)
    closest_type_column = lib.text.find_index(header, lib.headings.CLOSEST_RELATIVE)
    closest_tss_column = lib.text.find_index(header, lib.headings.CLOSEST_TSS_DISTANCE)
    
    
    mir_column = lib.text.find_index(header, lib.headings.MIR_SYMBOL)
    mir_type_column = lib.text.find_index(header, "Relative To miR")
    mss_column = lib.text.find_index(header, "mIR Start Distance")
    
    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
    
      tokens = ls.split("\t")

      location = tokens[location_column]
      p = float(tokens[p_column])
      score = float(tokens[score_column])
      
      centromere = tokens[centromere_column]
      
      if tokens[entrez_column] != lib.text.NA:
        # Preference is given to peaks within a gene
        entrezes = tokens[entrez_column].split(";")
        symbols = tokens[symbol_column].split(";")
        refseqs = tokens[refseq_column].split(";")
        tsses = tokens[tss_column].split(";")
        type = tokens[type_column]
      else:
        entrezes = tokens[closest_entrez_column].split(";")
        symbols = tokens[closest_symbol_column].split(";")
        refseqs = tokens[closest_refseq_column].split(";")
        tsses = tokens[closest_tss_column].split(";")
        type = tokens[closest_type_column]
      
      
        
      for i in range(0, len(refseqs)):
        # The core id identifies the unique genes on a peak, rather than
        #        
        entrez = entrezes[i]
        symbol = symbols[i]
        refseq = refseqs[i]
        tss = tsses[i]

        if refseq != lib.text.NA:
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
    
      if mir != lib.text.NA:
        self.mirs.add(mir)
        self.collapsed_types[mir].append(tokens[mir_type_column])
        self.collapsed_p[mir].append(p)
        self.collapsed_scores[mir].append(score)
        self.collapsed_locations[mir].append(location)
        self.collapsed_tss[mir].append(tokens[mss_column])
        self.collapsed_centromeres[mir].append(centromere)
    
    f.close()


class RefSeqGenes(lib.genes.RefSeqGenes):
    def __init__(self):
        super().__init__(lib.human.annotation.REFSEQ_FILE)


class AnnotatePeak(lib.genes.AnnotatePeak):
    """
    Core annotation for annotating peaks/regions
    """
  
    def __init__(self, 
                 type, 
                 prom_ext_5p=5000, 
                 prom_ext_3p=4000, 
                 bin_size=10000):
        super().__init__(type,
             lib.human.tss.RefSeqAnnotation(prom_ext_5p, prom_ext_3p, bin_size),
             RefSeqGenes(),
             lib.human.tss.RefSeqTss(prom_ext_5p, prom_ext_3p),
             lib.human.tss.RefSeqEnd(prom_ext_5p, prom_ext_3p),
             prom_ext_5p,
             prom_ext_3p,
             bin_size)
        
        # annotations specific to human
        self.annotation_modules.append(lib.human.genomic.SimpleTandemRepeats())
        self.annotation_modules.append(lib.human.genomic.EncodeBlacklist())
        self.annotation_modules.append(lib.human.genomic.GiuliaBlacklist())
