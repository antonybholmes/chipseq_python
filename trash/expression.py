import sys
import collections
import re


def gene_exp(genes_file, up_gene_file, down_gene_file):
  """
  Takes three gene lists for mapping genes to whether they are up
  or down regulated in a particular gene expression type e.g.
  """

  gene_exp = collections.defaultdict(str)
  
  load_gene_exp(genes_file, gene_exp, "not_moving")
  load_gene_exp(up_gene_file, gene_exp, "up")
  load_gene_exp(down_gene_file, gene_exp, "down")


def load_gene_exp(file, genes_map, type):
  
  f = open(file, "r")
  
  for line in f:
    ls = line.strip()
    
    if len(ls) == 0:
      continue
    
    tokens = ls.split("\t");
    
    if len(tokens[0]) == 0 or re.match(r'n\/a', tokens[0]):
      continue
    
    genes_map[tokens[0].lower()] = type
  
  f.close()


def load_mir_exp(file, genes_map, type):
  
  f = open(file, "r")
  
  for line in f:
    ls = line.strip()
    
    if len(ls) == 0:
      continue
    
    tokens = ls.split("\t");
    
    if len(tokens[0]) == 0 or re.match(r'n\/a', tokens[0]):
      continue
    
    genes_map[get_mir_id(tokens[0])] = type
  
  f.close()
  
  
class Expression(object):
  def __init__(self, dir):  
    self.expression_map = collections.defaultdict(str)
    
    load_gene_exp(dir + "not_moving_genes.txt",
      self.expression_map,
      "not_moving")
    
    # now mark as down-regulated
    load_gene_exp(dir + "down_genes.txt",
      self.expression_map,
      "down")
    
    # now mark as up-regulated
    load_gene_exp(dir + "up_genes.txt",
      self.expression_map,
      "up")
  
  def get_expression(self, id):
    s = id.lower()
    
    #sys.stderr.write("get exp for " + s + "\n")
    
    if s not in self.expression_map:
      return "n/a"
      
    return self.expression_map[s]
    

class MirExpression(object):
  def __init__(self):  
    self.expression_map = collections.defaultdict(str)
    self.previous_ids = MirPreviousIds()
  
  def get_expression(self, id):
    s = get_mir_id(id)

    if s in self.expression_map:
      return self.expression_map[s]
    
    previous_ids = self.previous_ids.get_previous_ids(id)

    if re.match(r'.*\-151.*', id):
      sys.stderr.write(id + " " + ";".join(previous_ids) + "\n")
    
    for pid in previous_ids:
      if pid in self.expression_map:
        return self.expression_map[pid]
    
    #s = get_mir_precursor_id(id)
    
    #if s in self.expression_map:
    #  return self.expression_map[s]
      
    return "n/a"
    

class AffyGeneCTRLvsSIExpression(Expression):
  """
  Annotate whether a gene is up or down regulated in affymetrix data
  """  
  
  def __init__(self):
    super(AffyGeneCTRLvsSIExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/bcl6_silencing/")
      

class AffyGeneCBvsMExpression(Expression):
  """
  Annotate whether a gene is up or down regulated in affymetrix data
  """  
  
  def __init__(self):
    super(AffyGeneCBvsMExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/hg-u133_plus_2_cb_vs_m/")


class AffyGeneCBvsNExpression(Expression):
  """
  Annotate whether a gene is up or down regulated in affymetrix data
  """  
  
  def __init__(self):
    super(AffyGeneCBvsNExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/hg-u133_plus_2_cb_vs_n/")
      
class RnaSeqExpression(Expression):
  def __init__(self, dir):
    super(RnaSeqExpression, self).__init__(dir)
    
    load_gene_exp(dir + "not_expressed_genes.txt",
      self.expression_map,
      "no_exp")

      
class RnaSeqGeneCBvsMExpression(RnaSeqExpression):
  def __init__(self):
    super(RnaSeqGeneCBvsMExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq_cb_vs_m/")


class RnaSeqGeneCBvsNExpression(RnaSeqExpression):
  def __init__(self):
    super(RnaSeqGeneCBvsNExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/rna_seq_cb_vs_n/")
      

      
class DzLzGeneExpression(Expression):
  """
  Annotate whether a gene is up or down regulated in affymetrix data
  """  
  
  def __init__(self):
    super(DzLzGeneExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/dark_zone_vs_light_zone/")


class SolidMirExpression(MirExpression):
  """
  Annotate whether a miR is up or down regulated in our SOLiD miR data
  """
  
  def __init__(self):
    super(SolidMirExpression, self).__init__()
    
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/solid_mir/not_moving_mirs.txt",
      self.expression_map,
      "not_moving");
  
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/solid_mir/down_mirs.txt",
      self.expression_map,
      "down");

    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/solid_mir/up_mirs.txt",
      self.expression_map,
      "up");


class SolidMirLog2Expression(MirExpression):
  def __init__(self):
    super(SolidMirLog2Expression, self).__init__()
    
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/solid_mir/not_moving_mirs_log2.txt",
      self.expression_map,
      "not_moving");
  
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/solid_mir/down_mirs_log2.txt",
      self.expression_map,
      "down");

    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/solid_mir/up_mirs_log2.txt",
      self.expression_map,
      "up");
      
      
class SolidMirCBvsNExpression(MirExpression):
  """
  Annotate whether a miR is up or down regulated in our SOLiD miR data
  """
  
  def __init__(self):
    super(SolidMirCBvsNExpression, self).__init__()
    
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/solid_mir/not_moving_mirs_cb_vs_n.txt",
      self.expression_map,
      "not_moving");
  
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/solid_mir/down_mirs_cb_vs_n.txt",
      self.expression_map,
      "down");

    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/solid_mir/up_mirs_cb_vs_n.txt",
      self.expression_map,
      "up");


class SolidMirCBvsNLog2Expression(MirExpression):
  def __init__(self):
    super(SolidMirCBvsNLog2Expression, self).__init__()
    
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/solid_mir/not_moving_mirs_cb_vs_n_log2.txt",
      self.expression_map,
      "not_moving");
  
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/solid_mir/down_mirs_cb_vs_n_log2.txt",
      self.expression_map,
      "down");

    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/solid_mir/up_mirs_cb_vs_n_log2.txt",
      self.expression_map,
      "up");


class AgilentMirCBvsMExpression(MirExpression):
  """
  Annotate whether a miR is up or down regulated in our SOLiD miR data
  """
  
  def __init__(self):
    super(AgilentMirCBvsMExpression, self).__init__()
    
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/agilent_mir/not_moving_mirs_cb_vs_m.txt",
      self.expression_map,
      "not_moving");
  
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/agilent_mir/down_mirs_cb_vs_m.txt",
      self.expression_map,
      "down");

    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/agilent_mir/up_mirs_cb_vs_m.txt",
      self.expression_map,
      "up");
      

class AgilentMirCBvsNExpression(MirExpression):
  """
  Annotate whether a miR is up or down regulated in our SOLiD miR data
  """
  
  def __init__(self):
    super(AgilentMirCBvsNExpression, self).__init__()
    
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/agilent_mir/not_moving_mirs_cb_vs_n.txt",
      self.expression_map,
      "not_moving");
  
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/agilent_mir/down_mirs_cb_vs_n.txt",
      self.expression_map,
      "down");

    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/agilent_mir/up_mirs_cb_vs_n.txt",
      self.expression_map,
      "up");
      

class AgilentMirCBvsMLog2Expression(MirExpression):
  def __init__(self):
    super(AgilentMirCBvsMLog2Expression, self).__init__()
    
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/agilent_mir/not_moving_mirs_cb_vs_m_log2.txt",
      self.expression_map,
      "not_moving");
  
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/agilent_mir/down_mirs_cb_vs_m_log2.txt",
      self.expression_map,
      "down");

    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/agilent_mir/up_mirs_cb_vs_m_log2.txt",
      self.expression_map,
      "up");
      

class AgilentMirCBvsNLog2Expression(MirExpression):
  def __init__(self):
    super(AgilentMirCBvsNLog2Expression, self).__init__()
    
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/agilent_mir/not_moving_mirs_cb_vs_n_log2.txt",
      self.expression_map,
      "not_moving");
  
    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/agilent_mir/down_mirs_cb_vs_n_log2.txt",
      self.expression_map,
      "down");

    load_mir_exp("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/agilent_mir/up_mirs_cb_vs_n_log2.txt",
      self.expression_map,
      "up");
      

