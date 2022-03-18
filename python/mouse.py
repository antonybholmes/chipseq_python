# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 12:22:21 2014

@author: Antony Holmes
"""

import annotation

class AffyMouseCBPKOvsWTCg1GeneExpression(annotation.Expression):
  def __init__(self):
    super(AffyMouseCBPKOvsWTCg1GeneExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/cbp_ko_vs_wt_cg1_mouse/")


class AffyMouseCBPKOvsHETCg1GeneExpression(annotation.Expression):
  def __init__(self):
    super(AffyMouseCBPKOvsHETCg1GeneExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/cbp_ko_vs_het_cg1_mouse/")


class AffyMouseCBPHETvsWTCg1GeneExpression(annotation.Expression):
  def __init__(self):
    super(AffyMouseCBPHETvsWTCg1GeneExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/cbp_het_vs_wt_cg1_mouse/")


class AffyMouseMLL2KOvsWTCD19GeneExpression(annotation.Expression):
  """
  Annotate whether a gene is up or down regulated in affymetrix data
  """  
  
  def __init__(self):
    super(AffyMouseMLL2KOvsWTCD19GeneExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/mll2_ko_vs_wt_cd19_mouse/")
      
      
class AffyMouseMLL2KOvsWTCG1GeneExpression(annotation.Expression):
  """
  Annotate whether a gene is up or down regulated in affymetrix data
  """  
  
  def __init__(self):
    super(AffyMouseMLL2KOvsWTCG1GeneExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/mll2_ko_vs_wt_cg1_mouse/")
      
      
class AffyMouseFoxo1GeneExpression(annotation.Expression):
  """
  Annotate whether a gene is up or down regulated in affymetrix data
  """  
  
  def __init__(self):
    super(AffyMouseFoxo1GeneExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/foxo1_wt_vs_ko_mouse/")


class AffyMouseFoxo1Log2GeneExpression(annotation.Expression):
  """
  Annotate whether a gene is up or down regulated in affymetrix data
  """  
  
  def __init__(self):
    super(AffyMouseFoxo1Log2GeneExpression, self).__init__("/ifs/scratch/cancer/Lab_RDF/abh2138/references/rdf/foxo1_wt_vs_ko_mouse/")