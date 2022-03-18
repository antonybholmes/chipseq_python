# -*- coding: utf-8 -*-
"""
Functions related to samples

Created on Sat Jan 31 16:46:42 2015

@author: antony
"""

import sys
import re

import lib.genomic


def get_sample_id(text):
  """
  Get the unique RK id of a sample.
  """
  
  sys.stderr.write(text + "\n")
  
  # Id consists of two letters and a three digit number
  return re.match(r'.*?_([A-Z]{2}\d{3}).*', text).group(1)
  

def parse_rdf_gene_id(text):
  """
  An RDF id contains a core gene id plus a decimal to indicate the variant
  """
  
  return re.match(r'(RDF\d+).*', text).group(1)
  
  
def parse_rdf_gene_variant_id(text):
  return int(re.match(r'.*(\d+)$', text).group(1))


def create_variant(id, chr, start, end):
  return id + "#" + chr + ":" + str(start) + "-" + str(end)


def parse_location_from_variant(location):
  matcher = re.match(r'.+(chr.+):(\d+)-(\d+)', location)
    
  chr = matcher.group(1)
  start = int(matcher.group(2))
  end = int(matcher.group(3))
  
  location = lib.genomic.Location(chr, start, end)
  
  return location
  

def parse_id_from_variant(variant_id):
  matcher = re.match(r'^([^#]+).*', variant_id)
    
  return matcher.group(1)
