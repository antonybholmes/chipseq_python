# -*- coding: utf-8 -*-
"""
Tally up the different peak types

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import collections

import lib.text


def peak_type_dist(file):
  
  counts = collections.Counter()
  
  f = open(file, 'r')

  # skip header
  header = f.readline().strip().split("\t")
  
  
  region_column = lib.text.find_index(header, "Peak Relative To Closest Gene")

  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    types = tokens[region_column].split(";")
    
    for type in types:
      counts[type] += 1
      
#      sub_types = type.split(",")
#      
#      for sub_type in sub_types:
#        if sub_type == "n/a":
#          continue
#  
#        counts[sub_type] += 1

    counts["peaks"] += 1        

  f.close()
  
  #sys.stdout.write("peaks\t" + str(counts["peaks"]) + "\n")

  fracs = []  
  labels = []
  
  for type in sorted(counts):
    if type == "peaks":
      continue
    
    sys.stdout.write(type + "\t" + str(counts[type]) + "\n")
    
    fracs.append(counts[type])
    labels.append(type)
  
  #path = '/usr/share/matplotlib/mpl-data/fonts/ttf/Arial.ttf'
  #prop = matplotlib.font_manager.FontProperties(fname=path)
  #matplotlib.rcParams['font.family'] = prop.get_name()

  #sys.stderr.write(prop.get_name() + "\n")

  #pylab.rc('font', family='Helvetica')
  #pylab.rc('text', usetex='false')
  #matplotlib.font_manager.findSystemFonts(fontpaths=['/usr/share/matplotlib/mpl-data/fonts/ttf'], fontext='ttf')
  #matplotlib.rcParams['font.family'] = 'Arial'
    
  #f = pylab.figure(num=1, figsize=(12, 12), dpi=300)
  
  #pylab.pie(fracs, labels=labels, shadow=False, startangle=90)
  
  #f.savefig("temp.pdf", dpi=300)
  
  
    
file = sys.argv[1]

peak_type_dist(file)
