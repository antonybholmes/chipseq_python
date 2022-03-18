# -*- coding: utf-8 -*-
"""
Find the common overlaps between a list of files reported as a common region

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re
import collections

import pychipseq.peaks
import pychipseq.genomic
import pychipseq.text

# special mode to read multiple files
if sys.argv[1] == "files.txt":
  files = []
  
  f = open(sys.argv[1], 'r')
  f.readline()
  
  for line in f:
    files.append(line.strip())
else:
  files = sys.argv[1:len(sys.argv)] #sorted(sys.argv[1:len(sys.argv)])
  
ids = []

pvalues = collections.defaultdict(float)
scores = collections.defaultdict(float)


for file in files:
  print(f'file: {file}', file=sys.stderr)
  
  matcher = re.match(r'.*TF_targets.(.+)_vs.*', file)
  
  if matcher != None:  
    id = matcher.group(1)
  else:
    matcher = re.match(r'.*Peaks_(.+?)(_vs)?\..*', file)
    
    if matcher != None:  
      id = matcher.group(1)
    else:
      print(f'Wrong type of file: {file}', file=sys.stderr)
      continue
    
  ids.append(id)
  
  print(f'ids {id}', file=sys.stderr)

  # now make a list of locations and best p-values  

  f = open(file, 'r')
  
  p_col = -1 #3
  score_col = -1 #4
  
  # Adjust colums to look it for peak files
  if "Peaks" in file:
    tokens = f.readline().strip().split("\t")
    
    p_col = pychipseq.text.find_index(tokens, "P-value")
    score_col = pychipseq.text.find_index(tokens, "Score")
    
  for line in f:
    tokens = line.strip().split("\t")
    
    #sys.stderr.write(tokens[0] + "\n")
    
    if pychipseq.genomic.is_location(tokens[0]):
      location = pychipseq.genomic.parse_location(tokens[0])
    else:
      if pychipseq.genomic.is_chr(tokens[0]):
        location = pychipseq.genomic.Location(tokens[0], int(tokens[1]), int(tokens[2]))
      else:
        print(f'Invalid line: {line}', file=sys.stderr)
        
        continue
    
    if file.endswith(".bed"):
      score = 0
      p = 0
    else:
      # account for inf values
      if re.match(r'.*inf.*', tokens[p_col]):    
        p = -1000
      else:
        if p_col != -1:
          p = float(tokens[p_col])
        else:
          p = 0
      
      if score_col != -1:
        score = float(tokens[score_col])
      else:
        score = 0
    
    location = f'{id}={location.chr}:{location.start}-{location.end}'
    
    pvalues[location] = p
    scores[location] = score
    
    #sys.stderr.write(location + " " + str(p) + "\n")
    
  f.close()

  
location_core_map = pychipseq.peaks.overlapping_peaks(files, ids)

# keep the ids sorted
#ids = sorted(ids)

sys.stdout.write("Genomic Location (hg19)\tPeak / Overlap Width\tBest P-value (ChIPseeqer)\tBest Score (ChIPseeqer)\tNumber Of Overlapping Peaks\t" + "\t".join(["Sample " + s for s in ids]) + "\n")

for core_location in sorted(location_core_map):
  location = pychipseq.genomic.parse_location(core_location)
  
  overlap = location.end - location.start + 1
  
  c = 0
  
  for id in location_core_map[core_location]:
    #for location in location_core_map[core_location][id]:
    c += 1
      
  p = 0

  score = 0  
  
  for id in ids:
    if id in location_core_map[core_location]:
      #for location in location_core_map[core_location][id]:
      location = location_core_map[core_location][id]
      
      if pvalues[location] < p:
        p = pvalues[location]

      if scores[location] > score:
        score = scores[location]
  
  sys.stdout.write("\t".join([core_location, str(overlap), str(p), str(score), str(c)]));

  # foreach RK id
  for id in ids:
    if id in location_core_map[core_location]:
      peaks = location_core_map[core_location][id] #";".join(sorted(location_core_map[core_location][id]))
      
      # get rid of internal id
      #lib.peaks = re.sub(r'^.+=', "", lib.peaks)
      
      sys.stdout.write("\t" + peaks)
    else:
      sys.stdout.write("\tn/a")
  
  sys.stdout.write("\n")
