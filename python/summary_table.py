# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 17:31:13 2014

@author: Antony Holmes
"""

# summary table of mappings and peaks at specific thresholds

import sys
import collections
import re

input_types = []
thresholds = []

for arg in sys.argv[1:]:
  sys.stderr.write(arg + "\n")
  
  input_type, p = arg.split(":", 1)
  
  input_types.append(input_type)
  thresholds.append(p)
  

# open alignment_qualities file

alignments = collections.defaultdict(lambda: collections.defaultdict(str))

f = open("../../../reads/alignment_qualities.txt", "r")

# skip header
f.readline()
  
for line in f:
  ls = line.strip()
    
  if len(ls) == 0:
    continue
    
  tokens = ls.split("\t")
  
  id = tokens[0]

  sys.stderr.write(id + "\n")
  
  matcher = re.search(r'(RK\d+)', id)
  
  seqid = matcher.group(1)
	 
  reads = tokens[1]
  mapped = tokens[2]
  p = tokens[3]
  
  alignments[seqid]["name"] = id
  alignments[seqid]["reads"] = reads
  alignments[seqid]["mapped"] = mapped
  alignments[seqid]["p"] = p

f.close()

peaks = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(str)))

f = open("peak_stats.txt", "r")

# skip header
f.readline()

for line in f:
  ls = line.strip()
    
  if len(ls) == 0:
    continue
    
  tokens = ls.split("\t")
  
  id = tokens[0]
  count = tokens[1]
  
  matcher = re.search(r'(RK\d+)', id)
  seqid = matcher.group(1)

  matcher = re.search(r'p(\d+)', id)
  p = matcher.group(1)
  
  input_type = ""

  
  
  if re.match(r'.*Ig.*', id):
    input_type = "ig"
  else:
    input_type = "input"
  
  sys.stderr.write(input_type + " " + id + "\n")
  
  peaks[seqid][input_type][p] = count

f.close()

f = open("summary_table.txt", "w")

f.write("\t".join(["SEQ ID", "# of reads", "Mapped reads (<= 2mm)", "% mapping"]))

for i in range(0, len(input_types)):
    input_type = input_types[i]
    threshold = thresholds[i]
    
    f.write("\t# of peaks (vs {0}) [p <10E-{1}]".format(input_type, threshold));

f.write("\n")

for seqid in sorted(alignments):
  reads = alignments[seqid]["reads"]
  mapped = alignments[seqid]["mapped"]
  p = float(alignments[seqid]["p"])
  
  f.write("\t".join([seqid, reads, mapped, "{0:.2f}".format(p)]))
  
  for i in range(0, len(input_types)):
    input_type = input_types[i].lower()
    threshold = thresholds[i]

    output = "n/a"
    
    if seqid in peaks:
      if input_type in peaks[seqid]:
        if threshold in peaks[seqid][input_type]:
          output = peaks[seqid][input_type][threshold]
            
    f.write("\t" + output);
    
    #sys.stderr.write(input_types[i] + "\n")
  
  f.write("\n")

f.close()