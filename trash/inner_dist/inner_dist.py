#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 12:24:32 2018

@author: antony
"""

import sys
import subprocess
import collections

#SAMTOOLS='/ifs/scratch/cancer/Lab_RDF/abh2138/tools/samtools-1.8/bin/samtools'
SAMTOOLS='/ifs/scratch/cancer/Lab_RDF/abh2138/tools/samtools-0.1.19/samtools'

read_map = collections.defaultdict(lambda: (str, str, int, int))

bam = sys.argv[1]

f = open('inner_dist.txt', 'w')

f.write('inner_dist\n')

cmd = [SAMTOOLS, 'view', '-f', '3', bam]

print(cmd)

stdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout

c = 1

pc = 0

for l in stdout:
  l = l.decode("utf-8")
    
  tokens = l.strip().split("\t")
  
  name = tokens[0]
  chr = tokens[2]
  start = int(tokens[3])
  seq = tokens[9]
  seqn = len(seq)
  
  d1 = (name, chr, start, seqn)
  
  #print(d1)
  
  if name in read_map:
    # found a pair so check the distance
    
    # since ordered by start, d2 comes first
    d2 = read_map[name]
    
    # end of first read
    s1 = d2[2] + d2[3]
    
    # start of second read
    s2 = d1[2]
    
    inner = s2 - s1
    
    #print(name, inner)
    
    del read_map[name]
    
    f.write('{}\n'.format(str(inner)))
    
    #if abs(inner) > 50000:
    #  print(name, s1, s2)
    
    #if pc == 10:
    #  break
    
    pc += 1
  else:
    read_map[name] = d1
    
  if c % 1000000 == 0:
    print('Processed', str(c), 'reads...', str(len(read_map)))
  
  c += 1
  
  
  
  
f.close()
stdout.close()
