# -*- coding: utf-8 -*-
"""
Rename files to conform with Katia's specification

Created on Wed Jul 30 15:41:33 2014

@author: Antony Holmes
"""

import sys
import re

file = sys.argv[1]

f = open(file, 'r')

header = f.readline().strip()

matcher = re.match(r'.*(RK\d+)_([^_]+)_([^_]+).*', header)

id = matcher.group(1)
type = matcher.group(2)
tf = matcher.group(3)
input = matcher.group(4)
p = matcher.group(5)

new_id = type + "_" + tf + "_" + id + "_vs_" + input + "_" + p

header = re.sub(r'name=\".+?\"', "name=\"" + new_id + "\"", header)
header = re.sub(r'description=\".+?\"', "description=\"" + new_id + "\"", header)

new_file = new_id + ".bed"

sys.stderr.write(file + " to " + new_file + "\n")

fout = open(new_file, 'w')

fout.write(header + "\n")

for line in f:
  line = line.strip()
  
  if len(line) == 0:
    continue

  fout.write(line + "\n")


f.close()
fout.close()