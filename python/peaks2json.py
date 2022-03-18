import sys
import collections
import re

file = sys.argv[1]

f = open(file, 'r')

peaks = []

for line in f:
  tokens = line.strip().split("\t")
  
  chr = tokens[0]
  start = int(tokens[1])
  end = int(tokens[2])
  
  #json = "{\"c\":\"" + chr + "\", \"s\":" + str(start) + ", \"e\":" + str(end) + "}";
  
  json = "\"" + chr + ":" + str(start) + "-" + str(end) + "\"";
  
  peaks.append(json)

f.close()


file = re.sub(r'TF_targets_','', re.sub(r'\.txt', '.json', file))

name = re.sub(r'\..+', '', re.sub(r'^.+\/', '', file))

file = name + ".json"

sys.stderr.write(file + " " + name + "\n")

f = open(file, 'w')

#f.write("{\"id\":\"" + name + "\",\"n\":" + str(len(peaks)) + ",\"l\":[")
f.write("[")


for i in range(0, len(peaks)):
  f.write(peaks[i])
  
  if i < len(peaks) - 1:
    f.write(",")
   
  # For testing
  #if i == 2:
  #  break
    
#f.write("]}")
f.write("]")

f.close()


