"""
Add histone mark overlaps to classify promoters
"""

import sys
import re
import time
import lib.text

def format(name, active, poised, silenced, other):
  """
  Output counts in a consistent manner
  """
  sys.stdout.write(name + "-active\t" + str(active) + "\n")
  sys.stdout.write(name + "-poised\t" + str(poised) + "\n")
  sys.stdout.write(name + "-silenced\t" + str(silenced) + "\n")
  sys.stdout.write(name + "-other\t" + str(other) + "\n")

#
# Begin
#

file = sys.argv[1]
p5 = int(sys.argv[2])
p3 = int(sys.argv[3])

sys.stdout.write(file + "\t" + str(p5) + "\t" + str(p3) + "\n")

out = re.sub(r'\.txt', '_classified_' + time.strftime('%Y%m%d') + '.txt', file)

f = open(file, 'r')
f2 = open(out, 'w')

header = f.readline().strip()

f2.write('{}\n'.format("\t".join([header, "Classification"])))
#f2.write('{}\n'.format("\t".join([header, "Classification", "Strict Classification"])))


tokens = header.split("\t")

tss_col = lib.text.find_index(tokens, "TSS Closest Distance")

if tss_col == -1:
  tss_col = lib.text.find_index(tokens, "Closest Region TSS Distance")
  
relative_col = lib.text.find_index(tokens, "Relative To Closest Gene")

if relative_col == -1:
  relative_col = lib.text.find_index(tokens, "Closest Region Relative To Gene")

#sys.stderr.write(tokens[relative_col] + "\n")

h3k4me3_col = lib.text.find_index(tokens, "H3K4me3")
h3k27me3_col = lib.text.find_index(tokens, "H3K27me3")
h3k4me1_col = lib.text.find_index(tokens, "H3K4me1")
h3k27ac_col = lib.text.find_index(tokens, "H3K27Ac")

promoter_active = 0
# Jiyuans less strict definition of promoter active
promoter_active_other = 0
promoter_poised = 0
promoter_silenced = 0
promoter_other = 0
intronic_active = 0
intronic_poised = 0
intronic_silenced = 0
intronic_other = 0
exonic_active = 0
exonic_poised = 0
exonic_silenced = 0
exonic_other = 0
intergenic_active = 0
intergenic_poised = 0
intergenic_silenced = 0
intergenic_other = 0


for line in f:
  line = line.strip()
  
  tokens = line.split("\t")
  
  relative = tokens[relative_col]
  
  tss = int(tokens[tss_col])
  
  h3k4me3 = tokens[h3k4me3_col]
  h3k27me3 = tokens[h3k27me3_col]
  h3k4me1 = tokens[h3k4me1_col]
  h3k27ac = tokens[h3k27ac_col]
  
  classification = "n/a"
  strict_classification = "n/a"
  
  if tss >= p5 and tss <= p3 and ("intergenic" not in relative):
    #if h3k4me3 != "n/a" and h3k27me3 == "n/a":
      #promoter_active += 1
      #classification = "promoter_active"
    #elif h3k4me3 != "n/a" and h3k27me3 != "n/a":
      #promoter_poised += 1
      #classification = "promoter_poised"
    #elif h3k4me3 == "n/a":
      #promoter_silenced += 1
      #classification = "promoter_silenced"
    #else:
      #promoter_other += 1
    
    # Jiyuan
    
    if h3k4me3 != "n/a":
      # active
      if h3k27me3 == "n/a":
        if h3k27ac != "n/a":
          # A stricter definition of promoter than Jiyuan
          promoter_active += 1
          classification = "promoter_active"
        else:
          promoter_other += 1
          classification = "promoter_other"
      elif h3k27me3 != "n/a":
        promoter_poised += 1
        classification = "promoter_poised"
      else:
        promoter_other += 1
    else:
      promoter_silenced += 1
      classification = "promoter_silenced"
  else:
    pass
  
  if tss > p3:
    if "intronic" in relative:
      if h3k4me3 == "n/a" and h3k4me1 != "n/a" and h3k27ac != "n/a":
        intronic_active += 1
        classification = "intronic_active"
      elif h3k4me3 == "n/a" and h3k4me1 != "n/a" and h3k27ac == "n/a":
        intronic_poised += 1
        classification = "intronic_poised"
      elif h3k4me3 == "n/a" and h3k4me1 == "n/a" and h3k27ac == "n/a":
        intronic_silenced += 1
        classification = "intronic_silenced"
      else:
        intronic_other += 1
    elif "exonic" in relative:
      if h3k4me3 == "n/a" and h3k4me1 != "n/a" and h3k27ac != "n/a":
        exonic_active += 1
        classification = "exonic_active"
      elif h3k4me3 == "n/a" and h3k4me1 != "n/a" and h3k27ac == "n/a":
        exonic_poised += 1
        classification = "exonic_poised"
      elif h3k4me3 == "n/a" and h3k4me1 == "n/a" and h3k27ac == "n/a":
        exonic_silenced += 1
        classification = "exonic_silenced"
      else:
        exonic_other += 1
        classification = "exonic_other"
    else:
      pass
  else:
    pass
    
  if ((tss < p5 or tss > p3) or (tss >= p5 and tss <= p3)) and "intergenic" in relative:
    if h3k4me3 == "n/a" and h3k4me1 != "n/a" and h3k27ac != "n/a":
      intergenic_active += 1
      classification = "intergenic_active"
    elif h3k4me3 == "n/a" and h3k4me1 != "n/a" and h3k27ac == "n/a":
      intergenic_poised += 1
      classification = "intergenic_poised"
    elif h3k4me3 == "n/a" and h3k4me1 == "n/a" and h3k27ac == "n/a":
      intergenic_silenced += 1
      classification = "intergenic_silenced"
    else:
      intergenic_other += 1
      classification = "intergenic_other"
  else:
    pass
    
  f2.write('{}\t{}\n'.format(line, classification))
  #f2.write('{}\n'.format("\t".join([classification, strict_classification])))
    
f.close()

format("promoter", promoter_active, promoter_poised, promoter_silenced, promoter_other)
format("exonic", exonic_active, exonic_poised, exonic_silenced, exonic_other)
format("intronic", intronic_active, intronic_poised, intronic_silenced, intronic_other)
format("intergenic", intergenic_active, intergenic_poised, intergenic_silenced, intergenic_other)
