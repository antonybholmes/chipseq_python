import sys
import lib.plot

file = sys.argv[1]
id1 = sys.argv[2]
id2 = sys.argv[3]
color = sys.argv[4] #'#2c5aa0'

sys.stderr.write(file + "\n")

lib.plot.plot_log10_tss(file, color, "Absolute distance from closest " + id2 + " region (bp)", id1 + " region count")
