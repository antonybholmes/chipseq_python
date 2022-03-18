import matplotlib
matplotlib.use('agg')

import sys
import lib.plot

file = sys.argv[1]
id = sys.argv[2]
color = sys.argv[3] #'#2c5aa0'

lib.plot.plot_log10_tss(file, color, "Absolute distance from TSS (bp)", id + " count")
