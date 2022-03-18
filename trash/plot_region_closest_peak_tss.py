import matplotlib
matplotlib.use('agg')

import sys

import lib.plot

file = sys.argv[1]
id1 = sys.argv[2]
id2 = sys.argv[3]
color = sys.argv[4]

lib.plot.plot_tss(file, color, "Distance from closest " + id2 + " region (kb)", id1 + " region count")
