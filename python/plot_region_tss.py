import matplotlib
matplotlib.use('agg')
import sys

import pychipseq.plot

file = sys.argv[1]
id = sys.argv[2]
color = sys.argv[3]

pychipseq.plot.plot_tss(file, color, "Distance from TSS (kb)", id + " count")
